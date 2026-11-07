"""Amplicon-level screening: which organisms carry sequence similar to a marker amplicon?

Two complementary checks on the reference amplicon of every selected locus:

* ``background_identity``: minimap2 of the amplicons against each background genome;
  reports the best identity per genome.
* ``remote_blast``: megablast against NCBI nt through the NCBI web service (Biopython
  ``qblast``); reports every hit with a guessed organism name and classifies it as
  target, known non-target, or unknown taxon.

An organism is a *risk* when its identity to the amplicon is at least the tolerated
identity of that amplicon, i.e. when it is at least as similar to the target as the
closest designed non-target is. Such taxa should be added as non-targets or background.
"""

from __future__ import annotations

import logging
import re
import subprocess
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import pandas as pd

from .loci import parse_paf

log = logging.getLogger(__name__)

STATUS_TARGET = "target"
STATUS_NONTARGET = "known non-target"
STATUS_UNKNOWN = "unknown taxon"


@dataclass
class BlastHit:
    query: str
    accession: str
    title: str
    organism: str
    pct_identity: float
    query_coverage: float
    evalue: float
    bitscore: float


_ORGANISM_RE = re.compile(r"^(?:PREDICTED:\s*|TPA(?:_asm)?:\s*|MAG:\s*|UNVERIFIED(?:_ORG)?:\s*)?(?:Candidatus\s+)?\[?([A-Z][a-z]+)\]?\s+([a-z][a-z\-]+)")


def guess_organism(title: str) -> str:
    """'Yersinia pestis strain X chromosome' -> 'Yersinia pestis'. Falls back to the first word."""
    m = _ORGANISM_RE.match(title.strip())
    if m:
        genus, species = m.group(1), m.group(2)
        if species in {"sp", "sp.", "strain", "isolate", "genome", "chromosome", "plasmid", "clone"}:
            return f"{genus} sp."
        return f"{genus} {species}"
    return title.split()[0] if title.split() else "unknown"


def parse_blast_xml(handle, min_identity: float = 0.0) -> Dict[str, List[BlastHit]]:
    """Parse classic BLAST XML; keep the best HSP per hit."""
    from Bio.Blast import NCBIXML

    hits_by_query: Dict[str, List[BlastHit]] = {}
    for record in NCBIXML.parse(handle):
        qid = (record.query or "").split()[0] if record.query else record.query_id
        qlen = record.query_length or 1
        out: List[BlastHit] = []
        for aln in record.alignments:
            best = max(aln.hsps, key=lambda h: h.bits)
            pct = 100.0 * best.identities / best.align_length if best.align_length else 0.0
            cov = 100.0 * (abs(best.query_end - best.query_start) + 1) / qlen
            if pct < min_identity:
                continue
            title = aln.hit_def or aln.title or ""
            out.append(BlastHit(qid, aln.accession or "", title, guess_organism(title), pct, cov, float(best.expect), float(best.bits)))
        hits_by_query[qid] = sorted(out, key=lambda h: -h.bitscore)
    return hits_by_query


def remote_blast(records: Sequence[Tuple[str, str]], db: str = "nt", megablast: bool = True,
                 hitlist_size: int = 100, entrez_query: Optional[str] = None, expect: float = 1e-20) -> str:
    """Submit all records in one request to NCBI BLAST and return the XML text."""
    from Bio.Blast import NCBIWWW

    fasta = "".join(f">{rid}\n{seq}\n" for rid, seq in records)
    kwargs = dict(hitlist_size=hitlist_size, expect=expect, format_type="XML")
    if megablast:
        kwargs["megablast"] = True
    if entrez_query:
        kwargs["entrez_query"] = entrez_query
    log.info("Submitting %d amplicon(s) to NCBI BLAST (%s, %s); this may take several minutes",
             len(records), "megablast" if megablast else "blastn", db)
    handle = NCBIWWW.qblast("blastn", db, fasta, **kwargs)
    return handle.read()


def _norm(name: str) -> str:
    return name.replace("_", " ").strip().lower()


def classify_organism(organism: str, target_species: Iterable[str], nontarget_species: Iterable[str]) -> str:
    org = _norm(organism)
    for sp in target_species:
        if org.startswith(_norm(sp)):
            return STATUS_TARGET
    for sp in nontarget_species:
        if org.startswith(_norm(sp)):
            return STATUS_NONTARGET
    return STATUS_UNKNOWN


def summarize_blast(
    hits_by_query: Mapping[str, Sequence[BlastHit]],
    target_species: Iterable[str],
    nontarget_species: Iterable[str],
    tolerated_identity: Mapping[str, float],
) -> pd.DataFrame:
    rows = []
    for q, hits in hits_by_query.items():
        tol = tolerated_identity.get(q, 100.0)
        for h in hits:
            status = classify_organism(h.organism, target_species, nontarget_species)
            rows.append({
                "Locus": q, "Organism": h.organism, "Accession": h.accession, "Pct_Identity": round(h.pct_identity, 2),
                "Query_Coverage": round(h.query_coverage, 1), "E_value": h.evalue, "Bitscore": h.bitscore,
                "Status": status, "Tolerated_Identity": round(tol, 2),
                "Risk": status == STATUS_UNKNOWN and h.pct_identity >= tol and h.query_coverage >= 50.0,
                "Title": h.title,
            })
    cols = ["Locus", "Organism", "Accession", "Pct_Identity", "Query_Coverage", "E_value", "Bitscore", "Status",
            "Tolerated_Identity", "Risk", "Title"]
    return pd.DataFrame(rows, columns=cols)


def tolerated_identity(amplicon_len: int, min_sites_to_nearest_nontarget: Optional[int], fallback: float) -> float:
    """Identity above which an unknown organism is closer to the target than its nearest designed non-target."""
    if not min_sites_to_nearest_nontarget or amplicon_len <= 0:
        return fallback
    return max(fallback, 100.0 * (1.0 - min_sites_to_nearest_nontarget / amplicon_len))


def background_identity(
    amplicons_fa: Path,
    background_files: Sequence[Path],
    work_dir: Path,
    threads: int,
    tolerated: Mapping[str, float],
    preset: str = "asm20",
) -> pd.DataFrame:
    """Best minimap2 identity of each amplicon in each background genome."""
    work_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for path in background_files:
        paf = work_dir / f"{path.stem}.paf"
        cmd = ["minimap2", "-c", "-x", preset, "--secondary=no", "-t", str(max(1, threads)), str(path), str(amplicons_fa)]
        with open(paf, "w") as out:
            subprocess.run(cmd, stdout=out, stderr=subprocess.DEVNULL, check=True)
        for q, hits in parse_paf(paf).items():
            best = max(hits, key=lambda h: h.nmatch)
            ident = 100.0 * best.identity
            rows.append({"Locus": q, "Background_Genome": path.stem, "Contig": best.tname, "Start": best.tstart,
                         "End": best.tend, "Pct_Identity": round(ident, 2), "Query_Coverage": round(100 * best.query_coverage, 1),
                         "Tolerated_Identity": round(tolerated.get(q, 100.0), 2),
                         "Risk": ident >= tolerated.get(q, 100.0) and best.query_coverage >= 0.5})
        paf.unlink(missing_ok=True)
    cols = ["Locus", "Background_Genome", "Contig", "Start", "End", "Pct_Identity", "Query_Coverage", "Tolerated_Identity", "Risk"]
    return pd.DataFrame(rows, columns=cols)
