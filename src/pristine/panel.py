"""Exports for PRISCREEN: the amplicon panel and the diagnostic-site tables.

Panel (format: PRISCREEN docs/priscreen_panel_format.md)::

    >locus|species=Genus_species|genome=Genome_id|role=target
    ACGT...

One record per distinct amplicon sequence of a (locus, species), in the orientation of
the reference amplicon, upper case. role is target, non-target or background (a
background genome that amplifies with a product distinguishable by sequence).

Diagnostic sites: per locus, every position at which the species differ, on the
reference amplicon. Substitutions are one row per position; an indel event (a run of
deleted reference bases, or bases inserted after one reference position, with the same
species pattern) is one row. Homopolymer context is reported, because 1-nt indels in
long homopolymers are the typical nanopore error: such an indel (homopolymer >=
HOMOPOLYMER_WEAK) is marked weak=yes and is not counted by PRISTINE as an informative site.
"""

from __future__ import annotations

import csv
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple

from .scoring import HOMOPOLYMER_WEAK
from .seqio import revcomp, write_fasta

ROLES = ("target", "non-target", "background")

# (genome_id, species, role, sequence)
PanelRow = Tuple[str, str, str, str]


def _aligner():
    from Bio.Align import PairwiseAligner

    return PairwiseAligner(mode="global", match_score=1, mismatch_score=-2, open_gap_score=-4, extend_gap_score=-1)


def orient(reference: str, seq: str) -> str:
    """`seq` in the orientation that aligns better to `reference` (upper case)."""
    seq = seq.upper()
    rc = revcomp(seq)
    if seq == reference or rc == reference:
        return reference
    al = _aligner()
    return seq if al.score(reference, seq) >= al.score(reference, rc) else rc


def panel_records(locus: str, rows: Sequence[PanelRow]) -> List[Tuple[str, str]]:
    """FASTA records: one per distinct sequence of each (locus, species), first genome kept.
    Order: targets, non-targets, background; within a role, input order."""
    seen = set()
    records = []
    for role in ROLES:
        for gid, species, r, seq in rows:
            if r != role:
                continue
            seq = seq.upper()
            if (species, seq) in seen:
                continue
            seen.add((species, seq))
            records.append((f"{locus}|species={species}|genome={gid}|role={role}", seq))
    return records


def _project(ref: str, seq: str):
    """Per reference position the aligned base or '-', and insertions {ref position after
    which they occur (-1 = before the first base): inserted bases}. `seq` must already be
    in reference orientation."""
    aln = _aligner().align(ref, seq)[0]
    states = ["-"] * len(ref)
    insertions: Dict[int, str] = {}
    ta, qa = aln.aligned
    prev_q = 0
    for (ts, te), (qs, qe) in zip(ta, qa):
        if qs > prev_q:
            insertions[ts - 1] = seq[prev_q:qs]
        for i, j in zip(range(ts, te), range(qs, qe)):
            states[i] = seq[j]
        prev_q = qe
    if prev_q < len(seq):
        insertions[len(ref) - 1] = seq[prev_q:]
    return states, insertions


def _run_length(ref: str, start: int, end: int, base: str) -> int:
    """Length of the run of `base` in `ref` that touches [start, end) (end exclusive)."""
    i, j = start, end
    while i > 0 and ref[i - 1] == base:
        i -= 1
    while j < len(ref) and ref[j] == base:
        j += 1
    return (j - i) - (end - start) + sum(1 for c in ref[start:end] if c == base)


def _homopolymer(ref: str, kind: str, start: int, end: int, alleles: Sequence[str]) -> Tuple[str, int]:
    """Homopolymer base and the longest run length among the alleles at an event, or ("", 1)."""
    if kind == "substitution":
        return "", 1
    if kind == "deletion":
        deleted = ref[start:end]
        if len(set(deleted)) != 1:
            return "", 1
        base = deleted[0]
        return base, _run_length(ref, start, end, base)
    # insertion after reference position `start`
    best = ("", 1)
    for ins in alleles:
        if not ins or len(set(ins)) != 1:
            continue
        base = ins[0]
        left = _run_length(ref, start + 1, start + 1, base)  # run touching the insertion point
        run = left + len(ins)
        if run > best[1]:
            best = (base, run)
    return best


def _species_state(values: Sequence[str]) -> str:
    vals = sorted(set(values))
    return vals[0] if len(vals) == 1 else "/".join(vals)


def diagnostic_sites(reference: str, rows: Sequence[PanelRow]) -> Tuple[List[str], List[Dict[str, object]]]:
    """Positions where the species' amplicons differ, on `reference` (all sequences in
    reference orientation). Returns (species columns, rows). A species whose sequences
    disagree at a position shows all its states, e.g. "A/G". A position is reported if two
    species share no state there; separates_target_from lists the species that share no
    state with any target genome."""
    reference = reference.upper()
    species_role: "OrderedDict[str, str]" = OrderedDict()
    for role in ROLES:
        for _, sp, r, _ in rows:
            if r == role:
                species_role.setdefault(sp, role)
    distinct: Dict[str, List[str]] = {sp: [] for sp in species_role}
    for _, sp, _, seq in rows:
        seq = seq.upper()
        if seq not in distinct[sp]:
            distinct[sp].append(seq)
    projected = {sp: [_project(reference, s) for s in seqs] for sp, seqs in distinct.items()}
    target_species = [sp for sp, r in species_role.items() if r == "target"]

    def state_at(sp: str, pos: int) -> str:
        return _species_state([st[pos] for st, _ in projected[sp]])

    def insertion_at(sp: str, pos: int) -> str:
        return _species_state([ins.get(pos, "") or "-" for _, ins in projected[sp]])

    def informative(states: Mapping[str, str]) -> bool:
        """Two species share no state here (their reads can be told apart at this site)."""
        sets = [set(v.split("/")) for v in states.values()]
        return any(not (a & b) for i, a in enumerate(sets) for b in sets[i + 1:])

    def separates(states: Mapping[str, str]) -> str:
        tset = {x for sp in target_species for x in states[sp].split("/")}
        return ",".join(sp for sp, r in species_role.items()
                        if r != "target" and not (tset & set(states[sp].split("/"))))

    out: List[Dict[str, object]] = []

    # insertions before the first base
    positions = list(range(-1, len(reference)))
    pending_del = None  # (start, pattern, states per species as list of strings)

    def flush_deletion():
        nonlocal pending_del
        if pending_del is None:
            return
        start, end, states = pending_del
        hb, hl = _homopolymer(reference, "deletion", start, end, [])
        out.append(_row(start + 1, "indel", end - start, reference[start:end], hb, hl, states, separates(states)))
        pending_del = None

    for pos in positions:
        if pos >= 0:
            states = {sp: state_at(sp, pos) for sp in species_role}
            if informative(states):
                has_gap = {sp: v == "-" for sp, v in states.items()}
                if any(has_gap.values()) and all(v in ("-", reference[pos]) for v in states.values()):
                    pattern = tuple(sorted(has_gap.items()))
                    if pending_del is not None and pending_del[1] == pos and \
                            tuple(sorted({sp: v == "-" for sp, v in pending_del[2].items()}.items())) == pattern:
                        start, _, prev = pending_del
                        merged = {sp: (prev[sp] + reference[pos] if states[sp] != "-" else "-") for sp in states}
                        pending_del = (start, pos + 1, merged)
                    else:
                        flush_deletion()
                        pending_del = (pos, pos + 1, {sp: (reference[pos] if v != "-" else "-") for sp, v in states.items()})
                    continue
                flush_deletion()
                kind = "indel" if any(v == "-" or "-" in v.split("/") for v in states.values()) else "substitution"
                hb, hl = _homopolymer(reference, "deletion", pos, pos + 1, []) if kind == "indel" else ("", 1)
                out.append(_row(pos + 1, kind, 1, reference[pos], hb, hl, states, separates(states)))
            else:
                flush_deletion()
        ins_states = {sp: insertion_at(sp, pos) for sp in species_role}
        if any(v != "-" for v in ins_states.values()) and informative(ins_states):
            flush_deletion()
            alleles = [a for v in ins_states.values() for a in v.split("/") if a != "-"]
            hb, hl = _homopolymer(reference, "insertion", pos, pos, alleles)
            length = max(len(a) for a in alleles)
            # position: the reference base after which the insertion lies (0 = before base 1)
            out.append(_row(pos + 1, "indel", length, "-", hb, hl, ins_states, separates(ins_states)))
    flush_deletion()
    return list(species_role), out


def _row(position, kind, length, ref, hb, hl, states, sep) -> Dict[str, object]:
    weak = kind == "indel" and length == 1 and hl >= HOMOPOLYMER_WEAK
    row: Dict[str, object] = {"position": position, "type": kind, "length": length, "ref": ref,
                              "homopolymer_base": hb, "homopolymer_length": hl, "weak": "yes" if weak else "no",
                              "separates_target_from": sep}
    row.update(states)
    return row


def write_diagnostic_sites(path: Path, reference: str, rows: Sequence[PanelRow]) -> int:
    """Write the diagnostic-site TSV of one locus; returns the number of sites."""
    species, sites = diagnostic_sites(reference, rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["position", "type", "length", "ref", "homopolymer_base", "homopolymer_length", "weak",
            "separates_target_from"] + species
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for s in sites:
            w.writerow(s)
    return len(sites)


def write_panel(path: Path, loci: Sequence[Tuple[str, Sequence[PanelRow]]]) -> None:
    records = [rec for locus, rows in loci for rec in panel_records(locus, rows)]
    write_fasta(path, records)


def merge_panels(group_panels: Mapping[str, Path], out_path: Path) -> None:
    """Merge group panels; a locus name occurring in two groups gets the group as prefix."""
    from .seqio import read_fasta

    per_group = {}
    owners: Dict[str, set] = {}
    for group, path in group_panels.items():
        if not Path(path).exists():
            continue
        recs = list(read_fasta(path))
        per_group[group] = recs
        for header, _ in recs:
            owners.setdefault(header.split("|", 1)[0], set()).add(group)
    if not per_group:
        return
    merged = []
    for group, recs in per_group.items():
        for header, seq in recs:
            locus, rest = header.split("|", 1)
            if len(owners[locus]) > 1:
                locus = f"{group}_{locus}"
            merged.append((f"{locus}|{rest}", seq))
    write_fasta(out_path, merged)
