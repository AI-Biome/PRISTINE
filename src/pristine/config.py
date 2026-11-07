"""Configuration loading for PRISTINE.

All relative paths in the YAML file are resolved relative to the directory that
contains the config file, never relative to the current working directory or to
hard-coded locations.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

log = logging.getLogger(__name__)

# Primer3 defaults. Anything the user puts under primer3.global_params overrides these.
DEFAULT_PRIMER3_PARAMS: Dict[str, Any] = {
    "PRIMER_NUM_RETURN": 5,
    "PRIMER_MIN_SIZE": 18,
    "PRIMER_OPT_SIZE": 20,
    "PRIMER_MAX_SIZE": 25,
    "PRIMER_MIN_TM": 57.0,
    "PRIMER_OPT_TM": 60.0,
    "PRIMER_MAX_TM": 63.0,
    "PRIMER_PAIR_MAX_DIFF_TM": 3.0,
    "PRIMER_MIN_GC": 35.0,
    "PRIMER_MAX_GC": 65.0,
    "PRIMER_MAX_POLY_X": 4,
    "PRIMER_PRODUCT_SIZE_RANGE": [[500, 1800]],
    # Primers must not contain ambiguous bases: the template is N-masked at every
    # position that is not conserved, so this is what keeps primers in conserved sites.
    "PRIMER_MAX_NS_ACCEPTED": 0,
    "PRIMER_LIBERAL_BASE": 1,
}

# Legacy keys from the prokka/panaroo based version of PRISTINE.
LEGACY_KEYS = {"input_type", "input_paths", "aligner", "primer3_config_file", "validation"}


@dataclass
class LocusDiscoveryParams:
    window_size: int = 2000
    step: int = 500
    minimap2_preset: str = "asm20"
    min_query_coverage: float = 0.9
    max_divergence: float = 0.15
    max_secondary_coverage: float = 0.5
    max_candidate_loci: int = 5000
    require_in_nontargets: bool = True
    max_n_fraction: float = 0.05
    refine_with_mafft: bool = False


@dataclass
class SNPPrimerDesignParams:
    snp_window_size: int = 250
    snp_top_n: int = 5
    min_snps: int = 3
    max_loci_to_try: int = 50
    require_snps_for_all_nontargets: bool = True


@dataclass
class PrimerSelectionParams:
    conservation: str = "all"  # "all" or "targets"
    use_thermodynamics: bool = True
    annealing_temp: Optional[float] = None  # default: PRIMER_MIN_TM - 3
    binding_tm_margin: float = 5.0  # a site binds if its duplex Tm >= annealing_temp - margin
    max_mismatches: int = 6  # prefilter: mismatches tolerated over the whole primer
    three_prime_region: int = 8  # 3'-terminal bases scrutinised separately
    three_prime_max_mismatches: int = 1  # mismatches tolerated inside the 3' region
    seed_length: int = 4  # exact k-mer used to find candidate sites
    max_offtarget_product_size: int = 4000
    background_dir: Optional[Path] = None
    # Background genomes may amplify if the amplicon is sequenced: "distinguish" accepts a
    # background product that differs from every target amplicon by >= min_snps events;
    # "reject" fails a pair on any background product.
    background_products: str = "distinguish"


@dataclass
class AmpliconScreenParams:
    background_identity: bool = True  # minimap2 identity of selected amplicons in background genomes
    remote_blast: bool = True  # megablast of selected amplicons against NCBI nt (network, minutes; set false offline)
    blast_db: str = "nt"
    megablast: bool = True
    hitlist_size: int = 100
    entrez_query: Optional[str] = None  # e.g. "bacteria[organism]"
    min_identity_report: float = 90.0  # hits below this identity are not reported
    fallback_tolerated_identity: float = 98.0  # used when an amplicon has no non-target sites to derive it from


@dataclass
class Config:
    input_dir: Path
    output_dir: Path
    max_cores: int = 1
    snp_avg_prop_threshold: float = 0.01
    primer3_params: Dict[str, Any] = field(default_factory=dict)
    locus_discovery: LocusDiscoveryParams = field(default_factory=LocusDiscoveryParams)
    snp_primer_design: SNPPrimerDesignParams = field(default_factory=SNPPrimerDesignParams)
    primer_selection: PrimerSelectionParams = field(default_factory=PrimerSelectionParams)
    amplicon_screen: AmpliconScreenParams = field(default_factory=AmpliconScreenParams)
    make_plots: bool = True

    def __post_init__(self) -> None:
        if self.max_cores < 1:
            raise ValueError("max_cores must be >= 1")
        if self.primer_selection.conservation not in {"all", "targets"}:
            raise ValueError("primer_selection.conservation must be 'all' or 'targets'")
        if self.primer_selection.background_products not in {"distinguish", "reject"}:
            raise ValueError("primer_selection.background_products must be 'distinguish' or 'reject'")
        ld = self.locus_discovery
        if ld.step < 1 or ld.window_size < 1:
            raise ValueError("locus_discovery.window_size and step must be positive")
        if ld.step > ld.window_size:
            raise ValueError("locus_discovery.step must not exceed window_size")
        max_product = max(hi for _, hi in self.product_size_ranges())
        if ld.window_size < max_product:
            raise ValueError(
                f"locus_discovery.window_size ({ld.window_size}) must be at least the largest "
                f"PRIMER_PRODUCT_SIZE_RANGE upper bound ({max_product})"
            )
        sel = self.primer_selection
        if sel.three_prime_region < (sel.three_prime_max_mismatches + 1) * sel.seed_length:
            # Pigeonhole: m mismatches in an R-base region leave an exact k-mer only if R >= (m+1)*k.
            raise ValueError(
                f"primer_selection.three_prime_region ({sel.three_prime_region}) must be at least "
                f"(three_prime_max_mismatches + 1) * seed_length "
                f"({(sel.three_prime_max_mismatches + 1) * sel.seed_length}); otherwise the seed scan "
                f"can miss binding sites with {sel.three_prime_max_mismatches} mismatch(es) in the 3' region"
            )
        if self.primer3_params.get("PRIMER_MAX_NS_ACCEPTED", 0) != 0:
            log.warning(
                "PRIMER_MAX_NS_ACCEPTED is forced to 0: primers must lie in fully conserved sites"
            )
            self.primer3_params["PRIMER_MAX_NS_ACCEPTED"] = 0

    @property
    def annealing_temp(self) -> float:
        sel = self.primer_selection
        return float(sel.annealing_temp) if sel.annealing_temp is not None else float(self.primer3_params["PRIMER_MIN_TM"]) - 3.0

    @property
    def min_site_tm(self) -> Optional[float]:
        """Duplex Tm above which a genomic site counts as a binding site (None: mismatch-only mode)."""
        if not self.primer_selection.use_thermodynamics:
            return None
        return self.annealing_temp - self.primer_selection.binding_tm_margin

    def product_size_ranges(self) -> List[List[int]]:
        return [list(map(int, r)) for r in self.primer3_params["PRIMER_PRODUCT_SIZE_RANGE"]]

    @property
    def min_product_size(self) -> int:
        return min(lo for lo, _ in self.product_size_ranges())

    @property
    def max_product_size(self) -> int:
        return max(hi for _, hi in self.product_size_ranges())


def parse_product_size_range(value: Any) -> List[List[int]]:
    """Accept "500-1800", "500-1800 1900-2000", [500, 1800], [[500, 1800], ...]."""
    if isinstance(value, str):
        ranges = []
        for token in value.replace(",", " ").split():
            lo, hi = token.split("-")
            ranges.append([int(lo), int(hi)])
    elif isinstance(value, (list, tuple)):
        if len(value) == 2 and all(isinstance(v, int) for v in value):
            ranges = [[int(value[0]), int(value[1])]]
        else:
            ranges = [[int(lo), int(hi)] for lo, hi in value]
    else:
        raise ValueError(f"Cannot parse PRIMER_PRODUCT_SIZE_RANGE: {value!r}")
    for lo, hi in ranges:
        if lo <= 0 or hi < lo:
            raise ValueError(f"Invalid product size range {lo}-{hi}")
    return ranges


def _resolve(base: Path, value: Optional[str]) -> Optional[Path]:
    if value is None:
        return None
    p = Path(os.path.expanduser(str(value)))
    return p if p.is_absolute() else (base / p).resolve()


def _dataclass_from_dict(cls, data: Optional[dict], name: str):
    data = dict(data or {})
    allowed = set(cls.__dataclass_fields__)
    unknown = set(data) - allowed
    if unknown:
        raise ValueError(f"Unknown keys in `{name}`: {sorted(unknown)}; allowed: {sorted(allowed)}")
    return cls(**data)


def load_config(path: os.PathLike | str) -> Config:
    """Load and validate a config.yaml. Relative paths are resolved against its directory."""
    cfg_path = Path(path).resolve()
    with open(cfg_path) as fh:
        raw = yaml.safe_load(fh) or {}
    base = cfg_path.parent

    legacy = LEGACY_KEYS & set(raw)
    input_dir = raw.get("input_dir")
    if input_dir is None and "input_paths" in raw:
        input_dir = (raw.get("input_paths") or {}).get("raw_dir")
        log.warning("`input_paths.raw_dir` is deprecated; use `input_dir`")
    for key in sorted(legacy - {"input_paths"}):
        log.warning("Config key `%s` is no longer used (prokka/panaroo/obipcr were removed) and is ignored", key)
    if input_dir is None:
        raise ValueError("`input_dir` must be set in the config")
    if "output_dir" not in raw:
        raise ValueError("`output_dir` must be set in the config")

    primer3_params = dict(DEFAULT_PRIMER3_PARAMS)
    user_p3 = (raw.get("primer3") or {}).get("global_params") or {}
    primer3_params.update(user_p3)
    primer3_params["PRIMER_PRODUCT_SIZE_RANGE"] = parse_product_size_range(
        primer3_params["PRIMER_PRODUCT_SIZE_RANGE"]
    )
    # Keep OPT values inside the user's MIN/MAX so primer3 does not abort.
    for opt, lo, hi in (
        ("PRIMER_OPT_TM", "PRIMER_MIN_TM", "PRIMER_MAX_TM"),
        ("PRIMER_OPT_SIZE", "PRIMER_MIN_SIZE", "PRIMER_MAX_SIZE"),
    ):
        if not (primer3_params[lo] <= primer3_params[opt] <= primer3_params[hi]):
            primer3_params[opt] = type(primer3_params[opt])((primer3_params[lo] + primer3_params[hi]) / 2)
            log.warning("%s was outside [%s, %s]; set to %s", opt, lo, hi, primer3_params[opt])

    sel_raw = dict(raw.get("primer_selection") or {})
    sel_raw["background_dir"] = _resolve(base, sel_raw.get("background_dir"))

    snp_raw = dict(raw.get("snp_primer_design") or {})

    cfg = Config(
        input_dir=_resolve(base, input_dir),
        output_dir=_resolve(base, raw["output_dir"]),
        max_cores=int(raw.get("max_cores", 1)),
        snp_avg_prop_threshold=float(raw.get("snp_avg_prop_threshold", 0.01)),
        primer3_params=primer3_params,
        locus_discovery=_dataclass_from_dict(LocusDiscoveryParams, raw.get("locus_discovery"), "locus_discovery"),
        snp_primer_design=_dataclass_from_dict(SNPPrimerDesignParams, snp_raw, "snp_primer_design"),
        primer_selection=_dataclass_from_dict(PrimerSelectionParams, sel_raw, "primer_selection"),
        amplicon_screen=_dataclass_from_dict(AmpliconScreenParams, raw.get("amplicon_screen"), "amplicon_screen"),
        make_plots=bool(raw.get("make_plots", True)),
    )
    if not cfg.input_dir.is_dir():
        raise FileNotFoundError(f"input_dir does not exist: {cfg.input_dir}")
    return cfg
