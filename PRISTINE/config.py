#!/usr/bin/env python3
"""
Configuration classes for PRISTINE pipeline.

This module contains all configuration-related dataclasses and the ConfigLoader
for parsing YAML configuration files.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, Dict
import yaml
import os
import sys


@dataclass
class InputPaths:
    """Container for input directory paths."""
    raw_dir: Optional[str] = None
    prokka_dir: Optional[str] = None
    panaroo_dir: Optional[str] = None


@dataclass
class Primer3Params:
    """Container for Primer3 global and design parameters."""
    global_params: Dict[str, object] = field(default_factory=dict)
    design_params: Dict[str, object] = field(default_factory=dict)


@dataclass
class SNPPrimerDesignParams:
    """Configuration for SNP-aware primer design strategy."""
    snp_window_size: int
    snp_top_n: int
    min_snps: int


@dataclass
class ValidationConfig:
    """Configuration for post-design validation using pBLAT."""
    perform: str
    database: str
    pblat_min_identity: float
    match_median_filter_tolerance: int


@dataclass
class Config:
    """
    Main configuration class for PRISTINE pipeline.

    Validates input_type and ensures required paths are set based on the
    selected input type (raw/prokka/panaroo).
    """
    input_type: str
    input_paths: InputPaths
    output_dir: str
    max_cores: int
    aligner: str
    snp_avg_prop_threshold: float
    primer3_config_file: Optional[str] = None
    primer3: Primer3Params = field(default_factory=Primer3Params)
    snp_primer_design: SNPPrimerDesignParams = field(default_factory=SNPPrimerDesignParams)
    validation: ValidationConfig = field(default_factory=ValidationConfig)

    def __post_init__(self):
        """Validate configuration after initialization."""
        allowed_inputs = {"raw", "prokka", "panaroo"}
        if self.input_type not in allowed_inputs:
            raise ValueError(f"`input_type` must be one of {allowed_inputs}, got: {self.input_type}")

        active_dir = {
            "raw": self.input_paths.raw_dir,
            "prokka": self.input_paths.prokka_dir,
            "panaroo": self.input_paths.panaroo_dir,
        }[self.input_type]
        if not active_dir:
            raise ValueError(f"{self.input_type}_dir must be set in `input_paths` for input_type = {self.input_type}")

        if self.primer3_config_file and self.primer3.global_params:
            print("Warning: primer3_config_file is set. Inline primer3 parameters will be ignored.")


class ConfigLoader:
    """Utility class for loading YAML configuration files."""

    @staticmethod
    def load(path: str) -> Config:
        """
        Load and parse a YAML configuration file.

        Args:
            path: Path to the YAML configuration file

        Returns:
            Config object with all nested dataclasses instantiated

        Raises:
            FileNotFoundError: If the configuration file is not found
            SystemExit: If the configuration file is not found (with helpful error message)
        """
        if not os.path.exists(path):
            print(f"Error: Configuration file not found: {path}", file=sys.stderr)
            print(f"Current working directory: {os.getcwd()}", file=sys.stderr)
            print(f"\nPlease ensure 'config.yaml' is in the same directory where you run the container.", file=sys.stderr)
            print(f"Or specify a custom path using: --config /path/to/config.yaml", file=sys.stderr)
            sys.exit(1)

        with open(path, "r") as f:
            raw = yaml.safe_load(f)

        # Parse nested dataclasses manually
        raw["input_paths"] = InputPaths(**raw.get("input_paths", {}))
        raw["primer3"] = Primer3Params(**raw.get("primer3", {}))
        raw["snp_primer_design"] = SNPPrimerDesignParams(**raw.get("snp_primer_design", {}))
        raw["validation"] = ValidationConfig(**raw.get("validation", {}))

        return Config(**raw)
