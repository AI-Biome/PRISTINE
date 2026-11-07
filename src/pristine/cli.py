"""Command line entry point: ``pristine [--config config.yaml] [--threads N] [--output DIR]``."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from . import __version__
from .config import load_config
from .pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="pristine", description="PRISTINE diagnostic primer design")
    p.add_argument("-c", "--config", default="config.yaml", help="path to config.yaml (default: ./config.yaml)")
    p.add_argument("-t", "--threads", type=int, default=None, help="override max_cores from the config")
    p.add_argument("-o", "--output", default=None, help="override output_dir from the config")
    p.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    p.add_argument("--version", action="version", version=f"pristine {__version__}")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=args.log_level, format="%(asctime)s %(levelname)s %(name)s: %(message)s",
                        datefmt="%H:%M:%S", stream=sys.stderr)
    cfg = load_config(args.config)
    if args.threads is not None:
        cfg.max_cores = max(1, args.threads)
    if args.output is not None:
        cfg.output_dir = Path(args.output).resolve()
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(cfg.output_dir / "pristine.log")
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s"))
    logging.getLogger().addHandler(fh)
    df = run_pipeline(cfg)
    print(df.to_string(index=False))
    return 1 if "Error" in df.columns and df["Error"].notna().any() else 0


if __name__ == "__main__":
    sys.exit(main())
