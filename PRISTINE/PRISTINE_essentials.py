#!/usr/bin/env python3
"""Backwards-compatible entry point. The pipeline now lives in the `pristine` package
(src/pristine); run `pristine --config config.yaml` or `pixi run pristine`."""

import sys

from pristine.cli import main

if __name__ == "__main__":
    sys.exit(main())
