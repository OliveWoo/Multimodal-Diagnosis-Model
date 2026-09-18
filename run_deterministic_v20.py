"""Canonical offline entry from standardized patient JSON to v20 pathogen results."""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
UPSTREAM = ROOT / "upstream"
if str(UPSTREAM) not in sys.path:
    sys.path.insert(0, str(UPSTREAM))

from tools.run_deterministic_v20_pipeline import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
