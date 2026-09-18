"""Two-hospital clinical handoff. Preserve original global patient numbers."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cohorts import main

if __name__ == '__main__':
    raise SystemExit(main('two_hospital'))
