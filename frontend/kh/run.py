"""KH-only clinical handoff. Local patient number means K case code."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cohorts import main

if __name__ == '__main__':
    raise SystemExit(main('KH'))
