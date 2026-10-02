from __future__ import annotations

from pathlib import Path

from pptx import Presentation


def collect_pptx_files(input_dir: Path) -> list[Path]:
    return sorted(
        path.resolve()
        for path in input_dir.rglob("*.pptx")
        if path.is_file() and not path.name.startswith("~$")
    )


def open_presentation(path: Path) -> Presentation:
    return Presentation(str(path))
