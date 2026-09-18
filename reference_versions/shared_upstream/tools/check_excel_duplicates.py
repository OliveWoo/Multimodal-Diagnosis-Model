"""Scan the raw_excel directory (or a user-provided path) for duplicate Excel files."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

DEFAULT_RAW_EXCEL_DIR = Path("raw_excel")
LOG_DIR = Path("outputs") / "logs"
REPORT_FILE = LOG_DIR / "duplicate_report.txt"
DEFAULT_CHUNK_SIZE = 131072

SkipInfo = Tuple[Path, str]

def iter_excel_files(base_dir: Path) -> Iterable[Path]:
    if not base_dir.exists():
        return []
    return (path for path in sorted(base_dir.rglob("*.xlsx")) if path.is_file())

def compute_md5(path: Path, chunk_size: int = DEFAULT_CHUNK_SIZE) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()

def group_by_digest(files: Iterable[Path]) -> Dict[str, List[Path]]:
    groups: Dict[str, List[Path]] = {}
    for file_path in files:
        digest = compute_md5(file_path)
        groups.setdefault(digest, []).append(file_path)
    return groups

def filter_duplicates(groups: Dict[str, List[Path]]) -> Dict[str, List[Path]]:
    return {digest: paths for digest, paths in groups.items() if len(paths) > 1}

def _safe_group_by_digest(files: Iterable[Path]) -> Tuple[Dict[str, List[Path]], List[SkipInfo]]:
    groups: Dict[str, List[Path]] = {}
    skipped: List[SkipInfo] = []
    for file_path in files:
        path = Path(file_path)
        if not path.exists():
            skipped.append((path, "檔案不存在"))
            continue
        if not path.is_file():
            skipped.append((path, "路徑不是可讀取的檔案"))
            continue
        try:
            digest = compute_md5(path)
        except OSError as exc:
            skipped.append((path, str(exc)))
            continue
        groups.setdefault(digest, []).append(path)
    return groups, skipped

def find_duplicates(files: Iterable[Path]) -> List[List[Path]]:
    """提供主流程呼叫的重複檔案檢查介面。"""
    groups, skipped = _safe_group_by_digest(files)
    if skipped:
        for path, reason in skipped:
            print(f"[duplicate-check] 無法讀取 {path}: {reason}")
    duplicates = filter_duplicates(groups)
    return [paths for _, paths in sorted(duplicates.items())]

def format_report(target_dir: Path, total_files: int, duplicates: Dict[str, List[Path]]) -> str:
    lines = []
    lines.append(f"Target directory: {target_dir.as_posix()}")
    lines.append(f"Scanned files: {total_files}")
    if not duplicates:
        lines.append("No duplicate Excel files were found.")
        return "\n".join(lines)
    lines.append("Duplicate Excel files detected:")
    for digest, paths in sorted(duplicates.items()):
        lines.append(f"Hash: {digest}")
        for path in paths:
            lines.append(f"  - {path.as_posix()}")
    return "\n".join(lines)

def write_report(report: str) -> None:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_FILE.write_text(report, encoding="utf-8")

def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Detect duplicate Excel files by MD5 hash.")
    parser.add_argument(
        "folder",
        nargs="?",
        default=None,
        help="Directory or subfolder name containing Excel files (default: raw_excel)",
    )
    return parser.parse_args(argv)

def resolve_excel_dir(folder_arg: str | None) -> Path:
    if not folder_arg:
        return DEFAULT_RAW_EXCEL_DIR
    candidate = Path(folder_arg)
    if candidate.exists():
        return candidate
    nested_candidate = DEFAULT_RAW_EXCEL_DIR / folder_arg
    if nested_candidate.exists():
        return nested_candidate
    return candidate

def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    excel_dir = resolve_excel_dir(args.folder)
    excel_files = list(iter_excel_files(excel_dir))
    if not excel_files:
        report = f"No Excel files were found in {excel_dir.as_posix()}."
        write_report(report)
        print(report)
        return
    groups, skipped = _safe_group_by_digest(excel_files)
    if skipped:
        for path, reason in skipped:
            print(f"Skip {path}: {reason}")
    duplicates = filter_duplicates(groups)
    report = format_report(excel_dir, len(excel_files), duplicates)
    write_report(report)
    print(report)

if __name__ == "__main__":
    main()
