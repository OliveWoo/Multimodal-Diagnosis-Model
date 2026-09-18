from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import date, datetime, time
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from tools.table1_sections import (  # noqa: E402
    GridRow,
    Table1Grid,
    TitledSection,
    find_row,
    load_table1_grid,
    normalize_text,
    split_titled_sections,
)
from utils import sanitize_filename  # noqa: E402


DEFAULT_SHEET_NAME = "Table 1"
SUPPORTED_EXTENSIONS = {".xlsx", ".xlsm"}

RAW_SHEETS: tuple[tuple[str, str], ...] = (
    ("underlying", "underlying"),
    ("admission diagnosis", "admission_diagnosis"),
    ("CBC", "CBC/otherLab"),
    ("other lab", "CBC/otherLab"),
    ("filmarray", "filmarray"),
    ("gm test", "filmarray"),
    ("molecular microbiology", "filmarray"),
    ("culture", "culture"),
    ("image", "image"),
)


def _normalize_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, datetime):
        if value.hour == value.minute == value.second == value.microsecond == 0:
            return value.strftime("%Y-%m-%d")
        if value.second == value.microsecond == 0:
            return value.strftime("%Y-%m-%d %H:%M")
        return value.strftime("%Y-%m-%d %H:%M:%S")
    if isinstance(value, date):
        return value.strftime("%Y-%m-%d")
    if isinstance(value, time):
        return value.strftime("%H:%M:%S" if value.second else "%H:%M")
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def _trim_row(values: list[str]) -> list[str]:
    while values and values[-1] == "":
        values.pop()
    return values


def _read_sheet_rows(sheet: Any) -> list[list[str]]:
    rows: list[list[str]] = [[]]
    for row in sheet.iter_rows():
        rows.append(_trim_row([normalize_text(_normalize_value(cell.value)) for cell in row]))
    return rows


def _raw_filename(workbook_stem: str, sheet_name: str) -> str:
    return f"{sanitize_filename(workbook_stem)}_{sanitize_filename(sheet_name)}_Raw.json"


def _default_output_dir(excel_path: Path) -> Path:
    return excel_path.with_name(f"{sanitize_filename(excel_path.stem)}_table1_json_Raw")


def _row_array(row: GridRow) -> list[str]:
    if not row.cells:
        return []
    values = [""] * (max(cell.col for cell in row.cells) + 1)
    for cell in row.cells:
        values[cell.col] = cell.text
    return _trim_row(values)


def _section_raw_rows(grid: Table1Grid, section: TitledSection) -> tuple[list[int], list[list[str]]]:
    if not section.row_numbers:
        return [], []
    mixed_row = any(grid.row(row).text.count("『") > 1 for row in section.row_numbers)
    if mixed_row:
        return [section.start_row], [[section.text]]
    pairs = [
        (row, _row_array(grid.row(row)))
        for row in range(section.start_row, section.end_row + 1)
        if grid.row(row).cells
    ]
    return [row for row, _ in pairs], [raw for _, raw in pairs]


def _classify_section(section: TitledSection) -> list[str]:
    title = section.title.casefold()
    output: list[str] = []
    if "一般血液檢驗" in title:
        output.append("CBC")
    elif any(
        token in title
        for token in (
            "床邊檢驗",
            "生化檢驗",
            "急診生化藥物檢驗",
            "血庫檢驗",
            "血液學檢驗",
            "臨床免疫檢驗",
            "風濕免疫內科檢驗",
            "生化電泳檢驗",
            "尿液檢驗",
            "ria 核醫檢驗",
            "糞便檢驗",
        )
    ):
        output.append("other lab")

    if "風濕免疫內科檢驗" in title:
        output.append("underlying")
    if "基因診斷檢查" in title:
        output.append("molecular microbiology")
    if "微生物分子生物學/核酸檢驗" in title or "各種體液檢驗" in title:
        output.append("filmarray")
    if any(token in title for token in ("細菌黴菌血清學檢驗", "小兒一般檢驗", "病毒血清學檢驗")):
        output.append("gm test")
    if re.search(r"培養檢驗|分離檢驗|抗酸菌直接鏡檢|革蘭氏染色鏡檢", title):
        output.append("culture")
    return output


def _append_rows(
    target: dict[str, dict[str, Any]],
    sheet_name: str,
    row_numbers: list[int],
    raw_rows: list[list[str]],
) -> None:
    payload = target[sheet_name]
    for row_number, raw_row in zip(row_numbers, raw_rows):
        marker = (row_number, json.dumps(raw_row, ensure_ascii=False))
        if marker in payload["seen"]:
            continue
        payload["seen"].add(marker)
        payload["source_rows"].append(row_number)
        payload["raw_data"].append(raw_row)


def _keyword_rows(grid: Table1Grid, patterns: tuple[str, ...], *, start: int = 1) -> tuple[list[int], list[list[str]]]:
    row_numbers: list[int] = []
    raw_rows: list[list[str]] = []
    for row in grid.iter_rows(start=start):
        if any(pattern.casefold() in row.text.casefold() for pattern in patterns):
            row_numbers.append(row.number)
            raw_rows.append(_row_array(row))
    return row_numbers, raw_rows


def parse_workbook(
    excel_path: Path,
    output_dir: Path | None = None,
    *,
    sheet_name: str = DEFAULT_SHEET_NAME,
) -> tuple[Path, dict[str, Any]]:
    excel_path = excel_path.expanduser().resolve()
    if not excel_path.exists():
        raise FileNotFoundError(f"Workbook not found: {excel_path}")
    if excel_path.suffix.lower() not in SUPPORTED_EXTENSIONS:
        raise ValueError(f"Unsupported workbook extension: {excel_path.suffix}")

    target_dir = (output_dir or _default_output_dir(excel_path)).expanduser().resolve()
    target_dir.mkdir(parents=True, exist_ok=True)
    grid = load_table1_grid(excel_path, sheet_name)
    image_start = find_row(grid, "Radiology / Special Examination Results") or find_row(grid, "放射報告與特殊檢查記錄")
    sections = split_titled_sections(grid, end_row=(image_start - 1 if image_start else None))

    grouped: dict[str, dict[str, Any]] = {
        sheet: {"raw_data": [], "source_rows": [], "seen": set()} for sheet, _ in RAW_SHEETS
    }
    for section in sections:
        row_numbers, raw_rows = _section_raw_rows(grid, section)
        for target_sheet in _classify_section(section):
            _append_rows(grouped, target_sheet, row_numbers, raw_rows)

    if image_start:
        pathology_row = find_row(grid, "Pathological Examination", start=image_start) or find_row(grid, "病理檢查報告", start=image_start)
        image_end = (pathology_row - 1) if pathology_row else len(grid.rows) - 1
        rows = list(range(image_start, image_end + 1))
        _append_rows(grouped, "image", rows, [_row_array(grid.row(row)) for row in rows])

    underlying_rows, underlying_data = _keyword_rows(
        grid,
        ("Systemic lupus erythematosus", "renal parenchymal disease", "Left renal cyst"),
    )
    _append_rows(grouped, "underlying", underlying_rows, underlying_data)
    admission_rows, admission_data = _keyword_rows(
        grid,
        ("Clinical Diagnosis", "診斷:", "診斷：", "Discharge Condition", "出院狀況", "死亡(expired)"),
    )
    _append_rows(grouped, "admission diagnosis", admission_rows, admission_data)

    manifest: dict[str, Any] = {
        "source_file": excel_path.name,
        "source_sheet": sheet_name,
        "output_dir": str(target_dir),
        "categories": {
            "CBC/otherLab": [],
            "filmarray": [],
            "culture": [],
            "image": [],
            "admission_diagnosis": [],
            "underlying": [],
        },
        "files": [],
    }

    for raw_sheet, category in RAW_SHEETS:
        data = grouped[raw_sheet]
        filename = _raw_filename(excel_path.stem, raw_sheet)
        payload = {"source_file": excel_path.name, "sheet": raw_sheet, "raw_data": data["raw_data"]}
        (target_dir / filename).write_text(
            json.dumps(payload, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        manifest["categories"][category].append(filename)
        manifest["files"].append(
            {
                "file": filename,
                "sheet": raw_sheet,
                "category": category,
                "row_count": len(data["raw_data"]),
                "source_rows": data["source_rows"],
            }
        )

    (target_dir / "table1_split_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return target_dir, manifest


def _iter_workbooks(target: Path) -> list[Path]:
    target = target.expanduser()
    if target.is_file():
        return [target]
    if target.is_dir():
        return sorted(path for path in target.iterdir() if path.suffix.lower() in SUPPORTED_EXTENSIONS)
    return []


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Split a one-sheet Table 1 workbook into section-aware legacy Raw JSON files.",
    )
    parser.add_argument("target", type=Path, help="An .xlsx/.xlsm file or a directory.")
    parser.add_argument("--output", type=Path, help="Output directory.")
    parser.add_argument("--sheet", default=DEFAULT_SHEET_NAME, help="Source sheet name.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    workbooks = _iter_workbooks(args.target)
    if not workbooks:
        print(f"No supported workbook found: {args.target}")
        return
    multiple = len(workbooks) > 1
    for workbook in workbooks:
        output_dir = None
        if args.output:
            output_dir = args.output / f"{sanitize_filename(workbook.stem)}_table1_json_Raw" if multiple else args.output
        target_dir, manifest = parse_workbook(workbook, output_dir, sheet_name=args.sheet)
        print(f"[Table1] {workbook.name} -> {target_dir}")
        for item in manifest["files"]:
            print(f"  - {item['sheet']} ({item['category']}): {item['row_count']} rows")


if __name__ == "__main__":
    main()
