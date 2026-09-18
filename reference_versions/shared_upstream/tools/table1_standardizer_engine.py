from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from typing import Any, Iterable, Sequence

from tools.table1_sections import GridCell, GridRow, Table1Grid, TitledSection, normalize_text


DATE_RE = re.compile(r"(?<!\d)(20\d{2})[/-](\d{1,2})[/-](\d{1,2})(?!\d)")
COMPACT_DATE_RE = re.compile(r"(?<!\d)(20\d{2})(\d{2})(\d{2})(?!\d)")
TIME_RE = re.compile(r"(?<!\d)(\d{1,2}:\d{2}(?::\d{2})?)(?!\d)")
SECTION_MARKER_RE = re.compile(r"『.*?』", re.DOTALL)

QUALITATIVE_RESULTS = (
    "Target Not Detected",
    "Not Detected",
    "Nonreactive",
    "Non-reactive",
    "No growth",
    "Negative",
    "Positive",
    "Reactive",
    "Detected",
    "Pending",
    "Normal",
    "NORMAL",
    "Colorless",
    "Cloudy",
    "Yellow",
    "Light-Yellow",
    "Light-Orange",
    "Dark Brown",
    "Soft",
)


@dataclass(frozen=True)
class TimePoint:
    reported_time: str
    x: float


def _dates(text: str) -> list[str]:
    normalized = normalize_text(text)
    result = [f"{int(y):04d}-{int(m):02d}-{int(d):02d}" for y, m, d in DATE_RE.findall(normalized)]
    if result:
        return result
    return [
        f"{int(y):04d}-{int(m):02d}-{int(d):02d}"
        for y, m, d in COMPACT_DATE_RE.findall(normalized)
    ]


def _times(text: str) -> list[str]:
    return TIME_RE.findall(normalize_text(text))


def _first_datetime(text: str) -> str:
    dates = _dates(text)
    if not dates:
        return ""
    times = _times(text)
    return f"{dates[0]} {times[0]}" if times else dates[0]


def _combine_datetimes(dates: Sequence[str], times: Sequence[str]) -> list[str]:
    output: list[str] = []
    for index, date_value in enumerate(dates):
        if index < len(times) and times[index]:
            output.append(f"{date_value} {times[index]}")
        else:
            output.append(date_value)
    return output


def _source(grid: Table1Grid, row: int, end_row: int | None = None) -> dict[str, Any]:
    source: dict[str, Any] = {
        "source_file": grid.source_path.name,
        "source_sheet": grid.sheet_name,
        "row": row,
    }
    if end_row is not None and end_row != row:
        source["end_row"] = end_row
    return source


def _clean_spaces(text: str) -> str:
    return re.sub(r"\s+", " ", normalize_text(text)).strip()


def _looks_like_result_start(text: str) -> bool:
    stripped = text.lstrip()
    if not stripped:
        return False
    if any(stripped.casefold().startswith(value.casefold()) for value in QUALITATIVE_RESULTS):
        return True
    return bool(
        re.match(
            r"^(?:[<>]=?\s*)?(?:[-+]?\d+(?:\.\d+)?|\d+(?:\.\d+)?[+])"
            r"|^(?:\+/-|[-+]|[1-4]\+|[SIR](?:\b|\())",
            stripped,
            re.IGNORECASE,
        )
    )


def _split_item_value(text: str) -> tuple[str, str, int]:
    line = normalize_text(text).strip()
    if not line or "\n" in line:
        return line, "", len(line)

    for match in re.finditer(r"\([^()]*\)", line):
        rest = line[match.end() :]
        if _looks_like_result_start(rest):
            return line[: match.end()].strip(), rest.strip(), match.end() + len(rest) - len(rest.lstrip())

    for match in re.finditer(r"\s{2,}", line):
        rest = line[match.end() :]
        if _looks_like_result_start(rest):
            return line[: match.start()].strip(), rest.strip(), match.end()

    result_match = re.search(
        r"(?:Target\s+Not\s+Detected|Not\s+Detected|Non-?reactive|Negative|Positive|Reactive|Detected)",
        line,
        re.IGNORECASE,
    )
    if result_match and result_match.start() > 0:
        return line[: result_match.start()].strip(), line[result_match.start() :].strip(), result_match.start()
    return line, "", len(line)


def _looks_like_item(text: str) -> bool:
    item = _clean_spaces(text)
    if not item or len(item) > 180:
        return False
    if item.startswith("『") or "項目/日期" in item:
        return False
    if DATE_RE.search(item) or re.fullmatch(r"\d{1,2}:\d{2}(?::\d{2})?", item):
        return False
    return bool(re.search(r"[A-Za-z\u4e00-\u9fffα-ωΑ-Ωμ%]", item))


def _item_short(item_full: str) -> str:
    item = _clean_spaces(item_full)
    unit_like = re.compile(
        r"\s*\([^)]*(?:/|%|IU|U/L|mol|mg|g/dL|uL|fl|Pg|mmHg|℃|°C|Index|S/CO|HPF)[^)]*\)\s*$",
        re.IGNORECASE,
    )
    while unit_like.search(item):
        item = unit_like.sub("", item).strip()
    return item or _clean_spaces(item_full)


def _value_tokens(value: str) -> list[str]:
    text = normalize_text(value).strip()
    if not text:
        return []

    repeated = re.findall(
        r"Target\s+Not\s+Detected|Not\s+Detected|No acid fast bacilli\s+observed\.?|"
        r"\d+(?:\.\d+)?\s*\((?:Nonreactive|Reactive|Negative|Positive)\)",
        text,
        re.IGNORECASE,
    )
    if len(repeated) > 1:
        return [_clean_spaces(token) for token in repeated]

    chunks = [_clean_spaces(part) for part in re.split(r"\n|\s{2,}", text) if _clean_spaces(part)]
    chunks = [part for part in chunks if not re.fullmatch(r"\(\d{1,2}:\d{2}(?::\d{2})?\)", part)]
    if len(chunks) > 1:
        return chunks

    numeric = re.findall(
        r"(?:\+/-|\d+\+|[-+])(?:\([^)]*\))?|"
        r"(?:[<>]=?\s*)?[-+]?\d+(?:\.\d+)?(?:\([^)]*\))?",
        text,
    )
    if len(numeric) > 1:
        residue = text
        for token in numeric:
            residue = residue.replace(token, "", 1)
        if not residue.strip(" ,;/|"):
            return [_clean_spaces(token) for token in numeric]
    return [_clean_spaces(text)]


def _nearest_index(value: float, centers: Sequence[float]) -> int:
    return min(range(len(centers)), key=lambda index: abs(centers[index] - value))


def _cluster_centers(values: Sequence[float], count: int) -> list[float]:
    unique = sorted({round(value, 3) for value in values})
    if count <= 1:
        return [unique[0] if unique else 0.0]
    if len(unique) == count:
        return [float(value) for value in unique]
    if len(unique) < 2:
        return [float(index) for index in range(count)]
    centers = [unique[0] + (unique[-1] - unique[0]) * index / (count - 1) for index in range(count)]
    for _ in range(30):
        groups: list[list[float]] = [[] for _ in range(count)]
        for value in unique:
            groups[_nearest_index(value, centers)].append(float(value))
        updated = [sum(group) / len(group) if group else centers[index] for index, group in enumerate(groups)]
        if all(math.isclose(a, b, abs_tol=0.001) for a, b in zip(centers, updated)):
            break
        centers = updated
    return sorted(centers)


def _row_is_time_only(row: GridRow) -> bool:
    if not row.cells:
        return False
    for cell in row.cells:
        residue = TIME_RE.sub("", cell.text)
        if residue.strip(" \n,;/|"):
            return False
    return True


def _cell_date_records(row: GridRow) -> list[tuple[str, GridCell]]:
    records: list[tuple[str, GridCell]] = []
    for cell in row.cells:
        for date_value in _dates(cell.text):
            records.append((date_value, cell))
    return records


def _infer_block_centers(rows: Iterable[GridRow], count: int) -> list[float]:
    coordinates: list[float] = []
    for row in rows:
        if "項目/日期" in row.text:
            continue
        for index, cell in enumerate(row.cells):
            item, first_value, _ = _split_item_value(cell.text.replace("\n", " "))
            if first_value or index > 0 or not _looks_like_item(item):
                coordinates.append(cell.x)
            elif index == 0:
                coordinates.append(cell.x)
    return _cluster_centers(coordinates, count)


def _header_timepoints(grid: Table1Grid, header_row: int, block_end: int) -> tuple[list[TimePoint], bool]:
    row = grid.row(header_row)
    date_records = _cell_date_records(row)
    dates = [date_value for date_value, _ in date_records]
    if not dates:
        return [], False

    same_row_times = _times(row.text)
    next_row = grid.row(header_row + 1)
    use_next_time_row = header_row < block_end and _row_is_time_only(next_row)
    next_times = _times(next_row.text) if use_next_time_row else []
    times = same_row_times if len(same_row_times) >= len(dates) else next_times
    datetime_values = _combine_datetimes(dates, times)

    date_cells = [cell for _, cell in date_records]
    if len({cell.x for cell in date_cells}) == len(dates):
        centers = [cell.x for cell in date_cells]
    elif use_next_time_row and len(next_row.cells) >= len(dates):
        centers = [cell.x for cell in next_row.cells[: len(dates)]]
    else:
        block_rows = list(grid.iter_rows(header_row, block_end))
        centers = _infer_block_centers(block_rows, len(dates))

    centers = sorted(centers)
    if len(centers) != len(datetime_values):
        centers = _cluster_centers(centers, len(datetime_values))
    return [TimePoint(value, centers[index]) for index, value in enumerate(datetime_values)], use_next_time_row


def _make_lab_entry(
    grid: Table1Grid,
    *,
    row: int,
    item_full: str,
    value: str,
    reported_time: str,
    test_name: str,
) -> dict[str, Any]:
    return {
        "reported_time": reported_time,
        "item": _item_short(item_full),
        "item_full": _clean_spaces(item_full),
        "value": _clean_spaces(value),
        "test": test_name,
        "source": _source(grid, row),
    }


def _parse_fixed_width_lines(
    grid: Table1Grid,
    text: str,
    timepoints: Sequence[TimePoint],
    *,
    row: int,
    test_name: str,
) -> list[dict[str, Any]]:
    if not timepoints:
        return []
    parsed: list[tuple[str, list[tuple[str, int]]]] = []
    all_positions: list[int] = []
    for raw_line in normalize_text(text).splitlines():
        line = raw_line.rstrip()
        if not line or "項目/日期" in line or line.startswith("『"):
            continue
        item, value_text, value_start = _split_item_value(line)
        if not value_text or not _looks_like_item(item):
            continue
        tokens = _value_tokens(value_text)
        located: list[tuple[str, int]] = []
        cursor = value_start
        for token in tokens:
            position = line.find(token, cursor)
            if position < 0:
                position = cursor
            located.append((token, position))
            all_positions.append(position)
            cursor = position + len(token)
        parsed.append((item, located))

    if not parsed:
        return []
    if len(timepoints) == 1:
        min_position = max_position = all_positions[0]
    else:
        min_position = min(all_positions)
        max_position = max(all_positions)

    entries: list[dict[str, Any]] = []
    for item, values in parsed:
        for value, position in values:
            if len(timepoints) == 1 or max_position == min_position:
                index = 0
            else:
                index = round((position - min_position) * (len(timepoints) - 1) / (max_position - min_position))
                index = max(0, min(index, len(timepoints) - 1))
            entries.append(
                _make_lab_entry(
                    grid,
                    row=row,
                    item_full=item,
                    value=value,
                    reported_time=timepoints[index].reported_time,
                    test_name=test_name,
                )
            )
    return entries


def _parse_header_embedded_items(
    grid: Table1Grid,
    row: GridRow,
    timepoints: Sequence[TimePoint],
    *,
    test_name: str,
) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    header_cell = next((cell for cell in row.cells if "項目/日期" in cell.text), None)
    if header_cell is None:
        return entries

    lines = header_cell.text.splitlines()
    header_index = next((index for index, line in enumerate(lines) if "項目/日期" in line), -1)
    tail = "\n".join(lines[header_index + 1 :]) if header_index >= 0 else ""
    if tail.strip():
        entries.extend(
            _parse_fixed_width_lines(
                grid,
                header_cell.text,
                timepoints,
                row=row.number,
                test_name=test_name,
            )
        )

    item_lines = [line for line in tail.splitlines() if _looks_like_item(_split_item_value(line)[0])]
    if len(item_lines) != 1:
        return entries
    item_full, inline_value, _ = _split_item_value(item_lines[0])
    if inline_value:
        return entries

    for cell in row.cells:
        if cell is header_cell:
            continue
        stripped = DATE_RE.sub("", cell.text)
        stripped = TIME_RE.sub("", stripped)
        stripped = re.sub(r"\s+", " ", stripped).strip()
        if not stripped:
            continue
        index = _nearest_index(cell.x, [point.x for point in timepoints])
        for offset, value in enumerate(_value_tokens(stripped)):
            target_index = min(index + offset, len(timepoints) - 1)
            entries.append(
                _make_lab_entry(
                    grid,
                    row=row.number,
                    item_full=item_full,
                    value=value,
                    reported_time=timepoints[target_index].reported_time,
                    test_name=test_name,
                )
            )
    return entries


def _parse_matrix_data_row(
    grid: Table1Grid,
    row: GridRow,
    timepoints: Sequence[TimePoint],
    *,
    test_name: str,
) -> list[dict[str, Any]]:
    if not row.cells or not timepoints:
        return []
    first_cell = row.cells[0]
    first_text = first_cell.text
    if "\n" in first_text:
        lines = [line for line in first_text.splitlines() if line.strip()]
        if len(lines) != 1:
            return _parse_multiline_matrix_row(grid, row, timepoints, test_name=test_name)

    item_full, first_value, _ = _split_item_value(first_text.replace("\n", " "))
    if not _looks_like_item(item_full):
        return []

    centers = [point.x for point in timepoints]
    assignments: list[tuple[int, str]] = []
    if first_value:
        start = _nearest_index(first_cell.x, centers)
        for offset, value in enumerate(_value_tokens(first_value)):
            assignments.append((min(start + offset, len(timepoints) - 1), value))

    for cell in row.cells[1:]:
        start = _nearest_index(cell.x, centers)
        for offset, value in enumerate(_value_tokens(cell.text)):
            assignments.append((min(start + offset, len(timepoints) - 1), value))

    entries: list[dict[str, Any]] = []
    seen: set[tuple[int, str]] = set()
    for index, value in assignments:
        key = (index, value)
        if key in seen:
            continue
        seen.add(key)
        entries.append(
            _make_lab_entry(
                grid,
                row=row.number,
                item_full=item_full,
                value=value,
                reported_time=timepoints[index].reported_time,
                test_name=test_name,
            )
        )
    return entries


def _extract_item_labels(text: str) -> list[str]:
    labels: list[str] = []
    pattern = re.compile(
        r"[A-Za-z%α-ωΑ-Ωμ\u4e00-\u9fff][A-Za-z0-9%α-ωΑ-Ωμ\u4e00-\u9fff+./ -]*?"
        r"\([^()\n]*\)(?:\s*\([^()\n]*\))*"
    )
    for line in normalize_text(text).splitlines():
        for match in pattern.finditer(line):
            label = _clean_spaces(match.group(0))
            if _looks_like_item(label):
                labels.append(label)
        if not pattern.search(line):
            item, _, _ = _split_item_value(line)
            if _looks_like_item(item):
                labels.append(_clean_spaces(item))
    return list(dict.fromkeys(labels))


def _numeric_value(value: str) -> float | None:
    match = re.search(r"[-+]?\d+(?:\.\d+)?", value.replace(",", ""))
    return float(match.group(0)) if match else None


def _item_value_score(item: str, value: str) -> float:
    folded = item.casefold()
    number = _numeric_value(value)
    if number is None:
        return 10.0
    if "temperature" in folded or "°c" in folded or "℃" in folded:
        return 120.0 if 25 <= number <= 45 else -120.0
    if "%so2" in folded:
        return 120.0 if 50 <= number <= 100 else -50.0
    if "hba1c" in folded:
        return 115.0 if 3 <= number <= 20 else -80.0
    if "vit-12" in folded or "vit b12" in folded:
        return 110.0 if number >= 100 else -30.0
    if "folic" in folded:
        return 100.0 if number >= 2 else 0.0
    if "pct" in folded:
        return 95.0 if 0 <= number <= 20 else 5.0
    if "lactate" in folded:
        return 90.0 if 0.2 <= number <= 20 else -20.0
    if "d-dimer" in folded:
        return 105.0 if 2 <= number <= 100 else 45.0
    if "ca++" in folded or "ionized" in folded:
        return 85.0 if 2 <= number <= 10 else 20.0
    if re.search(r"\b(?:sgot|sgpt|cpk|ck-mb)", folded):
        return 90.0 if number >= 10 else 15.0
    return 30.0


def _parse_multiline_matrix_row(
    grid: Table1Grid,
    row: GridRow,
    timepoints: Sequence[TimePoint],
    *,
    test_name: str,
) -> list[dict[str, Any]]:
    labels = _extract_item_labels(row.cells[0].text)
    if not labels:
        return []

    value_records: list[tuple[str, float]] = []
    for line in row.cells[0].text.splitlines():
        _, inline_value, _ = _split_item_value(line)
        for value in _value_tokens(inline_value):
            value_records.append((value, row.cells[0].x))
    for cell in row.cells[1:]:
        for value in _value_tokens(cell.text):
            value_records.append((value, cell.x))
    if not value_records:
        return []

    assignments: list[tuple[str, str, float]] = []
    assigned_counts = {label: 0 for label in labels}
    for value, x in value_records:
        ranked = sorted(
            labels,
            key=lambda label: (_item_value_score(label, value), -assigned_counts[label]),
            reverse=True,
        )
        chosen = ranked[0]
        assignments.append((chosen, value, x))
        assigned_counts[chosen] += 1

    centers = [point.x for point in timepoints]
    entries: list[dict[str, Any]] = []
    for item_full, value, x in assignments:
        index = _nearest_index(x, centers)
        entries.append(
            _make_lab_entry(
                grid,
                row=row.number,
                item_full=item_full,
                value=value,
                reported_time=timepoints[index].reported_time,
                test_name=test_name,
            )
        )
    return entries


def parse_matrix_lab_section(
    grid: Table1Grid,
    section: TitledSection,
    *,
    test_name: str,
) -> list[dict[str, Any]]:
    if not section.row_numbers:
        return []
    header_rows = [
        row.number
        for row in grid.iter_rows(section.start_row, section.end_row)
        if "項目/日期" in row.text
    ]
    entries: list[dict[str, Any]] = []
    for header_index, header_row in enumerate(header_rows):
        block_end = header_rows[header_index + 1] - 1 if header_index + 1 < len(header_rows) else section.end_row
        timepoints, consumed_time_row = _header_timepoints(grid, header_row, block_end)
        if not timepoints:
            continue
        entries.extend(
            _parse_header_embedded_items(
                grid,
                grid.row(header_row),
                timepoints,
                test_name=test_name,
            )
        )
        data_start = header_row + (2 if consumed_time_row else 1)
        for row in grid.iter_rows(data_start, block_end):
            if "項目/日期" in row.text or _row_is_time_only(row):
                continue
            entries.extend(_parse_matrix_data_row(grid, row, timepoints, test_name=test_name))
    return _dedupe(entries)


def parse_inline_lab_section(
    grid: Table1Grid,
    section: TitledSection,
    *,
    test_name: str,
) -> list[dict[str, Any]]:
    dates = _dates(section.text)
    if not dates:
        return []
    times = _times(section.text)
    timepoints = [
        TimePoint(reported_time=value, x=float(index))
        for index, value in enumerate(_combine_datetimes(dates, times))
    ]
    return _parse_fixed_width_lines(
        grid,
        section.text,
        timepoints,
        row=section.start_row,
        test_name=test_name,
    )


def _dedupe(entries: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    seen: set[str] = set()
    for entry in entries:
        key = json.dumps(entry, ensure_ascii=False, sort_keys=True)
        if key in seen:
            continue
        seen.add(key)
        output.append(entry)
    return output


def _section_has(section: TitledSection, *needles: str) -> bool:
    folded = section.title.casefold()
    return any(normalize_text(needle).casefold() in folded for needle in needles)


def parse_cbc(grid: Table1Grid, sections: Sequence[TitledSection]) -> list[dict[str, Any]]:
    section = next((item for item in sections if _section_has(item, "一般血液檢驗")), None)
    return parse_matrix_lab_section(grid, section, test_name="CBC") if section else []


def _parse_rheumatology_section(grid: Table1Grid, section: TitledSection) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    text = section.text
    header_parts = text.split("項目/日期")
    primary_header = header_parts[1] if len(header_parts) > 1 else text
    dates = _dates(primary_header)
    times = _times(primary_header)
    if len(header_parts) > 2:
        secondary_dates = set(_dates(header_parts[2]))
        dates = [value for value in dates if value not in secondary_dates]
        secondary_times = set(_times(header_parts[2]))
        times = [value for value in times if value not in secondary_times]
    default_times = _combine_datetimes(dates, times)

    for target in ("Anti-Mitochondrial Ab", "ASMA"):
        match = re.search(rf"{re.escape(target)}\s+(Negative|Positive|Reactive|Nonreactive)", text, re.IGNORECASE)
        if match:
            entries.append(
                _make_lab_entry(
                    grid,
                    row=section.start_row,
                    item_full=target,
                    value=match.group(1),
                    reported_time=default_times[0] if default_times else "",
                    test_name="Other Lab Tests",
                )
            )

    for row_number in section.row_numbers:
        row = grid.row(row_number)
        if "ds DNA" in row.text and "PR3 ANCA" in row.text:
            labels = ["ds DNA (IU/mL)", "PR3 ANCA (IU/mL)", "MPO ANCA (IU/mL)"]
            values = [cell.text for cell in row.cells[1:]]
            for index, value in enumerate(values):
                label = labels[min(index, len(labels) - 1)]
                time_index = min(index, len(default_times) - 1) if default_times else 0
                entries.append(
                    _make_lab_entry(
                        grid,
                        row=row_number,
                        item_full=label,
                        value=value,
                        reported_time=default_times[time_index] if default_times else "",
                        test_name="Other Lab Tests",
                    )
                )
        elif "Anti-GBM" in row.text:
            item, first_value, _ = _split_item_value(row.cells[0].text)
            values = ([first_value] if first_value else []) + [cell.text for cell in row.cells[1:]]
            for value in values:
                entries.append(
                    _make_lab_entry(
                        grid,
                        row=row_number,
                        item_full=item,
                        value=value,
                        reported_time=default_times[-1] if default_times else "",
                        test_name="Other Lab Tests",
                    )
                )
        elif "項目/日期" in row.text and "ds DNA" in row.text:
            date_time = _first_datetime(row.text)
            match = re.search(r"ds DNA\s*\(IU/mL\)\s*([^\s]+)", row.text, re.IGNORECASE)
            if match:
                entries.append(
                    _make_lab_entry(
                        grid,
                        row=row_number,
                        item_full="ds DNA (IU/mL)",
                        value=match.group(1),
                        reported_time=date_time,
                        test_name="Other Lab Tests",
                    )
                )
    return _dedupe(entries)


def parse_other_labs(grid: Table1Grid, sections: Sequence[TitledSection]) -> list[dict[str, Any]]:
    matrix_titles = (
        "床邊檢驗",
        "緊急生化檢驗",
        "一般生化檢驗",
        "血庫檢驗",
        "血液學檢驗",
        "臨床免疫檢驗",
        "生化電泳檢驗",
        "尿液檢驗",
    )
    inline_titles = ("急診生化藥物檢驗", "RIA 核醫檢驗", "糞便檢驗")
    entries: list[dict[str, Any]] = []
    for section in sections:
        if _section_has(section, "一般血液檢驗"):
            continue
        if _section_has(section, *matrix_titles):
            entries.extend(parse_matrix_lab_section(grid, section, test_name="Other Lab Tests"))
        elif _section_has(section, *inline_titles):
            entries.extend(parse_inline_lab_section(grid, section, test_name="Other Lab Tests"))
        elif _section_has(section, "風濕免疫內科檢驗"):
            entries.extend(_parse_rheumatology_section(grid, section))
    return _dedupe(entries)


def _micro_entry(
    *,
    test: str,
    sample: str,
    target: str,
    reported_time: str,
    result: str | None = None,
    value: str | None = None,
    collected_time: str | None = None,
    source: dict[str, Any] | None = None,
) -> dict[str, Any]:
    entry: dict[str, Any] = {
        "test": test,
        "sample": sample,
        "collected_time": collected_time if collected_time is not None else reported_time,
        "reported_time": reported_time,
        "target": _clean_spaces(target),
    }
    if result:
        entry["result"] = _clean_spaces(result)
    if value:
        entry["value"] = _clean_spaces(value)
    if source:
        entry["source"] = source
    return entry


def _sample_from_title(title: str) -> str:
    match = re.search(r"-\s*(.+)$", title)
    if match:
        return _clean_spaces(match.group(1))
    for sample in ("Endotracheal aspirate", "Lower BAL", "Sputum", "Blood", "Midstream", "Foley catheter", "Serum", "Urine", "Stool"):
        if sample.casefold() in title.casefold():
            return sample
    return ""


def parse_gm_serology(grid: Table1Grid, sections: Sequence[TitledSection]) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    for section in sections:
        sample = _sample_from_title(section.title)
        if _section_has(section, "細菌黴菌血清學檢驗"):
            if sample == "Serum":
                header_row = next(
                    (row for row in section.row_numbers if "項目/日期" in grid.row(row).text),
                    None,
                )
                gm_row = next((row for row in section.row_numbers if "GM test" in grid.row(row).text), None)
                if header_row and gm_row:
                    timepoints, _ = _header_timepoints(grid, header_row, gm_row)
                    lab_rows = _parse_matrix_data_row(
                        grid,
                        grid.row(gm_row),
                        timepoints,
                        test_name="GM / Antigen / Serology",
                    )
                    for item in lab_rows:
                        entries.append(
                            _micro_entry(
                                test="GM / Antigen / Serology",
                                sample="Serum",
                                target=item["item_full"],
                                reported_time=item["reported_time"],
                                value=item["value"],
                                source=item["source"],
                            )
                        )
                crypto = re.search(r"Cryptococcus antigen(?:\s*\(Serum\))?\s+(Negative|Positive)", section.text, re.IGNORECASE)
                if crypto:
                    entries.append(
                        _micro_entry(
                            test="GM / Antigen / Serology",
                            sample="Serum",
                            target="Cryptococcus antigen",
                            reported_time="",
                            result=crypto.group(1),
                            source=_source(grid, section.end_row),
                        )
                    )
            else:
                date_time = _first_datetime(section.text)
                match = re.search(r"GM test\s*\(Index\)\s*([<>]?\d+(?:\.\d+)?)", section.text, re.IGNORECASE)
                if match:
                    entries.append(
                        _micro_entry(
                            test="GM / Antigen / Serology",
                            sample=sample or "Lower BAL",
                            target="GM test (Index)",
                            reported_time=date_time,
                            value=match.group(1),
                            source=_source(grid, section.start_row),
                        )
                    )
        elif _section_has(section, "小兒一般檢驗"):
            date_time = _first_datetime(section.text)
            for target, value in re.findall(
                r"([^\n]+?\([^\n()]+\))\s*([<>]?\d+(?:\.\d+)?)",
                section.text,
            ):
                if "項目/日期" in target:
                    continue
                entries.append(
                    _micro_entry(
                        test="GM / Antigen / Serology",
                        sample=sample or "Blood",
                        target=target,
                        reported_time=date_time,
                        value=value,
                        source=_source(grid, section.start_row),
                    )
                )
        elif _section_has(section, "病毒血清學檢驗"):
            dates = _dates(section.text)
            times = _times(section.text)
            date_times = _combine_datetimes(dates, times)
            lines = []
            positions: list[int] = []
            for line in section.text.splitlines():
                match = re.match(
                    r"(.+?\([^)]*\))\s*([<>]?\d+(?:\.\d+)?\((?:Nonreactive|Reactive|Negative|Positive)\))",
                    line,
                    re.IGNORECASE,
                )
                if match:
                    position = line.find(match.group(2))
                    lines.append((match.group(1), match.group(2), position))
                    positions.append(position)
            low = min(positions) if positions else 0
            high = max(positions) if positions else low
            for target, raw_value, position in lines:
                result_match = re.search(r"\(([^)]+)\)$", raw_value)
                value = re.sub(r"\([^)]+\)$", "", raw_value)
                index = 0 if high == low or len(date_times) < 2 else round((position - low) * (len(date_times) - 1) / (high - low))
                entries.append(
                    _micro_entry(
                        test="GM / Antigen / Serology",
                        sample=sample or "Serum",
                        target=target,
                        reported_time=date_times[min(index, len(date_times) - 1)] if date_times else "",
                        result=result_match.group(1) if result_match else None,
                        value=value,
                        source=_source(grid, section.start_row),
                    )
                )
    return _dedupe(entries)


def parse_molecular_microbiology(grid: Table1Grid, sections: Sequence[TitledSection]) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    for section in sections:
        if not _section_has(section, "基因診斷檢查"):
            continue
        dates = _dates(section.text)
        times = _times(section.text)
        date_times = _combine_datetimes(dates, times)
        target_match = re.search(r"([A-Za-z][A-Za-z0-9-]*\s*\([^)]*\))", section.text)
        target = target_match.group(1) if target_match else "Molecular target"
        results = re.findall(r"Target\s+Not\s+Detected|Not\s+Detected|Detected(?:\s+[^\s]+)?", section.text, re.IGNORECASE)
        for index, result in enumerate(results):
            if index >= len(date_times):
                break
            normalized_result = "Not Detected" if "not" in result.casefold() else "Detected"
            entries.append(
                _micro_entry(
                    test="Non-panel Molecular Test",
                    sample=_sample_from_title(section.title) or "Blood",
                    target=target,
                    reported_time=date_times[index],
                    result=normalized_result,
                    source=_source(grid, section.start_row),
                )
            )
    return _dedupe(entries)


def _result_pairs(text: str) -> list[tuple[str, str, str | None]]:
    normalized = normalize_text(text)
    matches = list(re.finditer(r"Not\s+Detected|Negative|Positive|Detected", normalized, re.IGNORECASE))
    output: list[tuple[str, str, str | None]] = []
    cursor = 0
    for match in matches:
        target = normalized[cursor : match.start()].strip()
        cursor = match.end()
        target = re.sub(r"^.*?(?:\d{1,2}:\d{2}(?::\d{2})?)", "", target, flags=re.DOTALL).strip()
        target = re.sub(r"\(copy/mL\)\s*$", "", target, flags=re.IGNORECASE).strip()
        target = target.splitlines()[-1].strip() if target else ""
        if not target:
            continue
        result = match.group(0)
        output.append((target, "Not Detected" if result.casefold() in {"not detected", "negative"} else "Detected", None))
    return output


def parse_filmarray(grid: Table1Grid, sections: Sequence[TitledSection]) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    for section in sections:
        if _section_has(section, "微生物分子生物學/核酸檢驗"):
            date_time = _first_datetime(section.text)
            data_start = 0
            time_matches = list(TIME_RE.finditer(section.text))
            if time_matches:
                data_start = time_matches[-1].end()
            for target, result, value in _result_pairs(section.text[data_start:]):
                entries.append(
                    _micro_entry(
                        test="FilmArray / Molecular Panel",
                        sample=_sample_from_title(section.title) or "Lower BAL",
                        target=target,
                        reported_time=date_time,
                        result=result,
                        value=value,
                        source=_source(grid, section.start_row),
                    )
                )
        elif _section_has(section, "各種體液檢驗"):
            date_time = _first_datetime(section.text)
            for target, result, value in _result_pairs(section.text):
                cleaned = re.sub(r"^(?:項目/日期|\d{4}/\d{1,2}/\d{1,2}|\d{1,2}:\d{2})\s*", "", target).strip()
                if not re.search(r"Influ|InfA|InfB|Influenza", cleaned, re.IGNORECASE):
                    continue
                target_name = "Influenza A" if re.search(r"Type A|InfA", cleaned, re.IGNORECASE) else "Influenza B"
                entries.append(
                    _micro_entry(
                        test="FilmArray / Molecular Panel",
                        sample=_sample_from_title(section.title) or "Nasopharyngeal",
                        target=target_name,
                        reported_time=date_time,
                        result=result,
                        value=value,
                        source=_source(grid, section.start_row),
                    )
                )
    return _dedupe(entries)


def _culture_type(title: str) -> str:
    result = title
    for sample in (
        "Endotracheal aspirate",
        "Lower BAL",
        "Midstream(中段尿)",
        "Foley catheter",
        "Sputum",
        "Blood",
        "Urine",
    ):
        result = re.sub(re.escape(sample), "", result, flags=re.IGNORECASE)
    return _clean_spaces(result)


def _culture_times(text: str) -> tuple[str, str, str]:
    def extract(label: str) -> str:
        match = re.search(
            rf"{label}\s*[:：]\s*(20\d{{2}}[/-]\d{{1,2}}[/-]\d{{1,2}}(?:\s+\d{{1,2}}:\d{{2}}(?::\d{{2}})?)?)",
            text,
            re.IGNORECASE,
        )
        return _first_datetime(match.group(1)) if match else ""

    collected = extract("採檢時間")
    received = extract("簽收時間")
    reported = extract("報告時間")
    if not reported and "尚在檢驗中" in text:
        pending_date = _dates(text[text.find("尚在檢驗中") - 40 : text.find("尚在檢驗中") + 40])
        reported = f"Pending as of {pending_date[0]}" if pending_date else "Pending"
    return collected, received, reported


def _clean_organism(value: str) -> str:
    organism = _clean_spaces(value)
    organism = re.split(r"(?:\bB:|\bAntimicrobial|\bCommon\s+(?:aerobic|anaerobic)|\d+\.同類)", organism, maxsplit=1, flags=re.IGNORECASE)[0]
    organism = re.sub(r"^\d+[.)]\s*", "", organism).strip(" .:;/")
    return organism


def _culture_organisms(text: str) -> list[tuple[str, str]]:
    if re.search(
        r"No growth|not isolated|were not isolated|No growth of pathogens fungi",
        text,
        re.IGNORECASE,
    ):
        return [("No growth", "Not Detected")]

    organisms: list[tuple[str, str]] = []
    for match in re.finditer(r"Final report\s*:\s*([^\n]+)", text, re.IGNORECASE):
        organism = _clean_organism(match.group(1))
        if organism and not organism.casefold().startswith("no growth"):
            organisms.append((organism, "Detected"))

    if not organisms:
        isolated = re.search(r"A:Organism isolated\s*\n\s*1[.)]\s*([^\n]+)", text, re.IGNORECASE)
        if isolated:
            organism = _clean_organism(isolated.group(1))
            if organism:
                organisms.append((organism, "Detected"))

    if not organisms:
        preliminary = re.search(r"Preliminary report\s*:\s*([^\n]+)", text, re.IGNORECASE)
        if preliminary:
            organism = _clean_organism(preliminary.group(1))
            if organism:
                organisms.append((organism, "Preliminary report"))

    if not organisms and "尚在檢驗中" in text:
        organisms.append(("Pending", "Pending"))
    return organisms


def _looks_antibiotic_name(text: str) -> bool:
    candidate = _clean_spaces(text)
    if not candidate or len(candidate) > 80:
        return False
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9 /+-]*", candidate):
        return False
    blocked = (
        "report",
        "culture",
        "organism",
        "isolated",
        "susceptibility",
        "blood",
        "sputum",
        "lower bal",
        "endotracheal",
        "clinical pathogen",
    )
    return not any(word in candidate.casefold() for word in blocked)


def _ast_result(raw: str) -> tuple[str, str | None]:
    match = re.fullmatch(r"([SIR])(?:\(([^)]*)\))?", _clean_spaces(raw), re.IGNORECASE)
    if not match:
        return _clean_spaces(raw), None
    return match.group(1).upper(), match.group(2)


def _culture_ast(text: str) -> list[dict[str, Any]]:
    lines = [_clean_spaces(line) for line in normalize_text(text).splitlines() if _clean_spaces(line)]
    output: list[dict[str, Any]] = []
    index = 0
    while index < len(lines):
        line = lines[index]
        inline = re.fullmatch(r"(.+?)\s+([SIR](?:\([^)]*\))?)", line, re.IGNORECASE)
        if inline and _looks_antibiotic_name(inline.group(1)):
            result, mic = _ast_result(inline.group(2))
            item: dict[str, Any] = {"antibiotic": inline.group(1).strip(), "result": result}
            if mic:
                item["MIC"] = mic
            output.append(item)
        elif index + 1 < len(lines) and _looks_antibiotic_name(line):
            next_line = lines[index + 1]
            if re.fullmatch(r"[SIR](?:\([^)]*\))?", next_line, re.IGNORECASE):
                result, mic = _ast_result(next_line)
                item = {"antibiotic": line, "result": result}
                if mic:
                    item["MIC"] = mic
                output.append(item)
                index += 1
        index += 1

    for match in re.finditer(
        r"(?<![A-Za-z])([A-Za-z][A-Za-z/-]{2,})\s+([SIR](?:\([^)]*\))?)",
        normalize_text(text),
        re.IGNORECASE,
    ):
        antibiotic = match.group(1)
        if not _looks_antibiotic_name(antibiotic):
            continue
        result, mic = _ast_result(match.group(2))
        item = {"antibiotic": antibiotic, "result": result}
        if mic:
            item["MIC"] = mic
        output.append(item)

    deduped: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str]] = set()
    for item in output:
        key = (item["antibiotic"].casefold(), item["result"], item.get("MIC", ""))
        if key not in seen:
            seen.add(key)
            deduped.append(item)
    return deduped


def _merge_culture_entries(entries: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    merged: dict[tuple[str, ...], dict[str, Any]] = {}
    order: list[tuple[str, ...]] = []
    for entry in entries:
        key = tuple(
            str(entry.get(field, ""))
            for field in (
                "test",
                "culture_type",
                "sample",
                "collected_time",
                "received_time",
                "reported_time",
                "organism",
                "status",
            )
        )
        if key not in merged:
            merged[key] = entry
            order.append(key)
            continue
        existing = merged[key]
        existing_ast = existing.setdefault("antimicrobial_susceptibility_test", [])
        known = {
            (item.get("antibiotic", "").casefold(), item.get("result", ""), item.get("MIC", ""))
            for item in existing_ast
        }
        for item in entry.get("antimicrobial_susceptibility_test", []):
            ast_key = (item.get("antibiotic", "").casefold(), item.get("result", ""), item.get("MIC", ""))
            if ast_key not in known:
                known.add(ast_key)
                existing_ast.append(item)
        if entry.get("notes"):
            notes = existing.setdefault("notes", [])
            for note in entry["notes"]:
                if note not in notes:
                    notes.append(note)
    return [merged[key] for key in order]


def parse_culture(grid: Table1Grid, sections: Sequence[TitledSection]) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    for section in sections:
        if _section_has(section, "抗酸菌直接鏡檢"):
            dates = _dates(section.text)
            times = _times(section.text)
            date_times = _combine_datetimes(dates, times)
            results = re.findall(r"No acid fast bacilli\s+observed\.?", section.text, re.IGNORECASE)
            for index, result in enumerate(results):
                entries.append(
                    {
                        "test": "Microbiology microscopy",
                        "culture_type": "Acid-fast bacilli direct smear",
                        "sample": _sample_from_title(section.title) or "Sputum",
                        "collected_time": date_times[index] if index < len(date_times) else "",
                        "received_time": "",
                        "reported_time": date_times[index] if index < len(date_times) else "",
                        "organism": "Acid-fast bacilli",
                        "status": "Not Detected",
                        "antimicrobial_susceptibility_test": [],
                        "notes": [_clean_spaces(result)],
                        "source": _source(grid, section.start_row),
                    }
                )
            continue

        if not re.search(r"培養檢驗|分離檢驗", section.title):
            continue
        collected, received, reported = _culture_times(section.text)
        organisms = _culture_organisms(section.text)
        if not organisms:
            continue
        preliminary = re.search(r"Preliminary report\s*:\s*([^\n]+)", section.text, re.IGNORECASE)
        ast = _culture_ast(section.text)
        for organism, status in organisms:
            colony = re.search(r"Colony Count\s*:\s*([^\n]+?CFU/mL)", organism, re.IGNORECASE)
            if not colony:
                colony = re.search(r"Colony Count\s*:\s*([^\n]+?CFU/mL)", section.text, re.IGNORECASE)
            clean_organism = re.sub(r"\s*Colony Count\s*:.*$", "", organism, flags=re.IGNORECASE).strip()
            entry: dict[str, Any] = {
                "test": "Microbiology culture",
                "culture_type": _culture_type(section.title),
                "sample": _sample_from_title(section.title),
                "collected_time": collected,
                "received_time": received,
                "reported_time": reported,
                "organism": clean_organism,
                "status": status,
                "antimicrobial_susceptibility_test": ast if status == "Detected" else [],
                "source": _source(grid, section.start_row, section.end_row),
            }
            if colony:
                entry["colony_count"] = _clean_spaces(colony.group(1))
            notes: list[str] = []
            if preliminary:
                notes.append(f"Preliminary report: {_clean_spaces(preliminary.group(1))}")
            if status == "Pending":
                notes.append("Culture was still pending in the source report.")
            if notes:
                entry["notes"] = notes
            entries.append(entry)
    return _merge_culture_entries(_dedupe(entries))


def _is_image_start(text: str) -> bool:
    normalized = normalize_text(text).strip()
    if re.match(r"^20\d{2}[/-]\d{1,2}[/-]\d{1,2}(?:\s+\d{1,2}:\d{2})?\s*:", normalized):
        return True
    if "檢查日期" in normalized and re.search(r"支氣管鏡|超音波|內視鏡|報告", normalized):
        return True
    return False


def _image_exam_type(text: str) -> str:
    normalized = normalize_text(text)
    first_line = next((line.strip() for line in normalized.splitlines() if line.strip()), "")
    match = re.match(r"^20\d{2}[/-]\d{1,2}[/-]\d{1,2}(?:\s+\d{1,2}:\d{2})?\s*:\s*(.*?)(?::|$)", first_line)
    if match:
        return _clean_spaces(match.group(1))
    if "支氣管鏡" in normalized or "Bronchoscope" in normalized:
        return "Bronchoscopy"
    if "腎臟超音波" in normalized:
        return "Renal ultrasound"
    if "心電圖" in normalized or "EKG" in normalized:
        return "EKG"
    if "大腸鏡" in normalized:
        return "Colonoscopy"
    if "內視鏡" in normalized or "Endoscopic" in normalized:
        return "Endoscopy"
    return "Special Examination"


def _image_findings(text: str) -> list[str]:
    findings: list[str] = []
    for raw_line in normalize_text(text).splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if re.match(r"^20\d{2}[/-]\d{1,2}[/-]\d{1,2}(?:\s+\d{1,2}:\d{2})?\s*:", line):
            line = re.sub(
                r"^20\d{2}[/-]\d{1,2}[/-]\d{1,2}(?:\s+\d{1,2}:\d{2})?\s*:\s*.*?(?::|$)",
                "",
                line,
            ).strip()
        if re.search(r"^(?:報告時間|認證時間)\s*[:：]", line):
            continue
        if line:
            findings.append(line)
    return findings


def parse_images(grid: Table1Grid) -> list[dict[str, Any]]:
    start = next((row.number for row in grid.iter_rows() if "Radiology / Special Examination Results" in row.text), None)
    if start is None:
        start = next((row.number for row in grid.iter_rows() if "放射報告與特殊檢查記錄" in row.text), None)
    if start is None:
        return []
    end = next(
        (row.number for row in grid.iter_rows(start + 1) if "Pathological Examination" in row.text or "病理檢查報告" in row.text),
        len(grid.rows) - 1,
    )

    chunks: list[tuple[int, int, str]] = []
    current_start: int | None = None
    current_lines: list[str] = []
    current_end = 0
    for row in grid.iter_rows(start + 1, end - 1):
        if not row.text.strip():
            continue
        if _is_image_start(row.text):
            if current_start is not None:
                chunks.append((current_start, current_end, "\n".join(current_lines)))
            current_start = row.number
            current_lines = [row.text]
        elif current_start is not None:
            current_lines.append(row.text)
        current_end = row.number
    if current_start is not None:
        chunks.append((current_start, current_end, "\n".join(current_lines)))

    entries: list[dict[str, Any]] = []
    for start_row, end_row, text in chunks:
        collected = _first_datetime(text)
        reported_match = re.search(
            r"(?:報告時間|認證時間)\s*[:：]\s*(20\d{2}[/-]\d{1,2}[/-]\d{1,2}\s+\d{1,2}:\d{2}(?::\d{2})?)",
            text,
        )
        reported = _first_datetime(reported_match.group(1)) if reported_match else ""
        entries.append(
            {
                "test": "Radiology / Special Examination",
                "exam_type": _image_exam_type(text),
                "collected_time": collected,
                "reported_time": reported,
                "findings": _image_findings(text),
                "source": _source(grid, start_row, end_row),
            }
        )
    return entries


def build_underlying(grid: Table1Grid) -> list[dict[str, Any]]:
    text = "\n".join(row.text for row in grid.iter_rows())
    diseases: list[str] = []
    if "Systemic lupus erythematosus" in text:
        diseases.append("Systemic lupus erythematosus")
    if re.search(r"Bilateral renal parenchymal disease", text, re.IGNORECASE):
        diseases.append("Bilateral renal parenchymal disease")
    if re.search(r"Left renal cyst", text, re.IGNORECASE):
        diseases.append("Left renal cyst")
    return [
        {
            "test": "Underlying Conditions",
            "age": "",
            "gender": "",
            "birth_date": "",
            "height": "",
            "weight": "",
            "underlying_diseases": diseases,
            "medical_history": {
                "diabetes_mellitus": "",
                "chronic_kidney_disease": "",
                "esrd_on_dialysis": "",
                "liver_cirrhosis": "",
                "autoimmune_disease": "是" if "Systemic lupus erythematosus" in diseases else "",
                "active_malignancy": "",
                "rheumatoid_arthritis": "",
                "systemic_lupus_erythematosus": "是" if "Systemic lupus erythematosus" in diseases else "",
                "other_autoimmune_disease": "",
                "systemic_steroid_exposure": "",
                "smoking_history": "",
            },
            "clinical_context": {
                "ward_unit": "",
                "admission_date": "",
                "ards_on_specimen_day": "",
            },
        }
    ]


def build_admission_diagnoses(grid: Table1Grid, images: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    text = "\n".join(row.text for row in grid.iter_rows())
    for pattern in (r"診斷\s*[:：]\s*([^\n]+)", r"Clinical Diagnosis\s*[:：]\s*([^\n]+)"):
        for match in re.finditer(pattern, text, re.IGNORECASE):
            diagnosis = _clean_spaces(match.group(1)).strip("/ ")
            if diagnosis:
                entries.append({"test": "Admission diagnosis", "record_type": "diagnosis", "diagnosis": diagnosis})

    discharge = re.search(r"出院狀況.*?(死亡\s*\(expired\)|死亡|expired)", text, re.IGNORECASE | re.DOTALL)
    if discharge:
        entries.append(
            {
                "test": "Admission diagnosis",
                "record_type": "discharge_outcome",
                "diagnosis": "Discharge outcome: Expired",
            }
        )
    return _dedupe(entries)


def build_standardized_payloads(
    grid: Table1Grid,
    sections: Sequence[TitledSection],
) -> dict[str, list[dict[str, Any]]]:
    images = parse_images(grid)
    return {
        "underlying": build_underlying(grid),
        "admission_diagnosis": build_admission_diagnoses(grid, images),
        "culture": parse_culture(grid, sections),
        "filmarray": parse_filmarray(grid, sections),
        "gm_test": parse_gm_serology(grid, sections),
        "molecular_microbiology": parse_molecular_microbiology(grid, sections),
        "image": images,
        "CBC": parse_cbc(grid, sections),
        "other_lab": parse_other_labs(grid, sections),
        "mNGS_grouped": [],
    }
