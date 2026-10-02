from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, MutableMapping

CATEGORY_KEY_MAP = {
    "細菌": "bacterial",
    "病毒": "viral",
    "真菌": "fungal",
}
GROUP_KEYS = ("bacterial", "viral", "fungal", "others")


def _load_json(path: Path | str) -> Mapping[str, Any]:
    with Path(path).expanduser().resolve().open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _group_pathogens(suspected: Iterable[Mapping[str, Any]]) -> MutableMapping[str, list[dict[str, Any]]]:
    grouped: MutableMapping[str, list[dict[str, Any]]] = {key: [] for key in GROUP_KEYS}
    for organism in suspected or []:
        english_key = CATEGORY_KEY_MAP.get(str(organism.get("category")).strip(), "others")
        grouped[english_key].append(
            {
                "name": organism.get("name", ""),
                "reads": organism.get("count", 0),
            }
        )
    return grouped


def transform_mngs_payload(data: Mapping[str, Any]) -> dict[str, list[dict[str, Any]]]:
    transformed: dict[str, list[dict[str, Any]]] = {}
    for patient_id, entries in data.items():
        if not isinstance(entries, list):
            raise TypeError(f"Expected list for patient {patient_id}, got {type(entries)}")
        transformed_entries: list[dict[str, Any]] = []
        for entry in entries:
            # Accept either raw suspected_pathogens or already grouped pathogens.
            pathogens = entry.get("pathogens")
            if isinstance(pathogens, Mapping) and all(key in GROUP_KEYS for key in pathogens.keys()):
                grouped_pathogens = {k: list(v) for k, v in pathogens.items()}
            else:
                grouped_pathogens = _group_pathogens(entry.get("suspected_pathogens", []))

            specimen_site = entry.get("specimen_site")

            transformed_entry = {
                "specimen_code": entry.get("specimen_code", ""),
                "collected_time": entry.get("collected_time", ""),
            }
            # Keep specimen_site just after collected_time when present.
            if specimen_site:
                transformed_entry["specimen_site"] = specimen_site
            transformed_entry["pathogens"] = grouped_pathogens

            transformed_entries.append(transformed_entry)
        transformed[patient_id] = transformed_entries
    return transformed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Group mNGS suspected pathogens by category (bacterial / viral / fungal / others)."
    )
    parser.add_argument("input_json", type=Path, help="Path to the JSON produced by the raw mNGS parser.")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Where to store the transformed JSON (default: outputs/<input_stem>_mNGS_grouped.json).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    payload = _load_json(args.input_json)
    converted = transform_mngs_payload(payload)

    output_path = args.output
    if not output_path:
        output_dir = Path("outputs")
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{args.input_json.stem}_mNGS_grouped.json"
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(converted, handle, ensure_ascii=False, indent=2)


if __name__ == "__main__":
    main()
