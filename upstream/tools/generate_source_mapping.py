from pathlib import Path
import json
import argparse


def build_source_mapping(base_dir: Path):
    mapping = {
        "CBC/otherLab": [],
        "filmarray": [],
        "culture": [],
        "image": [],
        "admission_diagnosis": [],
        "underlying": []
    }

    for file in base_dir.glob("*.json"):
        name = file.name.lower()

        if "cbc" in name or "other_lab" in name:
            mapping["CBC/otherLab"].append(file.name)
        elif "filmarray" in name or "gm_test" in name:
            mapping["filmarray"].append(file.name)
        elif "culture" in name:
            mapping["culture"].append(file.name)
        elif "image" in name or "radiology" in name:
            mapping["image"].append(file.name)
        elif "admission" in name:
            mapping["admission_diagnosis"].append(file.name)
        elif "underlying" in name:
            mapping["underlying"].append(file.name)

    # 移除空模組（沒有資料的）
    mapping = {k: v for k, v in mapping.items() if v}
    return mapping


def main():
    parser = argparse.ArgumentParser(
        description="Build a source mapping JSON file for a patient folder."
    )
    parser.add_argument(
        "folder",
        type=Path,
        help="Path to the patient folder containing JSON files."
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional path for the output JSON file. Defaults to <folder>/source_mapping.json."
    )
    args = parser.parse_args()

    folder = args.folder.expanduser().resolve()
    if not folder.exists():
        raise FileNotFoundError(f"Folder not found: {folder}")

    mapping = build_source_mapping(folder)
    output_path = (
        args.output.expanduser().resolve() if args.output else folder / "source_mapping.json"
    )

    with output_path.open("w", encoding="utf-8") as fh:
        json.dump(mapping, fh, indent=2, ensure_ascii=False)

    print(f"Source mapping written to {output_path}")


if __name__ == "__main__":
    main()
