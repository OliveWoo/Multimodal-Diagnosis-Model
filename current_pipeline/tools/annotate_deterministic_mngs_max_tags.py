from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


DEFAULT_INPUT_SUFFIX = "mNGS_max_deterministic_opt_chosen_full"
DEFAULT_OUTPUT_SUFFIX = "mNGS_max_deterministic_opt_chosen_full_tag_explained"


VALUE_DESCRIPTIONS: dict[str, dict[Any, str]] = {
    "host_vulnerability_tier": {
        "V0": "未見明確免疫低下或重要宿主脆弱因子",
        "V1": "輕度宿主風險",
        "V2": "中度免疫低下或重要共病",
        "V3": "高度免疫低下",
        "Unknown": "資料不足",
    },
    "opportunistic_coverage_level": {
        "O0": "未見明確機會感染風險",
        "O1": "需要考慮機會感染",
        "O2": "高度需要考慮機會感染",
        "Unknown": "資料不足",
    },
    "host_level3_expansion": {
        True: "宿主狀態允許較寬鬆保留 Level 3 候選；不代表確診",
        False: "不因宿主狀態擴大保留 Level 3",
    },
    "mngs_signal_tier": {
        "M1_strong": "強 mNGS 訊號",
        "M2_moderate": "中等 mNGS 訊號",
        "M3_protected": "重要或較難檢出的病原被保留；不等於確診",
        "M4_weak": "弱 mNGS 訊號",
        "M5_background": "背景、定殖、污染或再活化傾向",
        "Unknown": "無法判定",
    },
    "reads_tier": {
        "R4_very_high": "reads 非常高",
        "R3_high": "reads 高",
        "R2_medium": "reads 中等",
        "R1_low": "reads 低",
        "R0_trace": "trace reads，訊號極低",
    },
    "rank_priority": {
        1: "背景控制條件最有利",
        2: "背景控制條件次佳",
        3: "背景控制條件中等",
        4: "背景控制條件較弱",
        5: "背景控制條件最弱",
        "1": "背景控制條件最有利",
        "2": "背景控制條件次佳",
        "3": "背景控制條件中等",
        "4": "背景控制條件較弱",
        "5": "背景控制條件最弱",
    },
    "integrated_causative_level": {
        "Level 1": "最強支持",
        "Level 2": "高可能性",
        "Level 3": "保留候選；需搭配臨床與其他證據",
        "Level 4": "弱支持或 fallback best available",
        "Level 5": "背景、定殖、污染或非主要病原傾向",
        "Unknown": "無法分級",
    },
    "basis_level": {
        "Level 1": "最強支持",
        "Level 2": "高可能性",
        "Level 3": "保留候選；需搭配臨床與其他證據",
        "Level 4": "弱支持或 fallback best available",
        "Level 5": "背景、定殖、污染或非主要病原傾向",
        "Unknown": "無法分級",
    },
    "observed_level": {
        "Level 1": "最強支持",
        "Level 2": "高可能性",
        "Level 3": "保留候選；需搭配臨床與其他證據",
        "Level 4": "弱支持或 fallback best available",
        "Level 5": "背景、定殖、污染或非主要病原傾向",
        "Unknown": "無法分級",
    },
    "specimen_alignment": {
        "Aligned": "檢體來源與推定感染來源相符",
        "Sterile_or_Systemic": "血液、無菌部位或全身性來源",
        "Non_aligned": "檢體來源與推定感染來源不相符",
        "Unknown": "無法判定檢體對位",
    },
    "is_protected_pathogen": {
        True: "protected pathogen；重要或較難檢出的病原，低 reads 時不直接排除",
        False: "未命中 protected pathogen 保留名單",
    },
    "protected_retention_condition": {
        True: "符合 protected pathogen 保留條件",
        False: "未符合 protected pathogen 保留條件",
    },
    "is_likely_colonizer_or_background": {
        True: "較像背景菌、定殖菌、污染或再活化訊號",
        False: "未命中背景或定殖 guardrail",
    },
    "selection_mode": {
        "Confirmed": "有 Level 1-3 或符合規則的 Best Available 候選",
        "Fallback": "無 Level 1-3 時，依規則輸出 Level 4/5 作為 best available",
        "No_high_priority_candidate": "沒有可輸出的候選",
    },
    "picked_role": {
        "Primary": "主要候選",
        "Secondary": "次要候選",
        "Protected_Level3": "因 protected/host 規則保留的 Level 3 候選",
        "BestAvailable": "無更高等級候選時的 best available",
    },
    "final_infection_likelihood": {
        "Severe": "高度支持感染或重症感染情境",
        "Likely": "較可能有感染相關病原",
        "Possible": "可能感染，但證據有限",
        "None": "未見足夠感染病原支持",
        "Unknown": "資料不足",
    },
    "dominant_source": {
        "Respiratory": "呼吸道來源",
        "Systemic": "全身性或血流來源",
        "Mixed": "多來源",
        "Unknown": "來源不明",
    },
    "dominant_pathogen_type": {
        "Bacterial": "細菌為主",
        "Viral": "病毒為主",
        "Fungal": "真菌為主",
        "Parasitic": "寄生蟲為主",
        "Mixed": "多種類型病原",
        "Unknown": "病原類型不明",
    },
}

CAUTION_FLAG_DESCRIPTIONS = {
    "Evidence_limited": "證據有限，需搭配臨床判斷",
    "Not_recommended_as_sole_treatment_basis": "不建議只依此結果作為治療依據",
    "Possible_reactivation": "可能是病毒再活化，不一定是主要感染原因",
    "Host_expansion_applied": "因宿主狀態放寬保留 Level 3",
    "Protected_low_read_pathogen": "低 reads 但因 protected pathogen 規則被保留",
}

RULE_DESCRIPTIONS = {
    "D-S1-CANDIDATE_FROM_RANKED": "候選病原來自 mNGS ranked candidates",
    "R-S3-M1": "mNGS 訊號判為 M1_strong：rank 1-2 且 reads 高，或 reads percentile >=0.85，或 dominant 且 rank <=3",
    "R-S3-M2": "mNGS 訊號判為 M2_moderate：中等訊號；reads_percentile 規則不適用於 R0_trace",
    "R-S3-M3": "mNGS 訊號判為 M3_protected：protected pathogen，或 culture/FilmArray/image 等非宿主輔助證據保留",
    "R-S3-M4": "mNGS 訊號判為 M4_weak：弱訊號",
    "R-S3-M5": "mNGS 訊號判為 M5_background：背景、定殖、污染或再活化傾向",
    "R-S4-BACKGROUND": "命中背景菌、定殖菌、污染或再活化 guardrail",
    "R-S4-CANDIDA-RESP-L3": "呼吸道 Candida 條件性保留到 Level 3；不代表確診",
    "R-S4-ORAL-UPPER_AIRWAY_HARD": "口咽或上呼吸道共生菌硬性 guardrail",
    "R-S5-L1": "整合分級為 Level 1：M1_strong 且檢體對位為 Aligned/Sterile_or_Systemic，且為 dominant 或 reads_tier 為 R4/R3",
    "R-S5-L2": "整合分級為 Level 2：M1_strong 但未達 Level 1 條件，或 M2_moderate 且檢體對位為 Aligned/Sterile_or_Systemic",
    "R-S5-L3": "整合分級為 Level 3：M2_moderate 但檢體未對位，或 M3_protected，或免疫低下時保留仍合理的 M4_weak",
    "R-S5-L4": "整合分級為 Level 4：M4_weak 且未被宿主條件提升，或口咽/上呼吸道共生菌 guardrail",
    "R-S5-L5": "整合分級為 Level 5：M5_background 或命中背景、定殖、污染、再活化 guardrail",
}


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Expand deterministic mNGS max tags into human-readable Traditional Chinese JSON."
    )
    parser.add_argument(
        "patient_root",
        nargs="?",
        type=Path,
        default=Path("outputs") / "patient_info_rich_normalized",
    )
    parser.add_argument("--patients", nargs="*", help="Optional patient IDs to process.")
    parser.add_argument("--input-json", type=Path, help="Annotate one deterministic max JSON file.")
    parser.add_argument("--output", type=Path, help="Output JSON path for --input-json mode.")
    parser.add_argument("--input-suffix", default=DEFAULT_INPUT_SUFFIX)
    parser.add_argument("--output-suffix", default=DEFAULT_OUTPUT_SUFFIX)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--skip-existing", action="store_true")
    return parser.parse_args(argv)


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as handle:
        return json.load(handle)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


def extract_patient_id(patient_dir: Path) -> str:
    parts = patient_dir.name.split("_")
    if len(parts) >= 3 and parts[0] == "NGS" and parts[1] == "patient":
        return parts[2]
    raise ValueError(f"Cannot extract patient ID from {patient_dir}")


def patient_dirs(patient_root: Path, patients: Sequence[str] | None) -> list[Path]:
    wanted = {str(pid) for pid in patients or []}
    dirs = sorted(
        (path for path in patient_root.glob("NGS_patient_*_json") if path.is_dir()),
        key=lambda path: int(extract_patient_id(path)),
    )
    if wanted:
        dirs = [path for path in dirs if extract_patient_id(path) in wanted]
    return dirs


def suffix_paths(patient_dir: Path, input_suffix: str, output_suffix: str) -> tuple[Path, Path]:
    patient_id = extract_patient_id(patient_dir)
    summary_dir = patient_dir / "summary_outputs"
    return (
        summary_dir / f"NGS_patient_{patient_id}_{input_suffix}.json",
        summary_dir / f"NGS_patient_{patient_id}_{output_suffix}.json",
    )


def explain_value(key: str, value: Any) -> Any:
    if key == "caution_flags" and isinstance(value, list):
        annotated_flags: list[Any] = []
        for item in value:
            if not isinstance(item, str):
                annotated_flags.append(item)
                continue
            description = CAUTION_FLAG_DESCRIPTIONS.get(item)
            annotated_flags.append(f"{item} {description}" if description else item)
        return annotated_flags
    if key == "applied_rules" and isinstance(value, list):
        annotated_rules: list[Any] = []
        for item in value:
            if not isinstance(item, str):
                annotated_rules.append(item)
                continue
            description = RULE_DESCRIPTIONS.get(item)
            annotated_rules.append(f"{item} {description}" if description else item)
        return annotated_rules
    if key in {"level_rule", "guardrail_rule"} and isinstance(value, str):
        description = RULE_DESCRIPTIONS.get(value)
        return f"{value} {description}" if description else value
    descriptions = VALUE_DESCRIPTIONS.get(key)
    if descriptions is None:
        return value
    if value in descriptions:
        display_value = str(value).lower() if isinstance(value, bool) else value
        return f"{display_value} {descriptions[value]}"
    value_text = str(value)
    if value_text in descriptions:
        return f"{value_text} {descriptions[value_text]}"
    return value


def annotate(obj: Any, parent_key: str = "") -> Any:
    if isinstance(obj, dict):
        return {key: annotate(explain_value(key, value), key) for key, value in obj.items()}
    if isinstance(obj, list):
        if parent_key in {"applied_rules", "caution_flags"}:
            return explain_value(parent_key, obj)
        return [annotate(item, parent_key) for item in obj]
    return obj


def annotate_file(
    input_path: Path,
    output_path: Path,
    *,
    overwrite: bool,
    skip_existing: bool,
) -> Path:
    if output_path.exists() and skip_existing:
        return output_path
    if output_path.exists() and not overwrite:
        raise FileExistsError(f"Output exists: {output_path}. Use --overwrite or --skip-existing.")
    if not input_path.exists():
        raise FileNotFoundError(input_path)
    payload = load_json(input_path)
    annotated = annotate(payload)
    if isinstance(annotated, dict):
        annotated.setdefault("_tag_explanation_metadata", {})
        annotated["_tag_explanation_metadata"] = {
            "source_json": str(input_path),
            "note": "此檔為人可讀 tag 展開版；不建議作為下游程式計算輸入。",
            "glossary_readme": "README_mngs_deterministic_tag_glossary.md",
        }
    write_json(output_path, annotated)
    return output_path


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    if args.input_json:
        output_path = args.output or args.input_json.with_name(f"{args.input_json.stem}_tag_explained.json")
        written = annotate_file(
            args.input_json,
            output_path,
            overwrite=args.overwrite,
            skip_existing=args.skip_existing,
        )
        print(f"wrote {written}")
        return

    written: list[Path] = []
    failures: list[str] = []
    for patient_dir in patient_dirs(args.patient_root, args.patients):
        input_path, output_path = suffix_paths(patient_dir, args.input_suffix, args.output_suffix)
        try:
            written.append(
                annotate_file(
                    input_path,
                    output_path,
                    overwrite=args.overwrite,
                    skip_existing=args.skip_existing,
                )
            )
            print(f"wrote {output_path}")
        except Exception as exc:  # noqa: BLE001
            failures.append(f"{input_path}: {exc}")
            print(f"failed {input_path}: {exc}", file=sys.stderr)
    print(f"written_count={len(written)}")
    if failures:
        print("failures:")
        for failure in failures:
            print(f"- {failure}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
