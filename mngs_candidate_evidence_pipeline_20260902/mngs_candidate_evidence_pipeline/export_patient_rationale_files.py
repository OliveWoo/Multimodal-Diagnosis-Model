#!/usr/bin/env python3
"""Split accepted rationale results into one readable Markdown file per patient.

Offline formatting only: no model calls or clinical text rewrites.
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from audit_llm_rationale_batch import link, render_patient
from generate_llm_rationales import load_json, write_json


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    source = load_json(args.input)
    if source.get("status") != "structurally_validated":
        raise SystemExit("Only structurally validated results may be exported")
    output_root = args.output_root.resolve()
    if output_root.exists() and any(output_root.iterdir()):
        raise SystemExit("Output folder must be new or empty")
    cases = source["patients"]
    seen = set()
    for entry in cases:
        cohort, patient_id = entry["cohort"], str(entry["result"]["patient_id"])
        if cohort not in ("KH", "two_hospital") or not patient_id.isdecimal():
            raise SystemExit("Invalid cohort or patient identifier")
        key = (cohort, patient_id)
        if key in seen:
            raise SystemExit(f"Duplicate patient: {key}")
        seen.add(key)
        if entry["result"] != load_json(Path(entry["output_file"])):
            raise SystemExit(f"Stored result differs from original output: {key}")
    output_root.mkdir(parents=True, exist_ok=True)
    index = [
        "# 每位病人一份：選菌理由", "",
        "每個 .md 檔只包含一位病人；按資料組分開，檔名即病例編號。", "",
        "本次只拆分既有輸出，未重跑 API、未更動選菌與理由文字。原本總檔與 JSON 仍保留。", "",
        "這是事後證據說明，不是新的診斷或原模型內部思考紀錄；引用與醫學敘述仍需人工核對。", "",
    ]
    manifest_rows = []
    counts = Counter()
    for cohort in ("KH", "two_hospital"):
        cohort_cases = sorted((entry for entry in cases if entry["cohort"] == cohort), key=lambda x: int(x["result"]["patient_id"]))
        index += [f"## {'高醫' if cohort == 'KH' else 'Two-hospital'}（{len(cohort_cases)} 位）", ""]
        for entry in cohort_cases:
            result = entry["result"]
            patient_id = str(result["patient_id"])
            generated = result["generated_rationale"]
            original_file = Path(entry["output_file"])
            output_file = output_root / cohort / f"P{patient_id}_reason_zh.md"
            body = render_patient(result, original_file)
            body[0] = body[0].replace("## ", "# ", 1)
            body = [line.replace("### ", "## ", 1) if line.startswith("### ") else line for line in body]
            notice = (
                "本頁整理既有選菌結果，不重新判定感染源；未經全面臨床審查。"
                "若未提供證據或檢驗編號，不代表醫院沒有該項資料。"
            )
            body[2:2] = [notice, ""]
            if generated["patient_level_evidence_ids"] and not generated["no_selected_pathogen"]:
                first_organism = next(i for i, line in enumerate(body) if line.startswith("## "))
                refs = ["病例總結引用（可包含不同候選菌的比較證據，不代表同種支持）：", ""]
                refs += [f"- {eid}：`{result['evidence_map'][eid]}`" for eid in generated["patient_level_evidence_ids"]]
                refs.append("")
                body[first_organism:first_organism] = refs
            output_file.parent.mkdir(parents=True, exist_ok=True)
            output_file.write_text("\n".join(body) + "\n", encoding="utf-8")
            rendered = output_file.read_text(encoding="utf-8")
            expected_texts = [generated["patient_result_summary"]]
            for item in generated["selected_pathogen_explanations"]:
                expected_texts += [item["organism"], item["concise_reason"], item["key_evidence_summary"], item["clinician_review_question"]]
                expected_texts += [lim["text"] for lim in item["limitations"]]
            if not all(text in rendered for text in expected_texts):
                raise SystemExit(f"Text lost during formatting: {cohort} P{patient_id}")
            if sum(line.startswith("# ") for line in rendered.splitlines()) != 1:
                raise SystemExit(f"Invalid patient heading count: {output_file}")
            organisms = "、".join(item["organism"] for item in sorted(generated["selected_pathogen_explanations"], key=lambda x: x["rank"])) or "無入選病原"
            index.append(f"- {link('P' + patient_id, output_file)}：{organisms}")
            counts[cohort] += 1
            manifest_rows.append({"cohort": cohort, "patient_id": patient_id, "selected_count": len(generated["selected_pathogen_explanations"]), "output_file": str(output_file), "original_json": str(original_file), "text_preserved": True})
        index.append("")
    (output_root / "README_zh.md").write_text("\n".join(index) + "\n", encoding="utf-8")
    write_json(output_root / "export_manifest.json", {"input": str(args.input.resolve()), "offline_only": True, "patient_count": len(manifest_rows), "cohort_counts": dict(counts), "patients": manifest_rows})
    print(f"Exported {len(manifest_rows)} patient files: {dict(counts)}; output={output_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
