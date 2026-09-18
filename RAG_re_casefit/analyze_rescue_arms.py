from __future__ import annotations

import argparse
import json
from pathlib import Path

from rag_re_casefit.rescue_evaluation import evaluate_rescue_arms, render_rescue_markdown


def main() -> int:
    parser = argparse.ArgumentParser(description="Evaluate immutable-picked A1/B_STRICT/C3 rescue arms")
    parser.add_argument("--casefit-dir", required=True, type=Path)
    parser.add_argument("--merge-dir", required=True, type=Path)
    parser.add_argument("--gold", required=True, type=Path)
    parser.add_argument("--clinical-config", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    parser.add_argument("--out-md", required=True, type=Path)
    args = parser.parse_args()

    report = evaluate_rescue_arms(
        args.casefit_dir, args.merge_dir, args.gold, args.clinical_config
    )
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_md.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    args.out_md.write_text(render_rescue_markdown(report), encoding="utf-8")
    print(json.dumps({
        "schema_version": report["schema_version"],
        "out_json": str(args.out_json),
        "out_md": str(args.out_md),
        "scope": report["scope"],
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
