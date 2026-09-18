from __future__ import annotations

import argparse
import json
from pathlib import Path

from rag_re_casefit.evaluation import evaluate, render_markdown


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--gold", type=Path, required=True)
    parser.add_argument("--json-out", type=Path, required=True)
    parser.add_argument("--md-out", type=Path, required=True)
    args = parser.parse_args()
    result = evaluate(args.output_dir, args.gold)
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    args.md_out.write_text(render_markdown(result), encoding="utf-8")
    print(json.dumps({
        "candidate_count": result["scope"]["candidate_count"],
        "gold_candidate_count": result["scope"]["gold_candidate_count"],
        "json": str(args.json_out.resolve()),
        "markdown": str(args.md_out.resolve()),
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
