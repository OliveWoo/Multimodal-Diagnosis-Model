from __future__ import annotations

import argparse
import json
from pathlib import Path

from rag_re_casefit.batch import run_batch


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--merge-dir", type=Path, required=True)
    parser.add_argument("--frozen-rag-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=Path("config/default.json"))
    parser.add_argument("--rag-re-config", type=Path, required=True)
    parser.add_argument("--rag-re-cache-dir", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8-sig"))
    rag_re_config = json.loads(args.rag_re_config.read_text(encoding="utf-8-sig"))
    manifest = run_batch(
        merge_dir=args.merge_dir,
        frozen_rag_dir=args.frozen_rag_dir,
        output_dir=args.output_dir,
        config=config,
        cache_dir=args.cache_dir,
        rag_re_config=rag_re_config,
        rag_re_cache_dir=args.rag_re_cache_dir,
    )
    print(json.dumps({
        "output_count": manifest["output_count"],
        "error_count": manifest["error_count"],
        "manifest": str((args.output_dir / "manifest.json").resolve()),
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
