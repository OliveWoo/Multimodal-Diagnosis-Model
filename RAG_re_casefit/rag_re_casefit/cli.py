from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List

from .case_card import build_case_card, find_candidate, validate_case_card
from .judge import CaseFitJudge
from .rules import aggregate_judgments


PROJECT = Path(__file__).resolve().parent.parent


def _read(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _write(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")
    temporary.replace(path)


def _normal_name(value: Any) -> str:
    return "".join(character for character in str(value or "").casefold() if character.isalnum())


def _articles(value: Any, organism: str | None = None) -> List[Dict[str, Any]]:
    rows = value.get("articles") if isinstance(value, dict) else value
    if rows is None and isinstance(value, dict) and isinstance(value.get("candidates"), list):
        target = _normal_name(organism)
        aliases = {
            "cmv": "humancytomegalovirus",
            "hsv1": "herpessimplexvirustype1",
        }
        target = aliases.get(target, target)
        matches = []
        for candidate in value["candidates"]:
            if not isinstance(candidate, dict):
                continue
            names = {
                _normal_name(candidate.get("organism_name")),
                _normal_name(candidate.get("canonical_organism_name")),
            }
            names |= {aliases.get(name, name) for name in names}
            if target in names:
                matches.append(candidate)
        if len(matches) != 1:
            raise ValueError(
                f"expected exactly one organism in RAG output, found {len(matches)}: {organism}"
            )
        evidence = matches[0].get("literature_evidence")
        retrieval = evidence.get("retrieval") if isinstance(evidence, dict) else None
        rows = retrieval.get("articles") if isinstance(retrieval, dict) else None
    if not isinstance(rows, list) or any(not isinstance(row, dict) for row in rows):
        raise ValueError(
            "articles input must be an array, an object with an articles array, "
            "or a RAG_re output containing the organism's frozen retrieval"
        )
    return rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="RAG_re parallel article-to-patient case-fit judge")
    sub = parser.add_subparsers(dest="command", required=True)

    card = sub.add_parser("build-card")
    card.add_argument("--input", type=Path, required=True)
    card.add_argument("--organism", required=True)
    card.add_argument("--output", type=Path, required=True)

    judge = sub.add_parser("judge")
    judge.add_argument("--case-card", type=Path, required=True)
    judge.add_argument("--organism", required=True)
    judge.add_argument("--articles", type=Path, required=True)
    judge.add_argument("--output", type=Path, required=True)
    judge.add_argument("--config", type=Path, default=PROJECT / "config" / "default.json")
    judge.add_argument("--cache-dir", type=Path)
    return parser


def main(argv: List[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "build-card":
        source_bytes = args.input.read_bytes()
        root = json.loads(source_bytes.decode("utf-8-sig"))
        candidate = find_candidate(root, args.organism)
        _write(args.output, build_case_card(root, candidate, source_bytes=source_bytes))
        return 0

    card = validate_case_card(_read(args.case_card))
    config = _read(args.config)
    llm_config = config.get("llm")
    if not isinstance(llm_config, dict):
        raise ValueError("config.llm must be an object")
    articles = _articles(_read(args.articles), args.organism)
    judge = CaseFitJudge(llm_config, cache_dir=args.cache_dir)
    judgments = judge.judge_many(args.organism, card, articles)
    output = {
        "schema_version": "rag_re_casefit.candidate_output.v1",
        "organism": args.organism,
        "case_card": card,
        "judgments": judgments,
        "aggregate": aggregate_judgments(judgments, config.get("aggregation")),
    }
    _write(args.output, output)
    return 0
