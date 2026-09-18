from __future__ import annotations

import argparse
import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence


TriState = bool | None


def tri_and(*values: TriState) -> TriState:
    if any(value is False for value in values):
        return False
    if all(value is True for value in values):
        return True
    return None


def tri_or(*values: TriState) -> TriState:
    if any(value is True for value in values):
        return True
    if all(value is False for value in values):
        return False
    return None


def normalize(value: Any) -> str:
    text = str(value or "").strip().lower()
    text = re.sub(r"\s+", " ", text)
    return text


BUILTIN_ALIASES = {
    "hsv": "hsv-1",
    "herpes simplex virus 1": "hsv-1",
    "human herpesvirus 1": "hsv-1",
    "pjp": "pneumocystis jirovecii",
    "pneumocystis jiroveci": "pneumocystis jirovecii",
    "candida albican": "candida albicans",
    "crkp": "klebsiella pneumoniae",
}


GENUS_STOPWORDS = {
    "human",
    "herpes",
    "herpesvirus",
    "influenza",
    "parainfluenza",
    "respiratory",
    "virus",
    "viral",
}


@dataclass(frozen=True)
class Formula:
    name: str
    expression: str
    # Disjunctive normal form: OR across terms, AND within each term.
    terms: tuple[tuple[str, ...], ...]

    @property
    def uses_a(self) -> bool:
        return any("A" in term for term in self.terms)

    @property
    def complexity(self) -> int:
        return sum(len(term) for term in self.terms) + len(self.terms) - 1

    def evaluate(self, values: Mapping[str, TriState]) -> TriState:
        return tri_or(*(tri_and(*(values[key] for key in term)) for term in self.terms))


FORMULAS = (
    Formula("A", "A", (("A",),)),
    Formula("B", "B", (("B",),)),
    Formula("C", "C", (("C",),)),
    Formula("A_OR_B", "A OR B", (("A",), ("B",))),
    Formula("A_OR_C", "A OR C", (("A",), ("C",))),
    Formula("B_OR_C", "B OR C", (("B",), ("C",))),
    Formula("A_OR_B_OR_C", "A OR B OR C", (("A",), ("B",), ("C",))),
    Formula("A_AND_B", "A AND B", (("A", "B"),)),
    Formula("A_AND_C", "A AND C", (("A", "C"),)),
    Formula("B_AND_C", "B AND C", (("B", "C"),)),
    Formula("A_AND_B_AND_C", "A AND B AND C", (("A", "B", "C"),)),
    Formula(
        "TWO_OF_THREE",
        "at least 2 of A/B/C",
        (("A", "B"), ("A", "C"), ("B", "C")),
    ),
    Formula("A_OR_B_AND_C", "A OR (B AND C)", (("A",), ("B", "C"))),
    Formula("B_OR_A_AND_C", "B OR (A AND C)", (("B",), ("A", "C"))),
    Formula("C_OR_A_AND_B", "C OR (A AND B)", (("C",), ("A", "B"))),
    Formula("A_AND_B_OR_C", "A AND (B OR C)", (("A", "B"), ("A", "C"))),
    Formula("B_AND_A_OR_C", "B AND (A OR C)", (("A", "B"), ("B", "C"))),
    Formula("C_AND_A_OR_B", "C AND (A OR B)", (("A", "C"), ("B", "C"))),
)

PRUNE_FORMULAS = (
    Formula("NO_PRUNE", "TRUE", ((),)),
    next(formula for formula in FORMULAS if formula.name == "B"),
    next(formula for formula in FORMULAS if formula.name == "C"),
    next(formula for formula in FORMULAS if formula.name == "B_OR_C"),
    next(formula for formula in FORMULAS if formula.name == "B_AND_C"),
)


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def canonical(value: Any, aliases: Mapping[str, str]) -> str:
    current = normalize(value)
    seen: set[str] = set()
    while current in aliases and current not in seen:
        seen.add(current)
        current = normalize(aliases[current])
    return current


def genus(value: str) -> str | None:
    match = re.match(r"^([a-z][a-z-]+)(?:\s+|$)", value)
    if not match or " " not in value:
        return None
    token = match.group(1)
    return None if token in GENUS_STOPWORDS else token


def names_match(left: str, right: str, *, genus_relaxed: bool) -> bool:
    if left == right:
        return True
    return bool(genus_relaxed and genus(left) and genus(left) == genus(right))


def maximum_match_count(
    truth: Sequence[str], predictions: Sequence[str], *, genus_relaxed: bool
) -> int:
    # One prediction may satisfy at most one gold label and vice versa.
    edges = [
        [index for index, predicted in enumerate(predictions) if names_match(gold, predicted, genus_relaxed=genus_relaxed)]
        for gold in truth
    ]
    assigned: dict[int, int] = {}

    def augment(gold_index: int, visited: set[int]) -> bool:
        for predicted_index in edges[gold_index]:
            if predicted_index in visited:
                continue
            visited.add(predicted_index)
            if predicted_index not in assigned or augment(assigned[predicted_index], visited):
                assigned[predicted_index] = gold_index
                return True
        return False

    return sum(augment(index, set()) for index in range(len(truth)))


def a_value(candidate: Mapping[str, Any], threshold: int) -> TriState:
    module = ((candidate.get("modules") or {}).get("A") or {})
    status = str(module.get("status") or "")
    if status != "ok":
        return None
    judgeable = int(module.get("judgeable_count") or 0)
    if judgeable < 5:
        return None
    return int(module.get("support_count") or 0) >= threshold


def abc_values(candidate: Mapping[str, Any], threshold: int) -> dict[str, TriState]:
    modules = candidate.get("modules") or {}
    return {
        "A": a_value(candidate, threshold),
        "B": (modules.get("B") or {}).get("positive"),
        "C": (modules.get("C") or {}).get("positive"),
    }


def fbeta(precision: float | None, recall: float | None, beta: float) -> float | None:
    if precision is None or recall is None:
        return None
    if precision == 0 and recall == 0:
        return 0.0
    beta_squared = beta * beta
    return (1 + beta_squared) * precision * recall / (beta_squared * precision + recall)


def evaluate_policy(
    *,
    outputs: Sequence[Mapping[str, Any]],
    gold_by_patient: Mapping[str, Sequence[str]],
    aliases: Mapping[str, str],
    selector: Callable[[Mapping[str, Any]], TriState],
    positive_patients_only: bool,
    genus_relaxed: bool,
) -> dict[str, Any]:
    tp = fp = fn = predictions_total = abstentions = 0
    patient_rows = []
    for output in outputs:
        patient_id = str(output.get("patient_id") or "")
        if patient_id not in gold_by_patient:
            continue
        truth = [canonical(value, aliases) for value in gold_by_patient[patient_id]]
        if positive_patients_only and not truth:
            continue
        selected_names: list[str] = []
        abstained_names: list[str] = []
        for candidate in output.get("candidates") or []:
            decision = selector(candidate)
            name = canonical(
                candidate.get("canonical_organism_name") or candidate.get("organism_name"),
                aliases,
            )
            if decision is True:
                selected_names.append(name)
            elif decision is None:
                abstained_names.append(name)
        selected_names = list(dict.fromkeys(selected_names))
        abstained_names = list(dict.fromkeys(abstained_names))
        matched = maximum_match_count(truth, selected_names, genus_relaxed=genus_relaxed)
        case_fp = len(selected_names) - matched
        case_fn = len(truth) - matched
        tp += matched
        fp += case_fp
        fn += case_fn
        predictions_total += len(selected_names)
        abstentions += len(abstained_names)
        patient_rows.append(
            {
                "patient_id": patient_id,
                "tp": matched,
                "fp": case_fp,
                "fn": case_fn,
                "predictions": len(selected_names),
            }
        )
    precision = tp / (tp + fp) if tp + fp else None
    recall = tp / (tp + fn) if tp + fn else None
    return {
        "patient_count": len(patient_rows),
        "truth_count": tp + fn,
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "predictions": predictions_total,
        "abstentions": abstentions,
        "precision": precision,
        "recall": recall,
        "f0_5": fbeta(precision, recall, 0.5),
        "f1": fbeta(precision, recall, 1.0),
        "f2": fbeta(precision, recall, 2.0),
    }


def selected_signature(
    outputs: Sequence[Mapping[str, Any]], selector: Callable[[Mapping[str, Any]], TriState]
) -> tuple[tuple[str, tuple[str, ...]], ...]:
    rows = []
    for output in outputs:
        selected = sorted(
            normalize(candidate.get("canonical_organism_name") or candidate.get("organism_name"))
            for candidate in output.get("candidates") or []
            if selector(candidate) is True
        )
        rows.append((str(output.get("patient_id") or ""), tuple(selected)))
    return tuple(rows)


def policy_rows(outputs: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    policies: list[dict[str, Any]] = []

    def add(
        *,
        family: str,
        name: str,
        expression: str,
        threshold: int | None,
        complexity: int,
        selector: Callable[[Mapping[str, Any]], TriState],
    ) -> None:
        policies.append(
            {
                "family": family,
                "name": name,
                "expression": expression,
                "a_threshold": threshold,
                "complexity": complexity,
                "selector": selector,
                "signature": selected_signature(outputs, selector),
            }
        )

    add(
        family="baseline",
        name="E0_BASELINE",
        expression="E0",
        threshold=None,
        complexity=1,
        selector=lambda candidate: bool(candidate.get("baseline_selected")),
    )

    for formula in FORMULAS:
        thresholds = range(3, 7) if formula.uses_a else (None,)
        for threshold in thresholds:
            effective_threshold = int(threshold or 3)
            suffix = f"__A_K{effective_threshold}" if formula.uses_a else ""

            def standalone(candidate: Mapping[str, Any], f: Formula = formula, k: int = effective_threshold) -> TriState:
                return f.evaluate(abc_values(candidate, k))

            add(
                family="standalone",
                name=f"{formula.name}{suffix}",
                expression=formula.expression,
                threshold=threshold,
                complexity=formula.complexity,
                selector=standalone,
            )

            def rescue(candidate: Mapping[str, Any], f: Formula = formula, k: int = effective_threshold) -> TriState:
                if candidate.get("baseline_selected") is True:
                    return True
                return f.evaluate(abc_values(candidate, k))

            add(
                family="rescue",
                name=f"E0_OR_{formula.name}{suffix}",
                expression=f"E0 OR ({formula.expression})",
                threshold=threshold,
                complexity=formula.complexity + 2,
                selector=rescue,
            )

    for prune in PRUNE_FORMULAS[1:]:
        def pruning(candidate: Mapping[str, Any], p: Formula = prune) -> TriState:
            if candidate.get("baseline_selected") is not True:
                return False
            return p.evaluate(abc_values(candidate, 3))

        add(
            family="pruning",
            name=f"E0_AND_{prune.name}",
            expression=f"E0 AND ({prune.expression})",
            threshold=None,
            complexity=prune.complexity + 2,
            selector=pruning,
        )

    # A is missing for all E0 candidates in this frozen evidence.  Therefore only
    # B/C-only filters are valid for the E0 side of a pruning or hybrid policy.
    for prune in PRUNE_FORMULAS:
        for rescue_formula in FORMULAS:
            thresholds = range(3, 7) if rescue_formula.uses_a else (None,)
            for threshold in thresholds:
                effective_threshold = int(threshold or 3)
                suffix = f"__A_K{effective_threshold}" if rescue_formula.uses_a else ""

                def hybrid(
                    candidate: Mapping[str, Any],
                    p: Formula = prune,
                    r: Formula = rescue_formula,
                    k: int = effective_threshold,
                ) -> TriState:
                    values = abc_values(candidate, k)
                    if candidate.get("baseline_selected") is True:
                        return p.evaluate(values)
                    return r.evaluate(values)

                add(
                    family="hybrid",
                    name=f"E0_{prune.name}__RESCUE_{rescue_formula.name}{suffix}",
                    expression=(
                        f"[E0 AND ({prune.expression})] OR "
                        f"[NOT E0 AND ({rescue_formula.expression})]"
                    ),
                    threshold=threshold,
                    complexity=prune.complexity + rescue_formula.complexity + 5,
                    selector=hybrid,
                )

    return policies


def pareto_frontier(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    frontier = []
    for row in rows:
        precision = row.get("precision")
        recall = row.get("recall")
        if precision is None or recall is None:
            continue
        dominated = any(
            other is not row
            and other.get("precision") is not None
            and other.get("recall") is not None
            and other["precision"] >= precision
            and other["recall"] >= recall
            and (other["precision"] > precision or other["recall"] > recall)
            for other in rows
        )
        if not dominated:
            frontier.append(dict(row))
    return sorted(frontier, key=lambda row: (-row["recall"], -row["precision"], row["complexity"]))


def percent(value: float | None) -> str:
    return "NA" if value is None else f"{value * 100:.1f}%"


def main() -> int:
    parser = argparse.ArgumentParser(description="Exhaustive monotone A/B/C factorial sweep")
    parser.add_argument("--outputs", type=Path, required=True)
    parser.add_argument("--gold", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()

    output_files = sorted(args.outputs.rglob("*.json")) if args.outputs.is_dir() else [args.outputs]
    outputs = [
        payload
        for path in output_files
        if isinstance((payload := load_json(path)), dict)
        and payload.get("schema_version") == "rag_re.output.v1"
    ]
    gold = load_json(args.gold)
    aliases = {normalize(key): normalize(value) for key, value in BUILTIN_ALIASES.items()}
    aliases.update(
        {normalize(key): normalize(value) for key, value in (gold.get("aliases") or {}).items()}
    )
    gold_by_patient = {
        str(case["patient_id"]): list(case.get("pathogens") or [])
        for case in gold.get("cases") or []
    }

    policies = policy_rows(outputs)
    scopes = {
        "upstream_answer_patients_genus_relaxed": {
            "positive_patients_only": True,
            "genus_relaxed": True,
        },
        "all_labeled_exact_synonym": {
            "positive_patients_only": False,
            "genus_relaxed": False,
        },
    }
    all_results: dict[str, list[dict[str, Any]]] = {}
    for scope_name, scope in scopes.items():
        rows = []
        for policy in policies:
            metrics = evaluate_policy(
                outputs=outputs,
                gold_by_patient=gold_by_patient,
                aliases=aliases,
                selector=policy["selector"],
                **scope,
            )
            rows.append(
                {
                    key: value
                    for key, value in policy.items()
                    if key not in {"selector", "signature"}
                }
                | metrics
                | {"signature": policy["signature"]}
            )

        # Keep the simplest label when multiple formulae produce exactly the same
        # patient-level predictions and metrics in the frozen dataset.
        deduped: dict[Any, dict[str, Any]] = {}
        for row in sorted(rows, key=lambda item: (item["complexity"], len(item["expression"]), item["name"])):
            signature = row.pop("signature")
            deduped.setdefault(signature, row)
        all_results[scope_name] = list(deduped.values())

    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    json_path = prefix.with_suffix(".json")
    csv_path = prefix.with_suffix(".csv")
    md_path = prefix.with_suffix(".md")

    json_payload = {
        "schema_version": "rag_re.factorial_sweep.v1",
        "exploratory_only": True,
        "warning": (
            "The clinical labels were used to rank policies. Results are hypothesis-generating, "
            "not an unbiased holdout estimate. A-containing E0 pruning was not evaluated because "
            "Module A was not attempted on E0 candidates."
        ),
        "evidence": {
            "outputs": str(args.outputs),
            "gold": str(args.gold),
            "output_count": len(outputs),
            "raw_policy_count": len(policies),
            "a_thresholds": [3, 4, 5, 6],
        },
        "scopes": {},
    }
    for scope_name, rows in all_results.items():
        json_payload["scopes"][scope_name] = {
            "unique_prediction_policy_count": len(rows),
            "pareto_frontier": pareto_frontier(rows),
            "policies": sorted(rows, key=lambda row: (-float(row["f0_5"] or -1), -float(row["recall"] or -1))),
        }
    with json_path.open("w", encoding="utf-8") as handle:
        json.dump(json_payload, handle, ensure_ascii=False, indent=2)

    fieldnames = [
        "scope", "family", "name", "expression", "a_threshold", "complexity",
        "patient_count", "truth_count", "tp", "fp", "fn", "predictions", "abstentions",
        "precision", "recall", "f0_5", "f1", "f2",
    ]
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for scope_name, rows in all_results.items():
            for row in rows:
                writer.writerow({"scope": scope_name, **{key: row.get(key) for key in fieldnames if key != "scope"}})

    upstream = all_results["upstream_answer_patients_genus_relaxed"]
    baseline = next(row for row in upstream if row["name"] == "E0_BASELINE")
    precision_focused = sorted(
        (row for row in upstream if row["recall"] is not None and row["recall"] >= 0.80),
        key=lambda row: (-row["precision"], -row["recall"], row["complexity"]),
    )[:12]
    moderate_floor = sorted(
        (row for row in upstream if row["recall"] is not None and row["recall"] >= 0.75),
        key=lambda row: (-row["precision"], -row["recall"], row["complexity"]),
    )[:12]
    frontier = pareto_frontier(upstream)

    lines = [
        "# A/B/C 全因子排列組合（探索性）",
        "",
        "## 重要限制",
        "",
        "- 這份資料的臨床答案已用來排序規則，因此只能用來產生假說，不能再宣稱為 untouched holdout 成效。",
        "- 凍結證據中，E0 候選沒有執行 Module A；因此 A 可用於非 E0 rescue，但不能公平評估 `E0 AND A` 類剪枝。",
        "- Primary view 採上游算法：只計有答案病人、凍結 synonym，加 genus-relaxed matching。",
        "",
        "## 搜尋空間",
        "",
        f"- 原始規則：{len(policies)} 組（A cutoff k=3/4/5/6）。",
        f"- 去除產生完全相同預測的重複規則後：{len(upstream)} 組。",
        "- 包含 standalone、E0-preserving rescue、E0 B/C pruning + non-E0 rescue hybrid。",
        "",
        "## E0 baseline",
        "",
        f"- TP/FP/FN = {baseline['tp']}/{baseline['fp']}/{baseline['fn']}；precision {percent(baseline['precision'])}；recall {percent(baseline['recall'])}。",
        "",
        "## Recall ≥ 80% 時，precision 最高",
        "",
        "| Rank | Policy | TP | FP | Precision | Recall | F0.5 |",
        "|---:|---|---:|---:|---:|---:|---:|",
    ]
    for index, row in enumerate(precision_focused, 1):
        lines.append(
            f"| {index} | `{row['expression']}`" + (f" (A k={row['a_threshold']})" if row.get("a_threshold") else "")
            + f" | {row['tp']} | {row['fp']} | {percent(row['precision'])} | {percent(row['recall'])} | {percent(row['f0_5'])} |"
        )
    lines.extend(
        [
            "",
            "## Recall ≥ 75% 時，precision 最高",
            "",
            "| Rank | Policy | TP | FP | Precision | Recall | F0.5 |",
            "|---:|---|---:|---:|---:|---:|---:|",
        ]
    )
    for index, row in enumerate(moderate_floor, 1):
        lines.append(
            f"| {index} | `{row['expression']}`" + (f" (A k={row['a_threshold']})" if row.get("a_threshold") else "")
            + f" | {row['tp']} | {row['fp']} | {percent(row['precision'])} | {percent(row['recall'])} | {percent(row['f0_5'])} |"
        )
    lines.extend(
        [
            "",
            "## Precision–recall Pareto frontier",
            "",
            "| Policy | TP | FP | Precision | Recall |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in frontier:
        lines.append(
            f"| `{row['expression']}`" + (f" (A k={row['a_threshold']})" if row.get("a_threshold") else "")
            + f" | {row['tp']} | {row['fp']} | {percent(row['precision'])} | {percent(row['recall'])} |"
        )
    lines.extend(
        [
            "",
            "## Files",
            "",
            f"- Full machine-readable results: `{json_path.name}`",
            f"- Full flat table: `{csv_path.name}`",
        ]
    )
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps({"json": str(json_path), "csv": str(csv_path), "markdown": str(md_path), "raw_policies": len(policies), "unique_upstream": len(upstream)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
