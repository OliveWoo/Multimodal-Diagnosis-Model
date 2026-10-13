"""Validate immutable policy bundles before comparative pipeline runs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BUNDLE = ROOT / "rules/current_pipeline_policy_bundle_v6_v9_calendar_v7_frozen_20261013.json"


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate(bundle_path: Path = DEFAULT_BUNDLE) -> dict[str, Any]:
    bundle_path = bundle_path.resolve()
    bundle = read_json(bundle_path)
    if bundle.get("status") != "frozen_reference":
        raise ValueError(f"Bundle is not frozen: {bundle_path}")
    checks: list[dict[str, Any]] = []
    for component in bundle.get("components") or []:
        relative = Path(str(component.get("path") or ""))
        path = relative if relative.is_absolute() else ROOT / relative
        expected_hash = str(component.get("sha256") or "").casefold()
        exists = path.is_file()
        actual_hash = sha256_file(path) if exists else ""
        payload = read_json(path) if exists else {}
        expected_policy = str(component.get("policy_id") or "")
        actual_policy = str(payload.get("policy_id") or "")
        passed = bool(
            exists
            and expected_hash
            and actual_hash == expected_hash
            and expected_policy
            and actual_policy == expected_policy
        )
        checks.append(
            {
                "role": component.get("role"),
                "path": str(path),
                "exists": exists,
                "expected_policy_id": expected_policy,
                "actual_policy_id": actual_policy,
                "expected_sha256": expected_hash,
                "actual_sha256": actual_hash,
                "passed": passed,
            }
        )
    if not checks:
        raise ValueError(f"Bundle has no components: {bundle_path}")
    failed = [row for row in checks if not row["passed"]]
    result = {
        "schema_version": "frozen_pipeline_policy_bundle_validation.v1",
        "bundle_id": bundle.get("bundle_id"),
        "bundle_path": str(bundle_path),
        "bundle_sha256": sha256_file(bundle_path),
        "component_count": len(checks),
        "passed_count": len(checks) - len(failed),
        "failed_count": len(failed),
        "passed": not failed,
        "checks": checks,
    }
    if failed:
        raise ValueError(json.dumps(result, ensure_ascii=False, indent=2))
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bundle", type=Path, default=DEFAULT_BUNDLE)
    args = parser.parse_args()
    print(json.dumps(validate(args.bundle), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
