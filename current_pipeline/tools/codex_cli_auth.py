"""Locate the Codex desktop CLI and run its authentication commands."""

from __future__ import annotations

import argparse
import os
import subprocess
from collections.abc import Sequence

from tools.draft_organism_taxonomy_reviews import find_codex_executable


def codex_environment() -> dict[str, str]:
    environment = os.environ.copy()
    for variable in ("OPENAI_API_KEY", "OPENAI_APIKEY", "OPENAI_BASE_URL", "OPENAI_API_BASE"):
        environment.pop(variable, None)
    return environment


def run(command: str) -> int:
    executable = find_codex_executable()
    if not executable:
        raise RuntimeError(
            "Codex CLI was not found on PATH or in the Codex desktop installation."
        )
    print(f"Using Codex CLI: {executable}", flush=True)
    arguments = [executable, "login"]
    if command == "login":
        arguments.append("--device-auth")
    else:
        arguments.append("status")
    completed = subprocess.run(
        arguments,
        env=codex_environment(),
        check=False,
    )
    return completed.returncode


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Authenticate the Codex CLI bundled with the desktop app without requiring PATH setup."
    )
    parser.add_argument("command", choices=("login", "status"))
    args = parser.parse_args(argv)
    raise SystemExit(run(args.command))


if __name__ == "__main__":
    main()
