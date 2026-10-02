from __future__ import annotations

import unittest
from unittest.mock import patch

from tools import codex_cli_auth


class CodexCliAuthTests(unittest.TestCase):
    @patch("tools.codex_cli_auth.subprocess.run")
    @patch("tools.codex_cli_auth.find_codex_executable", return_value=r"C:\Codex\codex.exe")
    def test_login_uses_bundled_executable(self, _find, run) -> None:
        run.return_value.returncode = 0
        self.assertEqual(codex_cli_auth.run("login"), 0)
        self.assertEqual(run.call_args.args[0], [r"C:\Codex\codex.exe", "login", "--device-auth"])

    @patch("tools.codex_cli_auth.subprocess.run")
    @patch("tools.codex_cli_auth.find_codex_executable", return_value=r"C:\Codex\codex.exe")
    def test_status_uses_bundled_executable(self, _find, run) -> None:
        run.return_value.returncode = 1
        self.assertEqual(codex_cli_auth.run("status"), 1)
        self.assertEqual(run.call_args.args[0], [r"C:\Codex\codex.exe", "login", "status"])


if __name__ == "__main__":
    unittest.main()
