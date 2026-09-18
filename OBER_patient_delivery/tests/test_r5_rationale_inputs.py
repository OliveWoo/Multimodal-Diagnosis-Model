import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


MODULE_FILE = Path(__file__).resolve().parents[1] / "prepare_r5_rationale_inputs.py"
SPEC = importlib.util.spec_from_file_location("prepare_r5_rationale_inputs", MODULE_FILE)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader
SPEC.loader.exec_module(MODULE)


class R5RationaleInputTests(unittest.TestCase):
    def test_apply_final_decision_adds_rescue_without_mutating_source(self):
        chain = {
        "病例識別": {"病人編號": "1"},
        "一眼看懂": {"入選病原數": 1, "最終入選病原": ["A"]},
        "步驟2至5_逐候選病原判讀": [
            {
                "病原": "A",
                "最終處置": "選入",
                "最終排序": 1,
                "為什麼選入或排除": ["原理由 A"],
                "為什麼會被考慮": [],
                "重要限制": [],
                "關鍵判斷資訊": {},
            },
            {
                "病原": "B",
                "最終處置": "未選入",
                "最終排序": None,
                "為什麼選入或排除": ["原理由 B"],
                "為什麼會被考慮": [],
                "重要限制": [],
                "關鍵判斷資訊": {},
            },
        ],
    }
        decision = {
        "patient_id": "1",
        "stage": "workflow_final",
        "status": "computational_final_pending_rationale_and_clinical_review",
        "selected": [
            {"organism": "A", "rank": 1, "selection_origin": "upstream_picked", "rule_path": "picked_immutable"},
            {"organism": "B", "rank": 2, "selection_origin": "OBER_R5_rescue", "rule_path": "A1_site_B"},
        ],
        "candidates": [
            {"organism": "A", "selection_origin": "upstream_picked", "reasons": ["保留 picked"]},
            {"organism": "B", "selection_origin": "OBER_R5_rescue", "reasons": ["R5 補入"], "R5_signals": {"rule_path": "A1_site_B"}},
        ],
    }
        with tempfile.TemporaryDirectory() as temp_dir:
            decision_file = Path(temp_dir) / "P1_final_decision.json"
            decision_file.write_text(json.dumps(decision), encoding="utf-8")

            output = MODULE.apply_final_decision(chain, decision, decision_file)

        self.assertEqual(chain["一眼看懂"]["入選病原數"], 1)
        self.assertEqual(output["一眼看懂"]["最終入選病原"], ["A", "B"])
        self.assertEqual(output["步驟2至5_逐候選病原判讀"][1]["最終處置"], "選入")
        self.assertEqual(output["步驟2至5_逐候選病原判讀"][1]["最終排序"], 2)
        self.assertEqual(output["步驟6_最終鑑別診斷排序"][1]["規則路徑"], "A1_site_B")


    def test_apply_final_decision_rejects_candidate_mismatch(self):
        chain = {
        "病例識別": {"病人編號": "1"},
        "步驟2至5_逐候選病原判讀": [{"病原": "A"}],
    }
        decision = {
        "patient_id": "1",
        "stage": "workflow_final",
        "selected": [],
        "candidates": [{"organism": "B", "reasons": []}],
    }
        with tempfile.TemporaryDirectory() as temp_dir:
            decision_file = Path(temp_dir) / "P1_final_decision.json"
            decision_file.write_text(json.dumps(decision), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "candidate mismatch"):
                MODULE.apply_final_decision(chain, decision, decision_file)


if __name__ == "__main__":
    unittest.main()
