from __future__ import annotations

import unittest
from datetime import datetime
from pathlib import Path
from tempfile import TemporaryDirectory

from tools import build_two_hospital_temporal_culture_overlay as overlay


class TemporalCultureOverlayTests(unittest.TestCase):
    def test_window_is_inclusive_and_monotonic(self) -> None:
        index = datetime(2025, 6, 4, 0, 1)
        early = overlay.alignment(index, datetime(2025, 6, 2, 0, 0))
        late_boundary = overlay.alignment(index, datetime(2025, 6, 6, 23, 59))
        self.assertTrue(early["event_aligned"])
        self.assertTrue(late_boundary["event_aligned"])
        self.assertEqual(late_boundary["signed_calendar_day_delta"], 2)
        self.assertEqual(
            late_boundary["event_alignment_status"],
            "within_plus_minus_two_calendar_dates",
        )
        late = overlay.alignment(index, datetime(2025, 6, 8))
        self.assertFalse(late["event_aligned"])
        self.assertEqual(late["signed_calendar_day_delta"], 4)
        self.assertEqual(late["event_alignment_status"], "post_index_outside_window")

    def test_crkp_is_identity_plus_resistance_phenotype(self) -> None:
        self.assertEqual(
            overlay.pathogen_normalization.canonical_key("CRKP"),
            "klebsiellapneumoniae",
        )
        self.assertEqual(
            overlay.pathogen_normalization.display_name("CRKP"),
            "Klebsiella pneumoniae",
        )
        self.assertEqual(
            overlay.pathogen_normalization.resistance_note("CRKP"),
            "CRKP phenotype",
        )

    def test_synthetic_workbook_preserves_three_late_resistance_links(self) -> None:
        with TemporaryDirectory() as temporary:
            path = Path(temporary) / "synthetic_cultures.xlsx"
            workbook = overlay.openpyxl.Workbook()
            sheet = workbook.active
            sheet.title = overlay.WORKSHEET_NAME
            sheet.append(["column"] * 75)
            values = [None] * 75
            values[0:2] = ["SYNTHETIC_CASE_1", "EXAMPLE_001"]
            for offset, column in enumerate(overlay.CULTURE_SLOT_START_COLUMNS[:3]):
                start = column - 1
                values[start] = "Sputum"
                values[start + 2] = "CRKP"
                values[start + 9] = datetime(2000, 1, 14 + offset)
                values[start + 10] = datetime(2000, 1, 17 + offset)
            sheet.append(values)
            workbook.save(path)
            workbook.close()
            workbook = overlay.openpyxl.load_workbook(path, read_only=True, data_only=True)
            try:
                row_number, row = overlay.find_source_row(
                    workbook[overlay.WORKSHEET_NAME],
                    hospital_case_id="SYNTHETIC_CASE_1",
                    specimen_id="EXAMPLE_001",
                )
                records = overlay.extract_culture_records(
                    row_number, row, index_time=datetime(2000, 1, 10), workbook_path=path
                )
            finally:
                workbook.close()
        crkp = [item for item in records if item.get("resistance_phenotype")]
        self.assertEqual(len(crkp), 3)
        self.assertEqual(
            {item["organism"] for item in crkp}, {"Klebsiella pneumoniae"}
        )
        self.assertTrue(all(not item["event_aligned"] for item in crkp))


if __name__ == "__main__":
    unittest.main()
