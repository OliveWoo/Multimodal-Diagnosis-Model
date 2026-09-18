import unittest

from rag_re.pubmed import build_query


class PubMedQueryTests(unittest.TestCase):
    def test_query_uses_locked_title_abstract_and_publication_date_fields(self):
        query = build_query(
            "Pneumocystis jirovecii",
            {
                "query_terms": "pneumonia OR lower respiratory tract infection",
                "target_site": "lower respiratory tract infection",
            },
            "2026-06-30",
        )
        self.assertIn('"Pneumocystis jirovecii"[Title/Abstract]', query)
        self.assertIn('"pneumonia"[Title/Abstract]', query)
        self.assertIn("1900/01/01:2026/06/30[dp]", query)


if __name__ == "__main__":
    unittest.main()
