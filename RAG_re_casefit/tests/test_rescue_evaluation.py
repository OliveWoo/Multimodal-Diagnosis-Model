import unittest

from rag_re_casefit.rescue_evaluation import _maximum_match_count


class RescueEvaluationTests(unittest.TestCase):
    def test_one_to_one_exact_matching(self):
        self.assertEqual(
            _maximum_match_count(
                ("candida albicans", "candida tropicalis"),
                ("candida albicans",),
                genus_relaxed=False,
            ),
            1,
        )

    def test_genus_relaxed_is_still_one_to_one(self):
        self.assertEqual(
            _maximum_match_count(
                ("candida albicans", "candida tropicalis"),
                ("candida glabrata",),
                genus_relaxed=True,
            ),
            1,
        )

    def test_virus_names_are_not_genus_relaxed(self):
        self.assertEqual(
            _maximum_match_count(
                ("human cytomegalovirus",),
                ("human alphaherpesvirus 1",),
                genus_relaxed=True,
            ),
            0,
        )


if __name__ == "__main__":
    unittest.main()
