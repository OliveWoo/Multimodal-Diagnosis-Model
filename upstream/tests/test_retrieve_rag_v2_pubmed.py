from __future__ import annotations

import unittest

from tools import retrieve_rag_v2_pubmed as retrieval


SAMPLE_XML = b"""<?xml version="1.0" encoding="UTF-8"?>
<PubmedArticleSet>
  <PubmedArticle>
    <MedlineCitation>
      <PMID>12345678</PMID>
      <Article>
        <Journal><JournalIssue><PubDate><Year>2024</Year></PubDate></JournalIssue><Title>Test Journal</Title></Journal>
        <ArticleTitle>Candida tropicalis in adult pulmonary infection</ArticleTitle>
        <Abstract><AbstractText Label="BACKGROUND">Respiratory colonization is common.</AbstractText></Abstract>
        <PublicationTypeList><PublicationType>Systematic Review</PublicationType></PublicationTypeList>
      </Article>
      <MeshHeadingList><MeshHeading><DescriptorName>Pneumonia</DescriptorName></MeshHeading></MeshHeadingList>
    </MedlineCitation>
    <PubmedData><ArticleIdList><ArticleId IdType="pubmed">12345678</ArticleId><ArticleId IdType="doi">10.1/test</ArticleId></ArticleIdList></PubmedData>
  </PubmedArticle>
</PubmedArticleSet>
"""


class PubMedRetrievalV2Tests(unittest.TestCase):
    def test_parse_pubmed_xml(self) -> None:
        records = retrieval.parse_pubmed_xml(SAMPLE_XML)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["pmid"], "12345678")
        self.assertEqual(records[0]["year"], 2024)
        self.assertEqual(records[0]["publication_types"], ["Systematic Review"])
        self.assertIn("Respiratory colonization", records[0]["abstract"])

    def test_query_plan_has_pathogenicity_and_counterevidence(self) -> None:
        request = {
            "organism": {
                "display_name": "Candida tropicalis",
                "pathogen_category": "Candida/yeast",
            }
        }
        tracks = retrieval.query_tracks(request)
        names = {track["track"] for track in tracks}
        self.assertIn("pulmonary_pathogenicity", names)
        self.assertIn("colonization_or_background", names)
        self.assertIn("guideline_or_review", names)

    def test_article_score_rewards_review_and_counterevidence(self) -> None:
        article = retrieval.parse_pubmed_xml(SAMPLE_XML)[0]
        score, reasons = retrieval.article_score(
            article,
            "Candida tropicalis",
            "Candida/yeast",
            {"pulmonary_pathogenicity", "colonization_or_background"},
        )
        self.assertGreater(score, 20)
        self.assertIn("species_name_in_title", reasons)
        self.assertIn("colonization_or_background_context_present", reasons)

    def test_non_pulmonary_background_article_fails_track_gate(self) -> None:
        article = {
            "title": "Ocular HSV-1 latency and reactivation",
            "abstract": "HSV-1 reactivation causes recurrent ocular disease.",
            "publication_types": ["Review"],
        }
        self.assertFalse(
            retrieval.article_matches_track(
                article,
                "HSV-1",
                "Herpesvirus/reactivation",
                "colonization_or_background",
                "species_specific",
            )
        )

    def test_animal_only_article_fails_track_gate(self) -> None:
        article = {
            "title": "Bacteroides fragilis lung abscess in rabbits",
            "abstract": "Bacteroides fragilis caused an experimental lung abscess.",
            "publication_types": ["Journal Article"],
            "mesh_terms": ["Animals", "Rabbits", "Lung Abscess"],
        }
        self.assertFalse(
            retrieval.article_matches_track(
                article,
                "Bacteroides fragilis",
                "Strict anaerobe/aspiration flora",
                "pulmonary_pathogenicity",
                "species_specific",
            )
        )

    def test_selection_does_not_pad_with_low_score_articles(self) -> None:
        articles = [
            {
                "pmid": "1",
                "retrieval_tracks": ["pulmonary_pathogenicity"],
                "retrieval_relevance_score": 24.0,
            },
            {
                "pmid": "2",
                "retrieval_tracks": ["pulmonary_pathogenicity"],
                "retrieval_relevance_score": 12.0,
            },
        ]
        selected = retrieval.select_diverse_articles(articles, max_articles=15)
        self.assertEqual([article["pmid"] for article in selected], ["1"])


if __name__ == "__main__":
    unittest.main()
