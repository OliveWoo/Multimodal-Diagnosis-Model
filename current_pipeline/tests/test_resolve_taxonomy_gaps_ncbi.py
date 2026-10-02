import unittest
import xml.etree.ElementTree as ET

from tools.resolve_taxonomy_gaps_ncbi import parse_taxa


class ParseNcbiTaxonomyTest(unittest.TestCase):
    def test_parse_taxa(self):
        root = ET.fromstring(
            """
            <TaxaSet><Taxon><TaxId>1</TaxId><ScientificName>Example species</ScientificName>
            <Rank>species</Rank><Division>Bacteria</Division><Lineage>cellular organisms; Bacteria</Lineage>
            <LineageEx><Taxon><TaxId>2</TaxId><ScientificName>Example</ScientificName><Rank>genus</Rank></Taxon></LineageEx>
            </Taxon></TaxaSet>
            """
        )
        record = parse_taxa(root)["1"]
        self.assertEqual(record["scientific_name"], "Example species")
        self.assertEqual(record["lineage_ex"][0]["rank"], "genus")


if __name__ == "__main__":
    unittest.main()
