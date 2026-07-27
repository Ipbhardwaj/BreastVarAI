"""
Unit Test Suite for Variant Parser Module
"""

import unittest
from src.breastvar_ai.variant_parser import parse_hgvs_protein, AA_3TO1


class TestVariantParser(unittest.TestCase):

    def test_parse_hgvs_protein_valid(self):
        ref, pos, alt = parse_hgvs_protein("NM_007294.4(BRCA1):c.5096G>A (p.Arg1699Gln)")
        self.assertEqual(ref, "R")
        self.assertEqual(pos, 1699)
        self.assertEqual(alt, "Q")

    def test_parse_hgvs_protein_invalid(self):
        ref, pos, alt = parse_hgvs_protein("Invalid_Variant_String")
        self.assertIsNone(ref)
        self.assertIsNone(pos)
        self.assertIsNone(alt)

    def test_amino_acid_dictionary(self):
        self.assertEqual(AA_3TO1["Ala"], "A")
        self.assertEqual(AA_3TO1["Val"], "V")


if __name__ == "__main__":
    unittest.main()