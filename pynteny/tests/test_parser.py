#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests for the filter module
"""

import unittest
from pynteny.src.filter import LabelParser, SyntenyParser


class TestLabelParser(unittest.TestCase):
    def test_parse(self):
        label = "Afifella_marina_BN_126_MMP03080610__FMVW01000002.1_552_0584494_0585393_pos"
        parsed_dict = LabelParser.parse(label)
        self.assertDictEqual(
            parsed_dict,
            {
                'full_label': 'Afifella_marina_BN_126_MMP03080610__FMVW01000002.1_552_0584494_0585393_pos',
                'gene_id': 'Afifella_marina_BN_126_MMP03080610',
                'contig': 'FMVW01000002.1',
                'gene_pos': 552,
                'locus_pos': (584494, 585393),
                'strand': 'pos'
                },
            "Failed to parse label"
            )

class TestSyntenyParser(unittest.TestCase):
    syn_struct = "<(TIGR00171.1|TIGR02084.1) 0 <(TIGR00170.1|TIGR02083.1) 1 <(NF002084.0|TIGR00973.1|TIGR00970.1)"
    def test_isValidStructure(self):
        self.assertTrue(
            SyntenyParser.isValidStructure(self.syn_struct),
            "Failed to assess correct synteny structure format"
        )
    def test_splitStrandFromLocus(self):
        locus_str = "<TIGR00171.1"
        self.assertEqual(
            SyntenyParser.splitStrandFromLocus(locus_str),
            ("neg", "TIGR00171.1")
        )
    def test_getHMMgroupsInStructure(self):
        self.assertEqual(
            SyntenyParser.getHMMgroupsInStructure(self.syn_struct),
            ["TIGR00171.1|TIGR02084.1", "TIGR00170.1|TIGR02083.1", "NF002084.0|TIGR00973.1|TIGR00970.1"],
            "Failed to retrieve HMM groups correctly"
        )
    def test_getMaximumDistancesInStructure(self):
        self.assertEqual(
            SyntenyParser.getMaximumDistancesInStructure(self.syn_struct),
            [0, 1],
            "Failed to retrieve max ORF distances correctly"
        )
    

if __name__ == '__main__':
    unittest.main()




