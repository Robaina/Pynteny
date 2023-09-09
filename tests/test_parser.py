#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests for the filter module
"""

import unittest
from pathlib import Path

import pynteny.parsers.labelparser as labelparser
import pynteny.parsers.syntenyparser as syntenyparser

this_file_dir = Path(__file__).parent


class TestLabelParser(unittest.TestCase):
    def test_parse(self):
        label = (
            "Afifella_marina_BN_126_MMP03080610__FMVW01000002.1_552_0584494_0585393_pos"
        )
        parsed_dict = labelparser.parse(label)
        self.assertDictEqual(
            parsed_dict,
            {
                "full_label": "Afifella_marina_BN_126_MMP03080610__FMVW01000002.1_552_0584494_0585393_pos",
                "gene_id": "Afifella_marina_BN_126_MMP03080610",
                "contig": "FMVW01000002.1",
                "gene_pos": 552,
                "locus_pos": (584494, 585393),
                "strand": "pos",
            },
            "Failed to parse label",
        )


class TestSyntenyParser(unittest.TestCase):
    syn_struct = "<(TIGR00171.1|TIGR02084.1) 0 <(TIGR00170.1|TIGR02083.1) 1 <(NF002084.0|TIGR00973.1|TIGR00970.1)"

    def test_is_valid_structure(self):
        self.assertTrue(
            syntenyparser.is_valid_structure(self.syn_struct),
            "Failed to assess correct synteny structure format",
        )

    def test_split_strand_from_locus(self):
        locus_str = "<TIGR00171.1"
        self.assertEqual(
            syntenyparser.split_strand_from_locus(locus_str), ("neg", "TIGR00171.1")
        )

    def test_get_HMM_groups_in_structure(self):
        self.assertEqual(
            syntenyparser.get_HMM_groups_in_structure(self.syn_struct),
            [
                "TIGR00171.1|TIGR02084.1",
                "TIGR00170.1|TIGR02083.1",
                "NF002084.0|TIGR00973.1|TIGR00970.1",
            ],
            "Failed to retrieve HMM groups correctly",
        )

    def test_get_maximum_distances_in_structure(self):
        self.assertEqual(
            syntenyparser.get_maximum_distances_in_structure(self.syn_struct),
            [0, 1],
            "Failed to retrieve max ORF distances correctly",
        )

    def test_parse_genes_in_synteny_structure(self):
        meta = Path(this_file_dir / "test_data/hmm_meta.tsv")
        parsed_struct, hmm_groups = syntenyparser.parse_genes_in_synteny_structure(
            synteny_structure="<leuD 0 <leuC 1 <leuA", hmm_meta=meta
        )
        self.assertEqual(
            parsed_struct,
            "<(TIGR00171.1|TIGR02084.1) 0 <(TIGR00170.1|TIGR02083.1) 1 <(TIGR00973.1|NF002084.0|TIGR00970.1)",
            "Failed to parse gene symbols in synteny structure",
        )

    def test_parse_synteny_structure(self):
        parsed_struct = syntenyparser.parse_synteny_structure("<leuD 0 <leuC 1 <leuA")
        self.assertEqual(
            parsed_struct["distances"], [0, 1], "Failed to parse synteny structure"
        )


if __name__ == "__main__":
    unittest.main()
