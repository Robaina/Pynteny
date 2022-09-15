#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests for the preprocessing module
"""

import unittest
from pynteny.src.preprocessing import RecordSequence


class TestRecordSequence(unittest.TestCase):
    def test_isLegitPeptideSequence(self):
        good_seq = "MSSLRQIAFYGKGGIGKSTTSQNTLAALVEMGQKILIVGCDPKADSTRLILNTKMQDTVL"
        bad_seq = "XXXXXXXXX12GKGGIGKSTTSQNTLAALVEMGQKILIVGCDPKADSTRLILNTKMQDTVL"
        self.assertTrue(
            RecordSequence.isLegitPeptideSequence(good_seq),
            "Failed to assert peptide sequence"
            )
        self.assertFalse(
            RecordSequence.isLegitPeptideSequence(bad_seq),
            "Failed to assert peptide sequence"
            )
    def test_isLegitDNASequence(self):
        good_seq = "ACATCGTTTCCCCTGTTTCCACAAGACCTACTACGGCTGTTTTCGTAGTTCTTTTAAGAG"
        bad_seq = "XXXXX12XXXCCCTGTTTCCACAAGACCTACTACGGCTGTTTTCGTAGTTCTTTTAAGAG"
        self.assertTrue(
            RecordSequence.isLegitDNAsequence(good_seq),
            "Failed to assert nucleotide sequence"
            )
        self.assertFalse(
            RecordSequence.isLegitDNAsequence(bad_seq),
            "Failed to assert nucleotide sequence"
            )


if __name__ == '__main__':
    unittest.main()