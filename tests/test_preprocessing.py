#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests for the preprocessing module
"""

import unittest
from pynteny.preprocessing import RecordSequence


class TestRecordSequence(unittest.TestCase):
    def test_is_legit_peptide_sequence(self):
        good_seq = "MSSLRQIAFYGKGGIGKSTTSQNTLAALVEMGQKILIVGCDPKADSTRLILNTKMQDTVL"
        bad_seq = "XXXXXXXXX12GKGGIGKSTTSQNTLAALVEMGQKILIVGCDPKADSTRLILNTKMQDTVL"
        self.assertTrue(
            RecordSequence.is_legit_peptide_sequence(good_seq),
            "Failed to assert peptide sequence",
        )
        self.assertFalse(
            RecordSequence.is_legit_peptide_sequence(bad_seq),
            "Failed to assert peptide sequence",
        )

    def test_is_legit_DNA_sequence(self):
        good_seq = "ACATCGTTTCCCCTGTTTCCACAAGACCTACTACGGCTGTTTTCGTAGTTCTTTTAAGAG"
        bad_seq = "XXXXX12XXXCCCTGTTTCCACAAGACCTACTACGGCTGTTTTCGTAGTTCTTTTAAGAG"
        self.assertTrue(
            RecordSequence.is_legit_DNA_sequence(good_seq),
            "Failed to assert nucleotide sequence",
        )
        self.assertFalse(
            RecordSequence.is_legit_DNA_sequence(bad_seq),
            "Failed to assert nucleotide sequence",
        )


if __name__ == "__main__":
    unittest.main()
