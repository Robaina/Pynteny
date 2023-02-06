#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests for the preprocessing module
"""

import tempfile
import unittest
from pathlib import Path

from pynteny.preprocessing import (
    FASTA,
    LabelledFASTA,
    is_legit_DNA_sequence,
    is_legit_peptide_sequence,
)

this_file_dir = Path(__file__).parent


class TestRecordSequence(unittest.TestCase):
    def test_is_legit_peptide_sequence(self):
        good_seq = "MSSLRQIAFYGKGGIGKSTTSQNTLAALVEMGQKILIVGCDPKADSTRLILNTKMQDTVL"
        bad_seq = "XXXXXXXXX12GKGGIGKSTTSQNTLAALVEMGQKILIVGCDPKADSTRLILNTKMQDTVL"
        self.assertTrue(
            is_legit_peptide_sequence(good_seq),
            "Failed to assert peptide sequence",
        )
        self.assertFalse(
            is_legit_peptide_sequence(bad_seq),
            "Failed to assert peptide sequence",
        )

    def test_is_legit_DNA_sequence(self):
        good_seq = "ACATCGTTTCCCCTGTTTCCACAAGACCTACTACGGCTGTTTTCGTAGTTCTTTTAAGAG"
        bad_seq = "XXXXX12XXXCCCTGTTTCCACAAGACCTACTACGGCTGTTTTCGTAGTTCTTTTAAGAG"
        self.assertTrue(
            is_legit_DNA_sequence(good_seq),
            "Failed to assert nucleotide sequence",
        )
        self.assertFalse(
            is_legit_DNA_sequence(bad_seq),
            "Failed to assert nucleotide sequence",
        )


class TestLabelledFASTA(unittest.TestCase):
    def test_from_genbank(self):
        with tempfile.NamedTemporaryFile() as tmp:
            labelled_FASTA = LabelledFASTA.from_genbank(
                this_file_dir / "test_data/MG1655.gb", output_file=tmp.name
            )
            self.assertIsInstance(
                labelled_FASTA,
                LabelledFASTA,
                "Failed to contruct labelled FASTA from GeneBank file",
            )


class TestFASTA(unittest.TestCase):
    def test_remove_duplicates(self):
        fasta = FASTA(this_file_dir / "test_data/MG1655.fasta")
        with tempfile.NamedTemporaryFile() as tmp:
            fasta.remove_duplicates(output_file=Path(tmp.name))
            self.assertEqual(
                fasta.file_path,
                Path(tmp.name),
                "Failed to remove duplicates in FASTA",
            )

    def test_remove_corrupted_sequences(self):
        fasta = FASTA(this_file_dir / "test_data/MG1655.fasta")
        with tempfile.NamedTemporaryFile() as tmp:
            fasta.remove_corrupted_sequences(output_file=Path(tmp.name))
            self.assertEqual(
                fasta.file_path,
                Path(tmp.name),
                "Failed to remove corrupted sequences in FASTA",
            )

    def test_filter_by_IDs(self):
        fasta = FASTA(this_file_dir / "test_data/MG1655.fasta")
        with tempfile.NamedTemporaryFile() as tmp:
            fasta.filter_by_IDs(
                record_ids=["b0001__U00096_0_189_255_pos"], output_file=Path(tmp.name)
            )
            self.assertEqual(
                fasta.file_path,
                Path(tmp.name),
                "Failed to filter records by ID in FASTA",
            )

    def test_filter_by_minimum_length(self):
        fasta = FASTA(this_file_dir / "test_data/MG1655.fasta")
        with tempfile.NamedTemporaryFile() as tmp:
            fasta.filter_by_minimum_length(min_length=50, output_file=Path(tmp.name))
            self.assertEqual(
                fasta.file_path,
                Path(tmp.name),
                "Failed to filter FASTA by length",
            )

    def test_split_by_contigs(self):
        fasta = FASTA(this_file_dir / "test_data/MG1655.fasta")
        with tempfile.TemporaryDirectory() as tmp:
            fasta.split_by_contigs(output_dir=tmp)
            self.assertGreater(
                len(list(Path(tmp).iterdir())), 1, "Failed to split contigs in FASTA"
            )


if __name__ == "__main__":
    unittest.main()
