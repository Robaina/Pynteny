#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests for the filter module
"""
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from pynteny.hmm import HMMER, PGAP

this_file_dir = Path(__file__).parent


class TestPGAP(unittest.TestCase):
    def test_remove_missing_HMMs_from_metadata(self):
        with tempfile.NamedTemporaryFile() as tmp:
            pgap = PGAP(this_file_dir / "test_data/hmm_meta.tsv")
            pgap.remove_missing_HMMs_from_metadata(
                this_file_dir / "test_data/hmms", tmp
            )
            new_labels = pd.read_csv(tmp.name, sep="\t").label.values.tolist()
            original_labels = pgap._meta.label.values.tolist()
        self.assertListEqual(
            original_labels, new_labels, "Failed to remove missing HMMs from metadata"
        )

    def test_get_HMM_group_for_gene_symbol(self):
        pgap = PGAP(this_file_dir / "test_data/hmm_meta.tsv")
        hmm_group = pgap.get_HMM_group_for_gene_symbol("leuD")
        self.assertEqual(
            hmm_group, "TIGR00171.1", "Failed to retrieve HMM group from gene symbol"
        )

    def test_get_meta_info_for_HMM(self):
        pgap = PGAP(this_file_dir / "test_data/hmm_meta.tsv")
        meta = pgap.get_meta_info_for_HMM("TIGR00170.1")
        self.assertEqual(
            meta["label"], "leuC", "Failed to retrieve meta info for HMM name"
        )

    def test_get_HMM_gene_ID(self):
        pgap = PGAP(this_file_dir / "test_data/hmm_meta.tsv")
        gene_id = pgap.get_HMM_gene_ID("TIGR00170.1")[0]
        self.assertEqual(gene_id, "leuC", "Failed to retrieve gene id from HMM name")


class TestHMMER(unittest.TestCase):
    def test_get_HMMER_tables(self):
        input_hmms = [file for file in (this_file_dir / "test_data/hmms").iterdir()]
        with tempfile.TemporaryDirectory() as tmp:
            hmmer = HMMER(
                input_hmms=input_hmms,
                hmm_output_dir=tmp,
                input_data=this_file_dir / "test_data/MG1655.fasta",
                additional_args=[None for _ in range(len(input_hmms))],
            )
            hits = hmmer.get_HMMER_tables(reuse_hmmer_results=True)
            self.assertGreaterEqual(
                list(hits.values())[0].shape[0], 1, "Failed to retrieve HMMER hits"
            )


if __name__ == "__main__":
    unittest.main()
