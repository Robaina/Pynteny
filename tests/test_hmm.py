#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests for the filter module
"""
import unittest
import tempfile
from pathlib import Path
import pandas as pd
from pynteny.hmm import PGAP, HMMER

this_file_dir = Path(Path(__file__).parent)


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


if __name__ == "__main__":
    unittest.main()
