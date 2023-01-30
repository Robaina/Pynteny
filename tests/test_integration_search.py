#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Integration tests for the search subcommand
"""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from pynteny.api import Search

tests_dir = Path(__file__).parent


class TestSyntenySearch(unittest.TestCase):
    def test_search(self):
        with TemporaryDirectory() as tempdir:
            search = Search(
                data=tests_dir / "test_data/MG1655.fasta",
                synteny_struc="<TIGR00171.1 0 <TIGR00170.1 1 <TIGR00973.1",
                hmm_dir=tests_dir / "test_data/hmms",
                hmm_meta=None,
                outdir=Path(tempdir),
                prefix="",
                hmmsearch_args=None,
                gene_ids=False,
                reuse=True,
                logfile=None,
                processes=None,
                unordered=False,
            )
            hmm_meta = tests_dir / "test_data/hmm_meta.tsv"
            synhits = search.run()
            synhits = synhits.add_HMM_meta_info_to_hits(hmm_meta)
            synhits_df = synhits.hits
            config = Path(tests_dir.parent) / "config.json"
            if config.exists():
                config.unlink()
        hit_labels = [
            "b0071__U00096_71_78847_79453_neg",
            "b0072__U00096_72_79463_80864_neg",
            "b0074__U00096_74_81957_83529_neg",
        ]
        self.assertListEqual(
            synhits_df.full_label.values.tolist(),
            hit_labels,
            "Failed finding synteny hits",
        )
