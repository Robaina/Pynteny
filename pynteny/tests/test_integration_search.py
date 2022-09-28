#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Integration tests for the search subcommand
"""

import os
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from pynteny.src.utils import CommandArgs
from pynteny.src.subcommands import synteny_search


this_file_dir = Path(Path(__file__).parent)

class TestSyntenySearch(unittest.TestCase):
    def test_search(self):
        with TemporaryDirectory() as tempdir:
            search_state = CommandArgs(
                data=this_file_dir / "test_data/MG1655.fasta",
                synteny_struc="<TIGR00171.1 0 <TIGR00170.1 1 <TIGR00973.1",
                hmm_dir=this_file_dir / "test_data/hmms",
                hmm_meta=None,
                outdir=Path(tempdir),
                prefix="",
                hmmsearch_args=None,
                gene_ids=False,
                logfile=None,
                processes=None,
                unordered=False,
                )
            synhits = synteny_search(search_state).getSyntenyHits()
            os.remove(os.path.join(this_file_dir.parent, "config.json"))
        hit_labels = [
            'b0071__U00096_71_78847_79453_neg',
            'b0072__U00096_72_79463_80864_neg',
            'b0074__U00096_74_81957_83529_neg'
            ]
        self.assertListEqual(
            synhits.full_label.values.tolist(),
            hit_labels,
            "Failed finding synteny hits"
        )
