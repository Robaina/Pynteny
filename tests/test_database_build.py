#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests for class Database
"""

import unittest
import tempfile
from pathlib import Path
from pynteny.subcommands import Database
from pynteny.api import Build
from pynteny.preprocessing import LabelledFASTA


this_file_dir = Path(Path(__file__).parent)


class TestDatabase(unittest.TestCase):
    def test_build_fasta(self):
        data = Path(this_file_dir / "test_data/test_assembly")
        database = Database(data)
        with tempfile.NamedTemporaryFile() as outfile:
            labelled_database = database.build(output_file=outfile.name)
            self.assertIsInstance(
                labelled_database,
                LabelledFASTA,
                "Failed to build database from assembly data",
            )

    def test_build(self):
        with tempfile.NamedTemporaryFile() as outfile:
            labelled_database = Build(
                data=Path(this_file_dir / "test_data/test_assembly"),
                outfile=outfile.name,
                logfile=None,
            ).run()

    def test_build_gbk(self):
        data = Path(this_file_dir / "test_data/MG1655.gb")
        database = Database(data)
        with tempfile.NamedTemporaryFile() as outfile:
            labelled_database = database.build(output_file=outfile.name)
            self.assertIsInstance(
                labelled_database,
                LabelledFASTA,
                "Failed to build database from GenBank data",
            )
