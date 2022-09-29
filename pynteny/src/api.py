#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Classes to facilitate usage within Python scripts / Notebooks
"""

from pathlib import Path

from pynteny.src.utils import CommandArgs
from pynteny.src.filter import SyntenyHits
from pynteny.src.subcommands import (
    synteny_search, build_database,
    download_hmms, parse_gene_ids, get_citation
    )


class Search():
    def __init__(
        self,
        data: Path,
        synteny_struc: str,
        gene_ids: bool = False,
        unordered: bool = False,
        hmm_dir: Path = None,
        hmm_meta: Path = None,
        outdir: Path = None,
        prefix: str = "",
        hmmsearch_args: str = None,
        logfile: Path = None,
        processes: int = None):

        """
        Query sequence database for HMM hits arranged in provided synteny structure.
        """

        self._args = CommandArgs(
            data=self.data,
            synteny_struc=self.synteny_struc,
            hmm_dir=self.hmm_dir,
            hmm_meta=self.hmm_meta,
            outdir=self.outdir,
            prefix=self.prefix,
            hmmsearch_args=self.hmmsearch_args,
            gene_ids=self.gene_ids,
            logfile=self.logfile,
            processes=self.processes,
            unordered=self.unordered,
            )

    def run(self) -> SyntenyHits:
        """
        Run pynteny search
        """
        return synteny_search(self._args)


class Build():
    def __init__(self):
        pass


class Download():
    def __init__(self):
        pass