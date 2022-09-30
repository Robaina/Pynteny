#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Classes to facilitate usage within Python scripts / Notebooks
"""

from pathlib import Path
from importlib import metadata

from pynteny.src.utils import CommandArgs
from pynteny.src.filter import SyntenyHits
from pynteny.cli import SubcommandParser
from pynteny.src.subcommands import (
    synteny_search, build_database,
    download_hmms, parse_gene_ids, get_citation
    )


meta = metadata.metadata("pynteny")
__version__ = meta["Version"]
__author__ = meta["Author"]


class Command():
    def __init__(self):
        """
        Parent class for Pynteny command
        """

    def _repr_html_(self):
        """
        This method is executed automatically by Jupyter to print html
        """
        return """
        <table>
            <tr>
                <td><strong>Pynteny version</strong></td><td>{version}</td>
            </tr><tr>
                <td><strong>Author</strong></td><td>{author}</td>
            </tr>
        </table>
        """.format(version=__version__,
                   author=__author__)

    def update(self, argname: str, value: str) -> None:
        """
        Update argument value in pynteny command
        """
        setattr(self._args, argname, value)
    
    @staticmethod
    def cite():
        """
        Display Pynteny citation
        """
        args = CommandArgs(
            version=__version__,
            author=__author__
            )
        citation = get_citation(args, silent=True)
        print(citation)


class Search(Command):
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
            data=data,
            synteny_struc=synteny_struc,
            hmm_dir=hmm_dir,
            hmm_meta=hmm_meta,
            outdir=outdir,
            prefix=prefix,
            hmmsearch_args=hmmsearch_args,
            gene_ids=gene_ids,
            logfile=logfile,
            processes=processes,
            unordered=unordered,
            )

    def parseGeneIDs(self, synteny_struc: str) -> str:
        """
        Parse gene IDs in synteny structure and find corresponding HMM names
        in provided HMM database
        """
        args = CommandArgs(
            synteny_struc=synteny_struc,
            hmm_meta=self._args.hmm_meta,
            logfile=self._args.logfile
        )
        return parse_gene_ids(args)

    def run(self) -> SyntenyHits:
        """
        Run pynteny search
        """
        return synteny_search(self._args)


class Build(Command):
    def __init__(
        self,
        data: Path,
        outfile: Path,
        logfile: Path = None,
        processes: int = None):

        """
        Translate nucleotide assembly file and assign contig and gene location 
        info to each identified ORF (using prodigal).
        """

        self._args = CommandArgs(
            data=data,
            outfile=outfile,
            logfile=logfile,
            processes=processes
            )

    def run(self) -> None:
        """
        Run pynteny search
        """
        return build_database(self._args)


class Download(Command):
    def __init__(
        self,
        outdir: Path,
        logfile: Path,
        unpack: bool = False):

        """
        Download HMM database from NCBI.
        """

        self._args = CommandArgs(
            outdir=outdir,
            logfile=logfile,
            unpack=unpack
            )

    def run(self) -> None:
        """
        Run pynteny download
        """
        return download_hmms(self._args)