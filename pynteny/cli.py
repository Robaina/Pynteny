#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations
import sys
import random
import tempfile
import argparse
from importlib import metadata
from pathlib import Path

import pynteny.src.subcommands as sub
from pynteny.src.utils import ConfigParser


meta = metadata.metadata("pynteny")
__version__ = meta["Version"]
__author__ = meta["Author"]


class Pynteny():
    """
    Based on:
    https://selvakumar-arumugapandian.medium.com/command-line-subcommands-with-pythons-argparse-4dbac80f7110

    """
    def __init__(self, subcommand: str, subcommand_args: list[str]):
        self._subcommand = subcommand
        self._subcommand_args = subcommand_args
        parser = argparse.ArgumentParser(
            description=(self._printLogo()),
            usage=("pynteny <subcommand> [-h] [args] \n"),
            epilog=(self._generateCoolQuotes()),
            formatter_class=argparse.RawTextHelpFormatter
            )
        parser._positionals.title = "subcommands"
        parser.add_argument(
            help=(
                "search \n"
                "build \n"
                "parse \n"
                "download \n"
                "app \n"
                "tests \n"
                "cite \n"
                ),
            dest="subcommand",
            metavar="",
            )
        parser.add_argument("-v","--version", help="show version and exit", action="version", version=__version__)
        if len(sys.argv) < 2:
            parser.print_help()
            sys.exit(1)
        args = parser.parse_args(self._subcommand)
        input_subcommand = getattr(args, "subcommand")
        ConfigParser.initialize_config_file()
        self._call_subcommand(subcommand_name=input_subcommand)

    def _printLogo(self):
        print(
            (
            """
    ____              __                  
   / __ \__  ______  / /____  ____  __  __
  / /_/ / / / / __ \/ __/ _ \/ __ \/ / / /
 / ____/ /_/ / / / / /_/  __/ / / / /_/ / 
/_/    \__, /_/ /_/\__/\___/_/ /_/\__, /  
      /____/                     /____/   

"""
            f"Synteny-based Hmmer searches made easy, v{__version__}\n"
            "Semidán Robaina Estévez (srobaina@ull.edu.es), 2022\n"
            " \n"
            )
            )
    
    def _generateCoolQuotes(self):
        quotes = [
            "May the force be with you (Yoda)",
            "This looks like a job for a computer (AI)",
            "This is such a wonderful day to do bioinformatics (SR)",
            "One does not simply walk into Mordor (J.R.R. Tolkien)",
            "Damn, looks like a rainy day, let's do bioiformatics! (SR)",
        ]
        return(
            f"{random.choice(quotes)}\n"
            " "
        )
        
    def _call_subcommand(self, subcommand_name: str) -> None: 
        subcommand = getattr(self, subcommand_name)
        subcommand()
    
    def search(self):
        parser = SubcommandParser.search()
        args = parser.parse_args(self._subcommand_args)
        sub.synteny_search(args)

    def build(self):
        parser = SubcommandParser.build()
        args = parser.parse_args(self._subcommand_args)
        sub.build_database(args)

    def parse(self):
        parser = SubcommandParser.parse()
        args = parser.parse_args(self._subcommand_args)
        sub.parse_gene_ids(args)

    def download(self):
        parser = SubcommandParser.download()
        args = parser.parse_args(self._subcommand_args)
        sub.download_hmms(args)
    
    def app(self):
        """
        Run pynteny app through Streamlit
        """
        parser = SubcommandParser.app()
        args = parser.parse_args(self._subcommand_args)
        sub.run_app()

    def tests(self):
        """
        Run pynteny unit and integration tests
        """
        parser = SubcommandParser.tests()
        args = parser.parse_args(self._subcommand_args)
        sub.run_tests()

    def cite(self):
        """
        Print pynteny's citation string
        """
        parser = SubcommandParser.cite()
        args = parser.parse_args(self._subcommand_args)
        args.version = __version__
        sub.get_citation(args)


class SubcommandParser():
    """
    Argparse parsers for Pynteny's subcommands
    """
    @staticmethod
    def getHelpStr(subcommand: str) -> str:
        parser = getattr(SubcommandParser, subcommand)()
        with tempfile.NamedTemporaryFile(mode="w+") as file:
            parser.print_help(file)
            file.flush()
            with open(file.name) as help_file:
                help_str = help_file.read()
        return help_str

    @staticmethod
    def search() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(
            description=(
                "Query sequence database for HMM hits arranged in provided synteny structure."
                ),
            usage=("pynteny search [-h] [args] \n"),
            epilog="  \n",
            formatter_class=argparse.RawTextHelpFormatter
            )

        optional = parser._action_groups.pop()
        required = parser.add_argument_group("required arguments")
        parser._action_groups.append(optional)

        required.add_argument("-s", "--synteny_struc", 
                            metavar="", dest="synteny_struc", 
                            type=str, required=True,
                            help=(
                                f"string displaying hmm structure to search for, such as: \n"
                                f" \n"
                                f"'>hmm_a n_ab <hmm_b n_bc hmm_c'\n"
                                f" \n"
                                f"where '>' indicates a hmm target located on the positive strand, \n"
                                f"'<' a target located on the negative strand, and n_ab cooresponds \n"
                                f"to the maximum number of genes separating matched genes a and b. \n" 
                                f"Multiple hmms may be employed. \n"
                                f"No order symbol in a hmm indicates that results should be independent \n"
                                f"of strand location. "
                                )
        )
        required.add_argument("-i", "--data", dest="data", metavar="", type=Path, required=True,
                             help="path to peptide database"
        )
        optional.add_argument("-d", "--hmm_dir", dest="hmm_dir", type=Path, metavar="",
                            required=False, default=None,
                            help=(
                                "path to directory containing hmm (i.e, tigrfam or pfam) models. \n"
                                "The directory can contain more hmm models than used in the synteny structure. \n"
                                "It may also be the path to a compressed (tar, tar.gz, tgz) directory. \n"
                                "If not provided, hmm models (PGAP database) will be downloaded from the NCBI.\n"
                                "(if not already downloaded)"
                                )
        )
        optional.add_argument("-o", "--outdir", dest="outdir", type=Path, metavar="",
                             help="path to output directory", default=None
        )
        optional.add_argument("-x", "--prefix", dest="prefix", type=str, metavar="",
                            default="",
                            help="prefix to be added to output files"
        )
        optional.add_argument("-p", "--processes", dest="processes", type=int, metavar="",
                            default=None,
                            help="maximum number of processes available to HMMER. Defaults to all but one."
        )
        optional.add_argument("-a", "--hmmsearch_args", dest="hmmsearch_args", type=str,
                            metavar="", default=None, required=False,
                            help=(
                                "list of comma-separated additional arguments to hmmsearch for each input hmm. \n"
                                "A single argument may be provided, in which case the same additional argument \n"
                                "is employed in all hmms.")
        )
        optional.add_argument("-g", "--gene_ids", dest="gene_ids",
                             default=False, action="store_true",
                             help=(
                                "use gene symbols in synteny structure instead of HMM names. \n"
                                "If set, a path to the hmm database metadata file must be provided \n"
                                "in argument '--hmm_meta'"
                                )
        )
        optional.add_argument("-u", "--unordered", dest="unordered",
                             default=False, action="store_true",
                             help=(
                                "whether the HMMs should be arranged in the exact same order displayed \n"
                                "in the synteny_structure or in  any order If ordered, the filters would \n"
                                "filter collinear rather than syntenic structures."
                                )
        )
        optional.add_argument("-m", "--hmm_meta", dest="hmm_meta", type=Path, default=None,
                             metavar="",
                             required=False, help="path to hmm database metadata file"
        )
        optional.add_argument("-l", "--log", dest="logfile", type=Path, default=None,
                             metavar="",
                             required=False, help="path to log file. Log not written by default."
        )
        return parser
    
    @staticmethod
    def build() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(
            description=(
                "Translate nucleotide assembly file and assign contig and gene location info \n"
                "to each identified ORF (using prodigal). Label predicted ORFs according to \n"
                "positional info and export a fasta file containing predicted and translated ORFs. \n"
                "Alternatively, extract peptide sequences from GenBank file containing ORF annotations \n"
                "and write labelled peptide sequences to a fasta file."
                ),
            usage=("pynteny build [-h] [args] \n"),
            epilog="  \n",
            formatter_class=argparse.RawTextHelpFormatter
            )

        optional = parser._action_groups.pop()
        required = parser.add_argument_group("required arguments")
        parser._action_groups.append(optional)

        required.add_argument("-i", "--data", dest="data", type=Path,
                            required=True, metavar="",
                            help=(
                                "path to assembly input nucleotide data or annotated GenBank file. \n"
                                "It can be a single file or a directory of files (either of FASTA or GeneBank format)."
                                )
        )
        optional.add_argument("-o", "--outfile", dest="outfile", type=Path, metavar="",
                             default=None, help=(
                                "path to output (labelled peptide database) file. Defaults to \n"
                                "file in directory of input data."
                                )
        )
        optional.add_argument("-n", "--processes", dest="processes", type=int, metavar="",
                            required=False, default=None, help=(
                                "set maximum number of processes. "
                                "Defaults to all but one."
                                )
        )
        optional.add_argument("-l", "--log", dest="logfile", type=Path, default=None,
                             metavar="",
                             required=False, help="path to log file. Log not written by default."
        )
        return parser
    
    @staticmethod
    def parse() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(
            description=(
                "Translate synteny structure with gene symbols into one with\n"
                "HMM groups, according to provided HMM database."
                ),        
            usage=("pynteny parse [-h] [args] \n"),
            epilog="  \n",
            formatter_class=argparse.RawTextHelpFormatter
            )

        optional = parser._action_groups.pop()
        required = parser.add_argument_group("required arguments")
        parser._action_groups.append(optional)

        required.add_argument("-s", "--synteny_struc", dest="synteny_struc", type=str,
                            required=True, metavar="",
                            help=(
                                "synteny structure containing gene symbols instead of HMMs"
                                )
        )
        optional.add_argument("-m", "--hmm_meta", dest="hmm_meta", type=Path, metavar="",
                             required=False, help=(
                                "path to hmm database metadata file. If already donwloaded with \n"
                                "pynteny downloaded, hmm meta file is retrieved from default location."
                                )
                                )
        optional.add_argument("-l", "--log", dest="logfile", type=Path, default=None,
                             metavar="",
                             required=False, help="path to log file. Log not written by default."
        )
        return parser
    
    @staticmethod
    def download() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(
            description=(
                "Download HMM database from NCBI."
                ),
            epilog="  \n",
            usage=("pynteny download [-h] [args] \n"),
            formatter_class=argparse.RawTextHelpFormatter
            )

        optional = parser._action_groups.pop()
        required = parser.add_argument_group("required arguments")
        parser._action_groups.append(optional)

        optional.add_argument("-o", "--outdir", dest="outdir", type=Path, metavar="",
                             required=False, default=None,
                             help=(
                                "path to directory where to download HMM database.\n"
                                "Defaults to pynteny installation directory."
                                )
                             )
        optional.add_argument("-u", "--unpack", dest="unpack",
                            default=False, action="store_true",
                            help="unpack originally compressed database files"
                            )
        optional.add_argument("-l", "--log", dest="logfile", type=Path, default=None,
                             metavar="",
                             required=False, help="path to log file. Log not written by default."
        )
        return parser
    
    @staticmethod
    def app() -> argparse.ArgumentParser:
        """
        Run pynteny app through Streamlit
        """
        parser = argparse.ArgumentParser(
            description=(
                "Run Pynteny in web browser."
                ),
            epilog="  \n",
            usage=("pynteny app [-h] \n"),
            formatter_class=argparse.RawTextHelpFormatter
            )
        optional = parser._action_groups.pop()
        required = parser.add_argument_group("required arguments")
        parser._action_groups.append(optional)
        return parser
    
    @staticmethod
    def tests() -> argparse.ArgumentParser:
        """
        Run pynteny unit and integration tests
        """
        parser = argparse.ArgumentParser(
            description=(
                "Run Pynteny unit and integration tests."
                ),
            epilog="  \n",
            usage=("pynteny tests [-h] \n"),
            formatter_class=argparse.RawTextHelpFormatter
            )
        optional = parser._action_groups.pop()
        required = parser.add_argument_group("required arguments")
        parser._action_groups.append(optional)
        return parser
    
    @staticmethod
    def cite() -> argparse.ArgumentParser:
        """
        Print pynteny's citation string
        """
        parser = argparse.ArgumentParser(
            description=(
                "Print pynteny's citation string."
                ),
            epilog="  \n",
            usage=("pynteny cite [-h] \n"),
            formatter_class=argparse.RawTextHelpFormatter
            )
        optional = parser._action_groups.pop()
        required = parser.add_argument_group("required arguments")
        parser._action_groups.append(optional)
        return parser


def main():
    subcommand, subcommand_args = sys.argv[1:2], sys.argv[2:]
    pynteny = Pynteny(subcommand, subcommand_args)
    
if __name__ == "__main__":
    main()
