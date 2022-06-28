#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
import argparse
from pathlib import Path

import pynteny.src.subcommands as sub
from pynteny.src.utils import ConfigParser



class Pynteny():
    """
    Temporary fix to enable subcommands without too much
    code rewriting. Based on:
    https://selvakumar-arumugapandian.medium.com/command-line-subcommands-with-pythons-argparse-4dbac80f7110

    """
    def __init__(self):
        ConfigParser.initialize_config_file()
        self._subcommand = sys.argv[1:2]
        self._subcommand_args = sys.argv[2:]
        parser = argparse.ArgumentParser(
            description=(
                """

    ____              __                  
   / __ \__  ______  / /____  ____  __  __
  / /_/ / / / / __ \/ __/ _ \/ __ \/ / / /
 / ____/ /_/ / / / / /_/  __/ / / / /_/ / 
/_/    \__, /_/ /_/\__/\___/_/ /_/\__, /  
      /____/                     /____/   

"""
                "synteny-based Hmmer searches made easy\n"
                " \n"
                ),
            epilog=(
                "Semidán Robaina Estévez (srobaina@ull.edu.es), 2022\n"
                " "
                ),
            formatter_class=argparse.RawTextHelpFormatter
            )
        parser._positionals.title = "subcommands"
        parser.add_argument(
            help=(
                f"search \n"
                f"preprocess \n"
                f"parse \n"
                f"download \n"
                ),
            dest="subcommand"
            )
        parser.add_argument("-v","--version", help="show version and exit", action="version", version="0.0.1")
        if len(sys.argv) < 2:
            parser.print_help()
            sys.exit(1)
        args = parser.parse_args(self._subcommand)
        input_subcommand = getattr(args, "subcommand")
        self._call_subcommand(subcommand_name=input_subcommand)
    
    def _call_subcommand(self, subcommand_name: str) -> None: 
        subcommand = getattr(self, subcommand_name)
        subcommand()

    def search(self):
        parser = argparse.ArgumentParser(
            description=(
                'Build peptide database from HMM synteny structure. \n'
                'The script outputs a main file containing sequences matching the provided \n'
                'hmm structure and corresponding to the main target indicated in the argument: \n'
                'target_hmm. These sequences are filtered by sequence length and by maximum \n'
                'number of total sequences. \n'
                'The script also outputs additional file containing the matched records for the \n'
                'other, non-target, hmms. However, these files are not processed any further.'
                ),
            usage=("pynteny search [-h] ̣̣̣̣̣̣̣̣̣̣̣̣[args] ̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣\n"),
            epilog='Semidán Robaina Estévez (srobaina@ull.edu.es), 2022',
            formatter_class=argparse.RawTextHelpFormatter
            )

        optional = parser._action_groups.pop()
        required = parser.add_argument_group('required arguments')
        parser._action_groups.append(optional)

        required.add_argument('--hmm_dir', dest='hmm_dir', type=Path,
                            required=True,
                            help=(
                                'path to directory containing hmm (i.e, tigrfam or pfam) models. \n'
                                'The directory can contain more hmm models than used in the synteny structure. \n'
                                'It may also be the path to a compressed (tar, tar.gz, tgz) directory.'
                                )
        )
        required.add_argument('--synteny_struc', dest='synteny_struc', type=str, required=True,
                            help=(
                                f'string displaying hmm sctructure to search for, such as: \n'
                                f" \n"
                                f'">hmm_a n_ab <hmm_b n_bc hmm_c"\n'
                                f" \n"
                                f'where ">" indicates a hmm target located on the positive strand, \n'
                                f'"<" a target located on the negative strand, and n_ab cooresponds \n'
                                f'to the maximum number of genes separating matched gene a and b. \n' 
                                f'Multiple hmms may be employed (limited by computational capabilities). \n'
                                f'No order symbol in a hmm indicates that results should be independent \n'
                                f'of strand location. '
                                )
        )
        required.add_argument('--in', dest='data', type=Path, required=True,
                            help='path to peptide database'
        )
        optional.add_argument('--outdir', dest='outdir', type=Path,
                            help='path to output directory'
        )
        optional.add_argument('--prefix', dest='prefix', type=str,
                            default='',
                            help='prefix to be added to output files'
        )
        optional.add_argument('--min_seq_length', dest='minseqlength',
                            default=None, type=int,
                            help=(
                                'minimum sequence length in reference database. '
                                'Defaults to zero'
                                )
        )
        optional.add_argument('--hmmsearch_args', dest='hmmsearch_args', type=str,
                            default=None, required=False,
                            help=(
                                'list of comma-separated additional arguments to hmmsearch for each input hmm. \n'
                                'A single argument may be provided, in which case the same additional argument \n'
                                'is employed in all hmms.')
        )
        optional.add_argument('--max_seq_length', dest='maxseqlength',
                            default=None, type=int,
                            required=False,
                            help=(
                                'maximum sequence length in reference database. '
                                'Defaults to inf'
                                )
        )
        optional.add_argument("--gene_ids", dest="gene_ids",
                             default=False, action="store_true",
                             help=(
                                "use gene symbols in synteny structure instead of HMM nanmes. \n"
                                "If set, a path to the hmm database metadata file must be provided \n"
                                "in argument '--hmm_meta'"
                                )
        )
        required.add_argument('--hmm_meta', dest='hmm_meta', type=Path,
                            required=False, help='path to hmm database metadata file'
        )
        args = parser.parse_args(self._subcommand_args)
        sub.synteny_search(args)
    
    def preprocess(self):
        parser = argparse.ArgumentParser(
            description=(
                "Translate nucleotide assembly file and assign contig and gene location info \n"
                "to each identified ORF (using prodigal). Then label predicted ORFs according to \n"
                "positional info and export a fasta file containing predicted and translated ORFs."
                ),
            usage=("pynteny preprocess [-h] ̣̣̣̣̣̣̣̣̣̣̣̣[args] ̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣\n"),
            epilog="Semidán Robaina Estévez (srobaina@ull.edu.es), 2022",
            formatter_class=argparse.RawTextHelpFormatter
            )

        optional = parser._action_groups.pop()
        required = parser.add_argument_group("required arguments")
        parser._action_groups.append(optional)

        required.add_argument("--assembly_fasta", dest="assembly_fasta", type=Path,
                            required=True,
                            help=(
                                "path to assembly input nucleotide data. It can be a single \n"
                                "FASTA file or a directory containing several FASTA files."
                                )
        )
        optional.add_argument("--outdir", dest="outdir", type=Path,
                            help="path to output directory"
        )
        optional.add_argument("--prefix", dest="prefix", type=str,
                            default="prodigal",
                            help="prefix to be added to output files"
        )
        optional.add_argument("--processes", "-p", dest="processes", type=int,
                            required=False, default=None, help=(
                                "set maximum number of processes. "
                                "Defaults to all but one."
                                )
                                )
        optional.add_argument("--split_contigs", dest="split",
                            default=False, action="store_true",
                            help="split assembly input file into files containing one contig each")
        optional.add_argument("--metagenome", dest="metagenome",
                            default=False, action="store_true",
                            help="treat input assembly as of type metagenome (for prodigal)")
        optional.add_argument("--prodigal_args", dest="prodigal_args", type=str,
                            default=None, required=False,
                            help=("additional arguments to prodigal provided as defined by it")
                            )
        args = parser.parse_args(self._subcommand_args)
        sub.translate_assembly(args)

    def parse(self):
        parser = argparse.ArgumentParser(
            description=(
                "Translate synteny structure with gene symbols into one with\n"
                "HMM groups, according to provided HMM database."
                ),
            usage=("pynteny parse [-h] ̣̣̣̣̣̣̣̣̣̣̣̣[args] ̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣\n"),
            epilog="Semidán Robaina Estévez (srobaina@ull.edu.es), 2022",
            formatter_class=argparse.RawTextHelpFormatter
            )

        optional = parser._action_groups.pop()
        required = parser.add_argument_group("required arguments")
        parser._action_groups.append(optional)

        required.add_argument("--synteny_struc", dest="synteny_struc", type=str,
                            required=True,
                            help=(
                                "synteny structure containing gene symbols instead of HMMs"
                                )
        )
        required.add_argument("--hmm_meta", dest="hmm_meta", type=Path,
                             required=True, help="path to hmm database metadata file")
        args = parser.parse_args(self._subcommand_args)
        sub.parse_gene_ids(args)

    def download(self):
        parser = argparse.ArgumentParser(
            description=(
                "Download HMM database from NCBI."
                ),
            epilog="Semidán Robaina Estévez (srobaina@ull.edu.es), 2022",
            usage=("pynteny download [-h] ̣̣̣̣̣̣̣̣̣̣̣̣[args] ̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣̣\n"),
            formatter_class=argparse.RawTextHelpFormatter
            )

        optional = parser._action_groups.pop()
        required = parser.add_argument_group("required arguments")
        parser._action_groups.append(optional)

        optional.add_argument("--dir", dest="dir", type=Path,
                             required=False, default=None,
                             help="path to directory where to download HMM database"
                             )
        optional.add_argument("--unpack", dest="unpack",
                            default=False, action="store_true",
                            help="unpack originally compressed database files"
                            )
        args = parser.parse_args(self._subcommand_args)
        sub.download_hmms(args)
        


def main():
    pynteny = Pynteny()
    
if __name__ == "__main__":
    main()
