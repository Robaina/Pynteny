#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
import argparse
from pathlib import Path

import pynteny.subcommands as sub



class Pynteny():
    """
    Temporary fix to enable subcommands without too much
    code rewriting. Based on:
    https://selvakumar-arumugapandian.medium.com/ \ 
    command-line-subcommands-with-pythons-argparse-4dbac80f7110

    """
    def __init__(self):
        self._subcommand = sys.argv[1:2]
        self._subcommand_args = sys.argv[2:]
        parser = argparse.ArgumentParser(
            description=("Pynteny: synteny-based hmm searches made easy in Python"),
            epilog="Semidán Robaina Estévez (srobaina@ull.edu.es), 2022",
            formatter_class=argparse.RawDescriptionHelpFormatter
            )
        parser.add_argument("subcommand", help="pynteny subcommand", choices=["search", "preprocess"])
        parser.add_argument("-v","--version", help="show version and exit", action="version", version="0.0.1")
        args = parser.parse_args(self._subcommand)
        self._call_subcommand(subcommand_name=args.subcommand)
    
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
            epilog='Semidán Robaina Estévez (srobaina@ull.edu.es), 2022',
            formatter_class=argparse.RawDescriptionHelpFormatter
            )

        optional = parser._action_groups.pop()
        required = parser.add_argument_group('required arguments')
        parser._action_groups.append(optional)

        required.add_argument('--hmm_dir', dest='hmm_dir', type=Path,
                            required=True,
                            help=(
                                'path to directory containing hmm (i.e, tigrfam or pfam) models. \n'
                                'The directory can contain more hmm models than used in the synteny structure. '
                                'It may also be the path to a compressed (tar, tar.gz, tgz) directory.'
                                )
        )
        required.add_argument('--synteny_struc', dest='synteny_struc', type=str, required=True,
                            help=(
                                'string displaying hmm sctructure to search for, such as: \n'
                                '">hmm_a n_ab <hmm_b n_bc hmm_c", \n'
                                'where ">" indicates a hmm target located on the positive strand, '
                                '"<" a target located on the negative strand, and n_ab cooresponds '
                                'to the maximum number of genes separating matched gene a and b. \n' 
                                'Multiple hmms may be employed (limited by computational capabilities).'
                                'No order symbol in a hmm indicates that results should be independent '
                                'of strand location. '
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
                                'list of comma-separated additional arguments to hmmsearch for each input hmm. '
                                'A single argument may be provided, in which case the same additional argument '
                                'is employed in all hmms.')
        )
        parser.add_argument('--max_seq_length', dest='maxseqlength',
                            default=None, type=int,
                            required=False,
                            help=(
                                'maximum sequence length in reference database. '
                                'Defaults to inf'
                                )
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
            epilog="Semidán Robaina Estévez (srobaina@ull.edu.es), 2022",
            formatter_class=argparse.RawDescriptionHelpFormatter
            )

        optional = parser._action_groups.pop()
        required = parser.add_argument_group("required arguments")
        parser._action_groups.append(optional)

        required.add_argument("--assembly_fasta", dest="assembly_fasta", type=Path,
                            required=True,
                            help=(
                                "path to assembly input nucleotide data. It can be a single FASTA file or "
                                "a directory containing several FASTA files."
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


def main():
    Pynteny()
    
if __name__ == "__main__":
    Pynteny()
