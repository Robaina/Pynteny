#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

from pynteny.src.utils import setDefaultOutputPath, isTarFile, extractTarFile, flattenDirectory
from pynteny.src.filter import filterFASTAbySyntenyStructure, SyntenyParser
import pynteny.src.subcommands as sub


parser = argparse.ArgumentParser(
    description=(
        'Build peptide database from HMM synteny structure. '
        'The script outputs a main file containing sequences matching the provided '
        'hmm structure and corresponding to the main target indicated in the argument: '
        'target_hmm. These sequences are filtered by sequence length and by maximum '
        'number of total sequences. '
        'The script also outputs additional file containing the matched records for the '
        'other, non-target, hmms. However, these files are not processed any further.'
        ),
    epilog='Semidán Robaina Estévez (srobaina@ull.edu.es), 2022'
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

args = parser.parse_args()

if __name__ == '__main__':
    sub.synteny_search(args)
