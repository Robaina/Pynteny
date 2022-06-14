#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse

from pynteny.utils import setDefaultOutputPath
from pynteny.filter import filterFASTABySyntenyStructure, SyntenyParser


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

required.add_argument('--hmm_dir', dest='hmm_dir', type=str,
                      required=True,
                      help=(
                          'path to directory containing hmm (i.e, tigrfam or pfam) models. '
                          'The directory can contain more hmm models than used in the synteny structure.'
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
required.add_argument('--in', dest='data', type=str, required=True,
                      help='path to peptide database'
)
optional.add_argument('--outdir', dest='outdir', type=str,
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
hmm_names = SyntenyParser.getHMMsInStructure(args.synteny_struc)
input_hmms = [
    os.path.join(args.hmm_dir, file)
    for file in os.listdir(args.hmm_dir)
    if any([hmm_name in file for hmm_name in hmm_names])
]
if len(input_hmms) < len(hmm_names):
    raise ValueError("Not all HMMs in synteny structure found in HMM directory")

if args.outdir is None:
    args.outdir = setDefaultOutputPath(args.data, only_dirname=True)
if not os.path.isdir(args.outdir):
    os.makedirs(args.outdir, exist_ok=True)
if args.hmmsearch_args is None:
    hmmsearch_args = ",".join(["None" for _ in input_hmms])
hmmsearch_args = list(map(lambda x: x.strip(), hmmsearch_args.split(",")))
hmmsearch_args = list(map(lambda x: None if x == 'None' else x, hmmsearch_args))
hmmer_output_dir = os.path.join(args.outdir, 'hmmer_outputs/')
    

def main():
    print(' 1. Searching database by synteny structure...')
    filterFASTABySyntenyStructure(
        synteny_structure=args.synteny_struc,
        input_fasta=args.data,
        input_hmms=input_hmms,
        output_dir=args.outdir,
        output_prefix=args.prefix,
        hmmer_output_dir=hmmer_output_dir,
        reuse_hmmer_results=True,
        method='hmmsearch',
        additional_args=hmmsearch_args
    )
    print('Finished!')

if __name__ == '__main__':
    main()
