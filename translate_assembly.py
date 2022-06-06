#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse

from pynteny.wrappers import runProdigal
from pynteny.preprocessing import assignGeneLocationToRecords


parser = argparse.ArgumentParser(
    description=(
        "Translate nucleotide assembly file and assign contig and gene location info "
        "to each identified ORF (using prodigal). Then label predicted ORFs according to "
        "positional info and export a fasta file containing predicted and translated ORFs."
        ),
    epilog='Semidán Robaina Estévez (srobaina@ull.edu.es), 2022'
    )

optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
parser._action_groups.append(optional)

required.add_argument('--assembly_fasta', dest='assembly_fasta', type=str,
                      required=True,
                      help=(
                          'path to directory containing hmm (i.e, tigrfam or pfam) models. '
                          'The directory can contain more hmm models than used in the synteny structure.'
                          )
)
optional.add_argument('--outdir', dest='outdir', type=str,
                      help='path to output directory'
)
optional.add_argument('--prefix', dest='prefix', type=str,
                      default='',
                      help='prefix to be added to output files'
)



args = parser.parse_args()

def main():
    
    runProdigal(
        input_file=args.assembly_fasta,
        output_prefix=args.prefix,
        output_dir=args.outdir,
        metagenome=False,
        additional_args=None
    )

    assignGeneLocationToRecords(
        gbk_file=os.path.join(args.outdir, f"{args.prefix}.gbk"),
        output_fasta=os.path.join(args.outdir, "ref_database.faa"),
        nucleotide=False
    )


if __name__ == '__main__':
    main()