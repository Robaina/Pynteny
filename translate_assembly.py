#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse

from pynteny.utils import parallelizeOverInputFiles, fullPathListDir
from pynteny.wrappers import runProdigal
from pynteny.preprocessing import parseProdigalOutput, mergeFASTAs


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
optional.add_argument("--processes", "-p", dest="processes", type=int,
                      required=False, default=None, help=(
                          "set maximum number of processes. "
                          "Defaults to all but one."
                          )
                          )

args = parser.parse_args()

if args.processes is None:
    args.processes = os.cpu_count() - 1

def main():


    if os.path.isdir(args.assembly_fasta):
        parallelizeOverInputFiles(
            runProdigal, 
            input_list=fullPathListDir(args.assembly_fasta),
            n_processes=args.processes,
            output_dir=args.outdir,
            metagenome=True,
            additional_args=None
        )
        mergeFASTAs(
            args.outdir,
            output_fasta=os.path.join()
        )
    else:
        runProdigal(
            input_file=args.assembly_fasta,
            output_prefix=args.prefix,
            output_dir=args.outdir,
            metagenome=True,
            additional_args=None
        )

    # assignGeneLocationToRecords(
    #     gbk_file=os.path.join(args.outdir, f"{args.prefix}.gbk"),
    #     output_fasta=os.path.join(args.outdir, "ref_database.faa"),
    #     nucleotide=False
    # )

    # parseProdigalOutput(
    #     prodigal_faa=os.path.join(args.outdir, f"{args.prefix}.faa"),
    #     output_file=None
    # )


if __name__ == '__main__':
    main()