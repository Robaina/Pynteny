#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse

from pynteny.utils import parallelizeOverInputFiles, fullPathListDir, setDefaultOutputPath
from pynteny.wrappers import runProdigal
from pynteny.preprocessing import parseProdigalOutput, splitFASTAbyContigs, mergeFASTAs


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
                          "path to assembly input nucleotide data. It can be a single FASTA file or "
                          "a directory containing several FASTA files."
                          )
)
optional.add_argument('--outdir', dest='outdir', type=str,
                      help='path to output directory'
)
optional.add_argument('--prefix', dest='prefix', type=str,
                      default='prodigal',
                      help='prefix to be added to output files'
)
optional.add_argument("--processes", "-p", dest="processes", type=int,
                      required=False, default=None, help=(
                          "set maximum number of processes. "
                          "Defaults to all but one."
                          )
                          )
optional.add_argument('--split_contigs', dest='split',
                      default=False, action='store_true',
                      help='split assembly input file into files containing one contig each')
optional.add_argument('--metagenome', dest='metagenome',
                      default=False, action='store_true',
                      help='treat input assembly as of type metagenome (for prodigal)')
optional.add_argument("--prodigal_args", dest="prodigal_args", type=str,
                      default=None, required=False,
                      help=(("additional arguments to prodigal provided as defined by it"))
                    )
args = parser.parse_args()

if args.processes is None:
    args.processes = os.cpu_count() - 1

def main():


    if args.split:
        print("1. Splitting assembly file...")
        split_dir = os.path.join(
            setDefaultOutputPath(args.assembly_fasta, only_dirname=True),
            f"split_{setDefaultOutputPath(args.assembly_fasta, only_basename=True)}"
        )
        os.makedirs(split_dir)
        splitFASTAbyContigs(
            input_fasta=args.assembly_fasta,
            output_dir=split_dir
        )
        input_assembly = split_dir
    else:
        input_assembly = args.assembly_fasta

    print("2. Running prodigal on assembly...")
    if os.path.isdir(input_assembly):
        split_prodigal_dir = os.path.join(args.outdir, "split_prodigal/")
        os.makedirs(split_prodigal_dir)
        parallelizeOverInputFiles(
            runProdigal, 
            input_list=fullPathListDir(input_assembly),
            n_processes=args.processes,
            output_dir=split_prodigal_dir,
            metagenome=args.metagenome,
            additional_args=args.prodigal_args
        )
        mergeFASTAs(
            split_prodigal_dir,
            output_fasta=os.path.join(
                args.outdir,
                f"{args.prefix}merged.faa"
                )
        )
    else:
        runProdigal(input_assembly,
            input_file=input_assembly,
            output_prefix=args.prefix,
            output_dir=args.outdir,
            metagenome=True,
            additional_args=None
        )
    os.remove(split_dir)
    
    print("3. Parsing prodigal output...")
    parseProdigalOutput(
        prodigal_faa=os.path.join(args.outdir, f"{args.prefix}merged.faa"),
        output_file=os.path.join(args.outdir, f"{args.prefix}positioned.faa")
    )

    print("Finished!")


if __name__ == '__main__':
    main()