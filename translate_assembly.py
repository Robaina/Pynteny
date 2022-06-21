#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil
import argparse
from pathlib import Path

from pynteny.utils import TemporaryFilePath, parallelizeOverInputFiles, fullPathListDir, setDefaultOutputPath
from pynteny.wrappers import runProdigal
from pynteny.preprocessing import FASTA, LabelledFASTA
    # (parseProdigalOutput, splitFASTAbyContigs,
    # mergeFASTAs, removeCorruptedSequences
    # )


parser = argparse.ArgumentParser(
    description=(
        "Translate nucleotide assembly file and assign contig and gene location info "
        "to each identified ORF (using prodigal). Then label predicted ORFs according to "
        "positional info and export a fasta file containing predicted and translated ORFs."
        ),
    epilog="Semidán Robaina Estévez (srobaina@ull.edu.es), 2022"
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
        split_dir = Path(args.assembly_fasta.parent) / f"split_{args.assembly_fasta.stem}"
        os.makedirs(split_dir)
        FASTA(args.assembly_fasta).splitByContigs(
            output_dir=split_dir
        )
        input_assembly = split_dir
    else:
        input_assembly = args.assembly_fasta

    print("2. Running prodigal on assembly...")
    if input_assembly.is_dir():
        split_prodigal_dir = args.outdir / "split_prodigal/"
        os.makedirs(split_prodigal_dir)
        parallelizeOverInputFiles(
            runProdigal, 
            input_list=list(input_assembly.iterdir()),
            n_processes=args.processes,
            output_dir=split_prodigal_dir,
            metagenome=args.metagenome,
            additional_args=args.prodigal_args
        )
        FASTA.mergeFASTAs(
            split_prodigal_dir,
            output_file=args.outdir / f"{args.prefix}merged.faa"
        )
    else:
        runProdigal(input_assembly,
            input_file=input_assembly,
            output_prefix=args.prefix,
            output_dir=args.outdir,
            metagenome=True,
            additional_args=None
        )
    shutil.rmtree(split_dir)
    
    print("3. Parsing prodigal output...")
    with TemporaryFilePath() as tempfasta:
        labelledfasta = LabelledFASTA.fromProdigalOutput(
            prodigal_faa=args.outdir / f"{args.prefix}merged.faa",
            output_file=Path(tempfasta)
        )
        labelledfasta.removeCorruptedSequences(
            output_file=args.outdir / f"{args.prefix}positioned.faa",
            is_peptide=True,
            keep_stop_codon=True
        )

    print("Finished!")


if __name__ == "__main__":
    main()