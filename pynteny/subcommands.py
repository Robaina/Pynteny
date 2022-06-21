#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions containing CLI subcommands
"""

import os
import shutil
from pathlib import Path

from pynteny.utils import setDefaultOutputPath, isTarFile, extractTarFile, flattenDirectory
from pynteny.filter import filterFASTABySyntenyStructure, SyntenyParser

from pynteny.utils import TemporaryFilePath, parallelizeOverInputFiles, setDefaultOutputPath
from pynteny.wrappers import runProdigal
from pynteny.preprocessing import FASTA, LabelledFASTA



def synteny_seach(args):

    if isTarFile(args.hmm_dir):
        print("0. Extracting hmm files to temporary directory...")
        temp_hmm_dir = Path(args.hmm_dir.parent) / "temp_hmm_dir"
        extractTarFile(
            tar_file=args.hmm_dir,
            dest_dir=temp_hmm_dir
        )
        flattenDirectory(
            temp_hmm_dir
        )
        hmm_dir = temp_hmm_dir
    else:
        hmm_dir = args.hmm_dir

    hmm_names = SyntenyParser.getHMMsInStructure(args.synteny_struc)
    input_hmms = [
        Path(os.path.join(hmm_dir, file))
        for file in os.listdir(hmm_dir)
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

    if isTarFile(args.hmm_dir):
        shutil.rmtree(temp_hmm_dir)

    print('Finished!')


def translate_assembly(args):

    if args.processes is None:
        args.processes = os.cpu_count() - 1

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