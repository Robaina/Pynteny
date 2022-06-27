#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Simple CLI wrappers to several tools
"""

import os
from pathlib import Path

from pynteny.pynteny.utils import terminalExecute, setDefaultOutputPath



def runSeqKitNoDup(input_fasta: Path, output_fasta: Path = None,
                   export_duplicates: bool = False):
    """
    Simpe CLI wrapper to seqkit rmdup
    """
    if output_fasta is None:
        output_fasta = setDefaultOutputPath(input_fasta, tag="_no_duplicates")
    if export_duplicates:
        dup_file = setDefaultOutputPath(input_fasta, tag="_duplicates", extension=".txt")
        dup_str = f"-D {dup_file}"
    else:
        dup_str = ""
    cmd_str = (
        f"seqkit rmdup {input_fasta} -s {dup_str} -o {output_fasta.as_posix()}"
    )
    terminalExecute(cmd_str)

def runProdigal(input_file: Path,
                output_prefix: str = None,
                output_dir: Path = None,
                metagenome: bool = False,
                additional_args: str = None):
    """
    Simple CLI wrapper to prodigal
    """
    if metagenome:
        procedure = 'meta'
    else:
        procedure = 'single'
    if output_dir is None:
        output_dir = setDefaultOutputPath(input_file, only_dirname=True)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    if output_prefix is None:
        output_prefix = setDefaultOutputPath(input_file, only_basename=True)
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ''
    output_gbk = output_dir / f"{output_prefix}.gbk"
    output_fasta = output_dir / f"{output_prefix}.faa"
    cmd_str = (
        f'prodigal -i {input_file} -o {output_gbk} -p {procedure} '
        f'-a {output_fasta} -q {args_str}'
        )
    terminalExecute(cmd_str, suppress_shell_output=False)

def runHMMsearch(hmm_model: Path, input_fasta: Path,
                 output_file: Path = None,
                 method: str = 'hmmsearch',
                 n_processes: int = None,
                 additional_args: str = None) -> None:
    """
    Simple CLI wrapper to hmmsearch or hmmscan
    --cut_nc, --cut_ga
    """
    if n_processes is None:
        n_processes = os.cpu_count() - 1
    if output_file is None:
        output_file = setDefaultOutputPath(input_fasta, '_hmmer_hits', '.txt')
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ''
    cmd_str = (
        f'{method} --tblout {output_file} {args_str} --cpu {n_processes} '
        f'{hmm_model} {input_fasta}'
        )
    terminalExecute(cmd_str, suppress_shell_output=True)

def runHMMbuild(input_aln: Path, output_hmm: Path = None,
                additional_args: str = None) -> None:
    """
    Simple CLI wrapper to hmmbuild (build HMM profile from MSA file)
    additional args: see hmmbuild -h
    """
    if output_hmm is None:
        output_hmm = setDefaultOutputPath(input_aln, extension='.hmm')
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ''
    cmd_str = f'hmmbuild {args_str} {output_hmm} {input_aln}'
    terminalExecute(cmd_str, suppress_shell_output=False)