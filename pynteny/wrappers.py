#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Simple CLI wrappers to several tools
"""

import os

from pynteny.utils import terminalExecute, setDefaultOutputPath


def runSeqKitNoDup(input_fasta: str, output_fasta: str = None,
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
        f"seqkit rmdup {input_fasta} -s {dup_str} -o {output_fasta}"
    )
    terminalExecute(cmd_str)

def runFastP(input_fastq1: str, input_fastq2: str = None,
             merge: bool = False, merged_out: str = None,
             output_dir: str = None,
             n_threads: int = None,
             additional_args: str = None) -> None:
    """
    Simple CLI wrapper to Fastp
    https://github.com/OpenGene/fastp
    """
    if input_fastq2 is None:
        input_str = f'-i {input_fastq1}'
    else:
        input_str = f'-i {input_fastq1} -I {input_fastq2}'
    if n_threads is not None:
        threads_str = f'--thread {n_threads}'
    else:
        threads_str = ''
    if output_dir is None:
        output_dir = setDefaultOutputPath(input_fastq1, only_dirname=True)
    if merge:
        merged_out = merged_out if merged_out is not None else os.path.join(output_dir, 'merged.fastq')
        out_str = f'--merge --merged_out {merged_out}'
    else:
        out_str = f'-o {"fastp_" + input_fastq1} -O {"fastp_" + input_fastq2}'
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ''
    cmd_str = f'fastp {input_str} {out_str} {threads_str} {args_str}'
    terminalExecute(cmd_str, suppress_shell_output=True)
        
def runProdigal(input_file: str, output_prefix: str = None,
                output_dir: str = None,
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
    output_gbk = os.path.join(output_dir, output_prefix + '.gbk')
    output_fasta = os.path.join(output_dir, output_prefix + '.faa')
    cmd_str = (
        f'prodigal -i {input_file} -o {output_gbk} -p {procedure} '
        f'-a {output_fasta} -q {args_str}'
        )
    terminalExecute(cmd_str, suppress_shell_output=False)

def runHMMsearch(hmm_model: str, input_fasta: str,
                 output_file: str = None,
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

def runHMMbuild(input_aln: str, output_hmm: str = None,
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
