#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Simple CLI wrappers to several tools
"""

import os
from pathlib import Path

from pynteny.utils import set_default_output_path, terminal_execute


def run_seqkit_nodup(
    input_fasta: Path, output_fasta: Path = None, export_duplicates: bool = False
):
    """Simpe CLI wrapper to seqkit rmdup to remove sequence duplicates
    in fasta file.

    Args:
        input_fasta (Path): path to input fasta.
        output_fasta (Path, optional): path to output fasta. Defaults to None.
        export_duplicates (bool, optional): whether to export a file containing
            duplicated sequences. Defaults to False.
    """
    input_fasta = Path(input_fasta)
    if output_fasta is None:
        output_fasta = set_default_output_path(input_fasta, tag="_no_duplicates")
    else:
        output_fasta = Path(output_fasta)
    if export_duplicates:
        dup_file = set_default_output_path(
            input_fasta, tag="_duplicates", extension=".txt"
        )
        dup_str = f"-D {dup_file}"
    else:
        dup_str = ""
    cmd_str = (
        f"seqkit rmdup {input_fasta} -s {dup_str} -o {output_fasta.as_posix()} --quiet"
    )
    terminal_execute(cmd_str, suppress_shell_output=True)


def run_prodigal(
    input_file: Path,
    output_file: Path = None,
    output_dir: Path = None,
    output_format: str = "fasta",
    metagenome: bool = False,
    additional_args: str = None,
):
    """Simple CLI wrapper to prodigal.

    Args:
        input_file (Path): path to input fasta file with nucleotide sequences.
        output_file (Path, optional): path to output file containing translated peptides.
            Defaults to None.
        output_dir (Path, optional): path to output directory (all prodigal output files).
            Defaults to None.
        output_format (str, optional): either 'gbk' or 'fasta'. Defaults to 'fasta'.
        metagenome (bool, optional): whether input fasta correspond to a metagenomic sample.
            Defaults to False.
        additional_args (str, optional): a string containing additional arguments to prodigal.
            Defaults to None.
    """
    input_file = Path(input_file)
    if metagenome:
        procedure = "meta"
    else:
        procedure = "single"
    if output_dir is None:
        output_dir = Path(input_file.parent)
    else:
        output_dir = Path(output_dir)
    if "fasta" in output_format.lower():
        if output_file is None:
            output_file = output_dir / f"{input_file.stem}prodigal_output.faa"
        out_str = f"-a {output_file}"
    else:
        if output_file is None:
            output_file = output_dir / f"{input_file.stem}prodigal_output.gbk"
        out_str = f"-o {output_file}"
    output_file = Path(output_file)
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ""
    cmd_str = f"prodigal -i {input_file} -p {procedure} " f"-q {out_str} {args_str}"
    terminal_execute(cmd_str, suppress_shell_output=True)


def run_HMM_search(
    hmm_model: Path,
    input_fasta: Path,
    output_file: Path = None,
    method: str = "hmmsearch",
    n_processes: int = None,
    additional_args: str = None,
) -> None:
    """Simple CLI wrapper to hmmsearch or hmmscan.

    Args:
        hmm_model (Path): path to profile HMM to be used.
        input_fasta (Path): path to fasta containing sequence database to be searched.
        output_file (Path, optional): path to prodigal output table file. Defaults to None.
        method (str, optional): either 'hmmsearch' or 'hmmscan'. Defaults to 'hmmsearch'.
        n_processes (int, optional): maximum number of threads. Defaults to all minus one.
        additional_args (str, optional): a string containing additional arguments to
            hmmsearch/scan. Defaults to None.
    """
    if n_processes is None:
        n_processes = os.cpu_count() - 1
    if output_file is None:
        output_file = set_default_output_path(input_fasta, "_hmmer_hits", ".txt")
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ""
    cmd_str = (
        f"{method} --tblout {output_file} {args_str} --cpu {n_processes} "
        f"{hmm_model} {input_fasta}"
    )
    terminal_execute(cmd_str, suppress_shell_output=True)
