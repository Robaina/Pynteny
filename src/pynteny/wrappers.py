#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
In-process replacements for the external CLI tools that Pynteny used to shell
out to (HMMER, Prodigal, seqkit). These now rely on the pure-Python,
pip-installable packages pyhmmer, pyrodigal and pyfastx so that Pynteny can be
installed without conda.
"""

from __future__ import annotations

import logging
import os
import shlex
from pathlib import Path

import pyfastx
import pyhmmer

from pynteny.utils import set_default_output_path

logger = logging.getLogger(__name__)


def run_seqkit_nodup(
    input_fasta: Path, output_fasta: Path = None, export_duplicates: bool = False
):
    """Remove duplicated sequences from a fasta file (by sequence content).

    Pure-Python replacement for ``seqkit rmdup -s``.

    Args:
        input_fasta (Path): path to input fasta.
        output_fasta (Path, optional): path to output fasta. Defaults to None.
        export_duplicates (bool, optional): whether to export a file containing
            the IDs of duplicated sequences. Defaults to False.
    """
    input_fasta = Path(input_fasta)
    if output_fasta is None:
        output_fasta = set_default_output_path(input_fasta, tag="_no_duplicates")
    else:
        output_fasta = Path(output_fasta)

    seen_sequences = set()
    duplicated_ids = []
    fasta = pyfastx.Fasta(input_fasta.as_posix(), build_index=False, full_name=True)
    with open(output_fasta, "w+", encoding="UTF-8") as outfile:
        for record_name, record_seq in fasta:
            seq_key = record_seq.upper()
            if seq_key in seen_sequences:
                duplicated_ids.append(record_name.split()[0])
                continue
            seen_sequences.add(seq_key)
            outfile.write(f">{record_name}\n{record_seq}\n")

    if export_duplicates:
        dup_file = set_default_output_path(
            input_fasta, tag="_duplicates", extension=".txt"
        )
        with open(dup_file, "w+", encoding="UTF-8") as f:
            f.write("\n".join(duplicated_ids))


def _parse_hmmsearch_args(additional_args: str) -> dict:
    """Translate a (subset of) hmmsearch/hmmscan CLI argument strings into
    keyword options understood by :func:`pyhmmer.hmmsearch`.

    Only the options commonly used with Pynteny are supported. Unknown flags
    (and ``--cpu``, which is handled separately) are ignored with a warning.

    Args:
        additional_args (str): CLI-style argument string, e.g. "--cut_ga" or
            "-E 1e-10".

    Returns:
        dict: keyword options for the pyhmmer Pipeline.
    """
    if not additional_args:
        return {}
    bit_cutoff_flags = {
        "--cut_ga": "gathering",
        "--cut_nc": "noise",
        "--cut_tc": "trusted",
    }
    float_flags = {
        "-E": "E",
        "--incE": "incE",
        "--domE": "domE",
        "-T": "T",
        "--incT": "incT",
        "--domT": "domT",
        "-Z": "Z",
        "--domZ": "domZ",
    }
    options = {}
    tokens = shlex.split(additional_args)
    i = 0
    while i < len(tokens):
        token = tokens[i]
        if token in bit_cutoff_flags:
            options["bit_cutoffs"] = bit_cutoff_flags[token]
            i += 1
        elif token in float_flags:
            if i + 1 >= len(tokens):
                logger.warning(f"Missing value for hmmsearch argument '{token}'")
                break
            options[float_flags[token]] = float(tokens[i + 1])
            i += 2
        elif token in ("--cpu",):
            # number of threads is handled through the 'processes' argument
            i += 2
        else:
            logger.warning(
                f"Ignoring unsupported hmmsearch argument for pyhmmer: '{token}'"
            )
            i += 1
    return options


def run_HMM_search(
    hmm_model: Path,
    input_fasta: Path,
    output_file: Path = None,
    method: str = "hmmsearch",
    processes: int = None,
    additional_args: str = None,
) -> None:
    """Run an HMM search with pyhmmer and write the results as an HMMER3
    tabular (``--tblout``) output file.

    In-process replacement for the ``hmmsearch``/``hmmscan`` CLI tools. The
    on-disk format is unchanged (standard HMMER3 tabular output), so downstream
    parsing and result reuse keep working as before.

    Args:
        hmm_model (Path): path to profile HMM to be used.
        input_fasta (Path): path to fasta containing the sequence database to be searched.
        output_file (Path, optional): path to HMMER tabular output file. Defaults to None.
        method (str, optional): either 'hmmsearch' or 'hmmscan'. Defaults to 'hmmsearch'.
        processes (int, optional): maximum number of threads. Defaults to all minus one.
        additional_args (str, optional): a CLI-style string with additional
            arguments for hmmsearch/hmmscan. Defaults to None.
    """
    hmm_model = Path(hmm_model)
    input_fasta = Path(input_fasta)
    if processes is None:
        processes = os.cpu_count() - 1
    if output_file is None:
        output_file = set_default_output_path(input_fasta, "_hmmer_hits", ".txt")
    else:
        output_file = Path(output_file)

    options = _parse_hmmsearch_args(additional_args)

    with pyhmmer.easel.SequenceFile(input_fasta.as_posix(), digital=True) as seq_file:
        sequences = seq_file.read_block()
    with pyhmmer.plan7.HMMFile(hmm_model.as_posix()) as hmm_file:
        hmms = list(hmm_file)

    search = pyhmmer.hmmscan if method == "hmmscan" else pyhmmer.hmmsearch
    if method == "hmmscan":
        # hmmscan queries sequences against the HMM database
        all_hits = list(search(sequences, hmms, cpus=processes, **options))
    else:
        all_hits = list(search(hmms, sequences, cpus=processes, **options))

    with open(output_file, "wb") as fh:
        for i, hits in enumerate(all_hits):
            hits.write(fh, format="targets", header=(i == 0))
