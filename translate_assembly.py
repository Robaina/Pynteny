#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil
import argparse
from pathlib import Path

from pynteny.utils import TemporaryFilePath, parallelizeOverInputFiles, fullPathListDir, setDefaultOutputPath
from pynteny.wrappers import runProdigal
from pynteny.preprocessing import FASTA, LabelledFASTA
import pynteny.subcommands as sub


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

if __name__ == "__main__":
    sub.translate_assembly(args)