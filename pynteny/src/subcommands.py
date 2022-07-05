#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions containing CLI subcommands
"""

import os
import sys
import shutil
import logging
from pathlib import Path
import wget

from pynteny.src.utils import ConfigParser, setDefaultOutputPath, isTarFile, extractTarFile, flattenDirectory
from pynteny.src.filter import filterFASTABySyntenyStructure, SyntenyParser

from pynteny.src.utils import TemporaryFilePath, parallelizeOverInputFiles, setDefaultOutputPath
from pynteny.src.wrappers import runProdigal
from pynteny.src.preprocessing import FASTA, LabelledFASTA


logging.basicConfig(format='%(asctime)s | %(levelname)s: %(message)s', level=logging.NOTSET)
logger = logging.getLogger(__name__)


def synteny_search(args):
    """
    Search peptide database by synteny structure containing HMMs.
    """
    config = ConfigParser.get_default_config()
    if args.hmm_dir is None:
        if not config.get_field("data_downloaded"):
            logger.error("Please download hmm database first or provide path to hmm directory.")
            sys.exit(1)
        else:
            args.hmm_dir = Path(config.get_field("PGAP_file"))
    if args.gene_ids:
        if args.hmm_meta is None:
            if not config.get_field("data_downloaded"):
                logger.error("Please download hmm database first or provide path to hmm metadata file.")
                sys.exit(1)
            else:
                args.hmm_meta = Path(config.get_field("PGAP_meta_file"))
        logger.info("Finding matching HMMs for gene symbols")
        # pgap = PGAP(args.hmm_meta)
        gene_synteny_struc = SyntenyParser.parseGenesInSyntenyStructure(
            synteny_structure=args.synteny_struc,
            hmm_meta=args.hmm_meta
        )
        args.synteny_struc = gene_synteny_struc
        logger.info(f"Found the following HMMs in database for given structure:\n{gene_synteny_struc}")

    if isTarFile(args.hmm_dir):
        logger.info("Extracting hmm files to temporary directory")
        temp_hmm_dir = Path(args.hmm_dir.parent) / "temp_hmm_dir"
        if temp_hmm_dir.exists():
            shutil.rmtree(temp_hmm_dir)
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

    hmm_names = SyntenyParser.getAllHMMsInStructure(args.synteny_struc)
    input_hmms = [
        file for file in hmm_dir.iterdir()
        if any([hmm_name in file.as_posix() for hmm_name in hmm_names])
    ]
    if len(input_hmms) < len(hmm_names):
        logger.error(
            "Not all HMMs in synteny structure found in HMM directory. "
            "Remember to include '--gene_ids' option if you want to search by gene symbols."
            )
        sys.exit(1)

    if args.outdir is None:
        args.outdir = setDefaultOutputPath(args.data, only_dirname=True)
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir, exist_ok=True)
    if args.hmmsearch_args is None:
        hmmsearch_args = ",".join(["None" for _ in input_hmms])
    hmmsearch_args = list(map(lambda x: x.strip(), hmmsearch_args.split(",")))
    hmmsearch_args = list(map(lambda x: None if x == 'None' else x, hmmsearch_args))
    hmmer_output_dir = os.path.join(args.outdir, 'hmmer_outputs/')
        
    logger.info('Searching database by synteny structure')
    filterFASTABySyntenyStructure(
        synteny_structure=args.synteny_struc,
        input_fasta=args.data,
        input_hmms=input_hmms,
        hmm_meta=args.hmm_meta,
        output_dir=args.outdir,
        output_prefix=args.prefix,
        hmmer_output_dir=hmmer_output_dir,
        reuse_hmmer_results=True,
        method='hmmsearch',
        additional_args=hmmsearch_args
    )
    if temp_hmm_dir.exists():
        shutil.rmtree(temp_hmm_dir)
    logger.info('Finished!')

def translate_assembly(args):
    """
    Preprocess assembly FASTA file and translate it to protein FASTA file.
    Add positional information to sequence headers.
    """
    if args.processes is None:
        args.processes = os.cpu_count() - 1

    if args.split:
        logger.info("Splitting assembly file")
        split_dir = os.path.join(
            setDefaultOutputPath(args.assembly_fasta, only_dirname=True),
            f"split_{setDefaultOutputPath(args.assembly_fasta, only_basename=True)}"
        )
        split_dir = Path(args.assembly_fasta.parent) / f"split_{args.assembly_fasta.stem}"
        os.makedirs(split_dir, exist_ok=True)
        FASTA(args.assembly_fasta).splitByContigs(
            output_dir=split_dir
        )
        input_assembly = split_dir
    else:
        input_assembly = args.assembly_fasta

    logger.info("Running prodigal on assembly")
    if input_assembly.is_dir():
        split_prodigal_dir = args.outdir / "split_prodigal/"
        os.makedirs(split_prodigal_dir, exist_ok=True)
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
    
    logger.info("Parsing prodigal output")
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
    logger.info("Finished!")

def parse_gene_ids(args):
    """
    Convert gene symbols to hmm names.
    """
    gene_synteny_struc = SyntenyParser.parseGenesInSyntenyStructure(
        synteny_structure=args.synteny_struc,
        hmm_meta=args.hmm_meta
        )
    logger.info(f"Found the following HMMs in database for given structure:\{gene_synteny_struc}")

def download_hmms(args):
    """
    Download HMM (PGAP) database from NCBI.
    """
    module_dir = Path(__file__).parent
    config_path = Path(module_dir.parent) / "config.json"
    config = ConfigParser(config_path)
    if args.dir is None:
        download_dir = Path(module_dir.parent) / "data"
    else:
        download_dir = Path(args.dir).absolute()
    if not download_dir.exists():
        os.makedirs(download_dir, exist_ok=True)

    if not config.get_field("data_downloaded"):
        config.update_config("database_dir", download_dir.as_posix())
        config.update_config("upack_PGAP_database", args.unpack)

        data_url = "https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz"
        meta_url = "https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv"
        logger.info("Downloading PGAP database")
        try:
            PGAP_file = download_dir / "hmm_PGAP.HMM.tgz"
            meta_file = download_dir / "hmm_PGAP.tsv"
            wget.download(data_url, PGAP_file.as_posix())
            wget.download(meta_url, meta_file.as_posix())
            print(" \n")
            logger.info("Database dowloaded successfully\n")
            config.update_config("data_downloaded", True)
            config.update_config("PGAP_file", PGAP_file.as_posix())
            config.update_config("PGAP_meta_file",  meta_file.as_posix())
        except Exception as e:
            logger.exception("Failed to download PGAP database. Please check your internet connection.")
            sys.exit(1)
