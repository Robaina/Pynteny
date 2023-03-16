#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions containing CLI subcommands
"""

import logging
import os
import shutil
import sys
from argparse import ArgumentParser
from pathlib import Path
from typing import Union

import numpy as np

import pynteny.parsers.syntenyparser as syntenyparser
from pynteny.filter import SyntenyHits, filter_FASTA_by_synteny_structure
from pynteny.hmm import PGAP
from pynteny.preprocessing import Database
from pynteny.utils import (
    CommandArgs,
    ConfigParser,
    download_file,
    is_tar_file,
)


def init_logger(args: Union[CommandArgs, ArgumentParser]) -> logging.Logger:
    """Initialize logger object

    Args:
        args (Union[CommandArgs, ArgumentParser]): arguments object

    Returns:
        logging.Logger: initialized logger object
    """
    if args.logfile is None:
        args.logfile = Path(os.devnull)
    elif not Path(args.logfile.parent).exists():
        Path(args.logfile.parent).mkdir(parents=True)
    logging.basicConfig(
        format="%(asctime)s | %(levelname)s: %(message)s",
        handlers=[logging.FileHandler(args.logfile), logging.StreamHandler(sys.stdout)],
        level=logging.NOTSET,
    )
    logger = logging.getLogger(__name__)
    return logger


def synteny_search(args: Union[CommandArgs, ArgumentParser]) -> SyntenyHits:
    """Search peptide database by synteny structure containing HMMs.

    Args:
        args (Union[CommandArgs, ArgumentParser]): arguments object.

    Returns:
        SyntenyHits: instance of SyntenyHits.
    """
    logger = init_logger(args)
    if not Path(args.data).exists():
        logger.error("Sequence data file does not exist")
        sys.exit(1)

    config = ConfigParser.get_default_config()
    args.synteny_struc = syntenyparser.reformat_synteny_structure(args.synteny_struc)
    if not syntenyparser.is_valid_structure(args.synteny_struc):
        logger.error(
            (
                f"Invalid synteny structure format: {args.synteny_struc}. "
                "Execute pynteny search --help to see the right format."
            )
        )
        sys.exit(1)
    if args.hmm_dir is None:
        if not config.get_field("data_downloaded"):
            logger.warning(
                "HMM database not found. Downloading PGAP database from NCBI"
            )
            down_args = CommandArgs(unpack=True, outdir=None, logfile=None)
            download_hmms(down_args)
        else:
            args.hmm_dir = Path(config.get_field("PGAP_database"))
    if args.gene_ids:
        if args.hmm_meta is None:
            if not config.get_field("data_downloaded"):
                logger.error(
                    "Please download hmm database first or provide path to hmm metadata file."
                )
                sys.exit(1)
            else:
                args.hmm_meta = Path(config.get_field("PGAP_meta_file"))
        logger.info("Finding matching HMMs for gene symbols")
        (
            gene_synteny_struc,
            gene_to_hmm_group,
        ) = syntenyparser.parse_genes_in_synteny_structure(
            synteny_structure=args.synteny_struc, hmm_meta=args.hmm_meta
        )
        args.synteny_struc = gene_synteny_struc
        logger.info(
            f"Found the following HMMs in database for given structure:\n{gene_synteny_struc}"
        )

    temp_hmm_dir = Path(args.hmm_dir.parent) / "temp_hmm_dir"
    if is_tar_file(args.hmm_dir):
        PGAP.extract_PGAP_to_directory(args.hmm_dir, temp_hmm_dir)
        hmm_dir = temp_hmm_dir
    else:
        hmm_dir = args.hmm_dir

    hmm_names = syntenyparser.get_all_HMMs_in_structure(args.synteny_struc)
    unique_hmm_names = np.unique(hmm_names)
    input_hmms = [
        file
        for file in hmm_dir.iterdir()
        if any([hmm_name in file.as_posix() for hmm_name in hmm_names])
    ]
    if len(input_hmms) < len(unique_hmm_names):
        logger.error(
            "Not all HMMs in synteny structure found in HMM directory. "
            "Remember to include '--gene_ids' option if you want to search by gene symbols."
        )
        sys.exit(1)
    if args.outdir is None:
        args.outdir = Path(args.data.parent)
    if not args.outdir.exists():
        args.outdir.mkdir(parents=True, exist_ok=True)
    if args.hmmsearch_args is None:
        hmmsearch_args = ",".join(["None" for _ in input_hmms])
    else:
        hmmsearch_args = args.hmmsearch_args
    hmmsearch_args = list(map(lambda x: x.strip(), hmmsearch_args.split(",")))
    hmmsearch_args = list(map(lambda x: None if x == "None" else x, hmmsearch_args))
    hmmer_output_dir = os.path.join(args.outdir, "hmmer_outputs/")
    synteny_table = args.outdir / f"{args.prefix}synteny_matched.tsv"

    logger.info("Searching database by synteny structure")
    synteny_hits = filter_FASTA_by_synteny_structure(
        synteny_structure=args.synteny_struc,
        unordered=args.unordered,
        input_fasta=args.data,
        input_hmms=input_hmms,
        hmm_meta=args.hmm_meta,
        hmmer_output_dir=hmmer_output_dir,
        reuse_hmmer_results=args.reuse,
        method="hmmsearch",
        processes=args.processes,
        additional_args=hmmsearch_args,
    )
    synteny_hits.write_to_TSV(synteny_table)
    logger.info("Writing matching sequences to FASTA files")
    synteny_hits.write_hit_sequences_to_FASTA_files(
        sequence_database=args.data, output_dir=args.outdir, output_prefix=args.prefix
    )
    if temp_hmm_dir.exists():
        shutil.rmtree(temp_hmm_dir)
    logger.info("Finished!")
    logging.shutdown()
    return synteny_hits


def build_database(args: Union[CommandArgs, ArgumentParser]) -> None:
    """Build annotated peptide database from input assembly
    or GenBank data.

    Args:
        args (Union[CommandArgs, ArgumentParser]): arguments object.
    """
    logger = init_logger(args)
    logger.info("Building annotated peptide database")
    if args.processes is None:
        args.processes = os.cpu_count() - 1
    if args.data.is_dir() and args.prepend:
        prepend_file_name = True
    else:
        prepend_file_name = False

    database = Database(args.data)
    database.build(
        output_file=args.outfile,
        prepend_file_name=prepend_file_name,
        processes=args.processes,
        tmpdir=args.tmpdir,
    )
    logger.info("Database built successfully!")
    logging.shutdown()


def parse_gene_ids(args: Union[CommandArgs, ArgumentParser]) -> str:
    """Convert gene symbols to hmm names.

    Args:
        args (Union[CommandArgs, ArgumentParser]): arguments object.

    Returns:
        str: synteny structure where gene symbols are replaced
            by HMM names.
    """
    logger = init_logger(args)
    config = ConfigParser.get_default_config()
    if args.hmm_meta is None:
        if not config.get_field("data_downloaded"):
            logger.error(
                "Please download hmm database meta file or provide path to existing one first."
            )
            sys.exit(1)
        else:
            args.hmm_meta = Path(config.get_field("PGAP_meta_file"))
    (
        gene_synteny_struc,
        gene_to_hmm_group,
    ) = syntenyparser.parse_genes_in_synteny_structure(
        synteny_structure=args.synteny_struc, hmm_meta=args.hmm_meta
    )
    logger.info(
        f'Translated \n "{args.synteny_struc}" \n to \n "{gene_synteny_struc}" \n according to provided HMM database metadata'
    )
    logging.shutdown()
    return gene_synteny_struc


def download_hmms(args: Union[CommandArgs, ArgumentParser]) -> None:
    """Download HMM (PGAP) database from NCBI.

    Args:
        args (Union[CommandArgs, ArgumentParser]): arguments object.
    """
    logger = init_logger(args)
    config = ConfigParser.get_default_config()
    if (config.get_field("data_downloaded")) and (not args.force):
        logger.info("PGAP database already downloaded. Skipping download")
        sys.exit(1)
    if args.outdir is None:
        logger.error("Please, provide directory in which to download PGAP database.")
        sys.exit(1)
    else:
        download_dir = Path(args.outdir).absolute()
    if not download_dir.exists():
        download_dir.mkdir(parents=True, exist_ok=True)

    config.update_config("database_dir", download_dir.as_posix())
    config.update_config("upack_PGAP_database", args.unpack)

    data_url = "https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz"
    meta_url = "https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv"
    logger.info("Downloading PGAP database")
    try:
        PGAP_file = download_dir / "hmm_PGAP.HMM.tgz"
        meta_file = download_dir / "hmm_PGAP.tsv"
        download_file(data_url, PGAP_file)
        download_file(meta_url, meta_file)
        logger.info("Database dowloaded successfully\n")
        config.update_config("data_downloaded", True)
        config.update_config("PGAP_database", PGAP_file.as_posix())
        config.update_config("PGAP_meta_file", meta_file.as_posix())
    except Exception:
        logger.exception(
            "Failed to download PGAP database. Please check your internet connection."
        )
        sys.exit(1)
    logger.info("Removing missing entries from PGAP metadata file")
    PGAP(meta_file).remove_missing_HMMs_from_metadata(PGAP_file, meta_file)
    if args.unpack:
        logger.info("Unpacking PGAP database")
        unpacked_PGAP_dir = download_dir / "hmm_PGAP"
        PGAP.extract_PGAP_to_directory(PGAP_file, output_dir=unpacked_PGAP_dir)
        os.remove(PGAP_file)
        config.update_config("PGAP_database", unpacked_PGAP_dir.as_posix())
        logger.info("PGAP database unpacked successfully")
    logging.shutdown()


def get_citation(args: Union[CommandArgs, ArgumentParser], silent: bool = False) -> str:
    """Get Pynteny citation string.

    Args:
        args (argparse.ArgumentParser): arguments object.
        silent (bool, optional): do not print to terminal.
            Defaults to False.

    Returns:
        str: Pyntey citation text.
    """
    citation = (
        "Semidán Robaina Estévez (2022). Pynteny: synteny-aware hmm searches made easy"
        f"(Version {args.version}). Zenodo. https://doi.org/10.5281/zenodo.7048685"
    )
    if not silent:
        print("If you use this software, please cite it as below: ")
        print(citation)
    return citation
