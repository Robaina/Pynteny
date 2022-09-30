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

from pynteny.src.filter import SyntenyHits, SyntenyParser, filterFASTAbySyntenyStructure
from pynteny.src.hmm import PGAP
from pynteny.src.utils import CommandArgs, ConfigParser, isTarFile, terminalExecute
from pynteny.src.preprocessing import Database



def synteny_search(args) -> SyntenyHits:
    """
    Search peptide database by synteny structure containing HMMs.
    """
    logger = logging.getLogger(__name__)
    if not Path(args.data).exists():
        logger.error("Sequence data file does not exist")
        sys.exit(1)
    if args.logfile is None:
        args.logfile = Path(os.devnull)
    elif not Path(args.logfile.parent).exists():
        Path(args.logfile.parent).mkdir(parents=True)
    logging.basicConfig(format='%(asctime)s | %(levelname)s: %(message)s',
                        handlers=[
                            logging.FileHandler(args.logfile),
                            logging.StreamHandler(sys.stdout)
                            ],
                        level=logging.NOTSET)
    config = ConfigParser.get_default_config()
    args.synteny_struc = SyntenyParser.reformatSyntenyStructure(args.synteny_struc)
    if not SyntenyParser.isValidStructure(args.synteny_struc):
        logger.error(
            (
                f"Invalid synteny structure format: {args.synteny_struc}. "
                "Execute pynteny search --help to see the right format."
                )
        )
        sys.exit(1)
    if args.hmm_dir is None:
        if not config.get_field("data_downloaded"):
            logger.warning("HMM database not found. Downloading PGAP database from NCBI")
            down_args = CommandArgs(unpack=True, outdir=None, logfile=None)
            download_hmms(down_args)
        else:
            args.hmm_dir = Path(config.get_field("PGAP_database"))
    if args.gene_ids:
        if args.hmm_meta is None:
            if not config.get_field("data_downloaded"):
                logger.error("Please download hmm database first or provide path to hmm metadata file.")
                sys.exit(1)
            else:
                args.hmm_meta = Path(config.get_field("PGAP_meta_file"))
        logger.info("Finding matching HMMs for gene symbols")
        gene_synteny_struc, gene_to_hmm_group = SyntenyParser.parseGenesInSyntenyStructure(
            synteny_structure=args.synteny_struc,
            hmm_meta=args.hmm_meta
        )
        args.synteny_struc = gene_synteny_struc
        logger.info(f"Found the following HMMs in database for given structure:\n{gene_synteny_struc}")
    
    temp_hmm_dir = Path(args.hmm_dir.parent) / "temp_hmm_dir"
    if isTarFile(args.hmm_dir):
        PGAP.extractPGAPtoDirectory(args.hmm_dir, temp_hmm_dir)
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
        args.outdir = Path(args.data.parent)
    if not args.outdir.exists():
        args.outdir.mkdir(parents=True, exist_ok=True)
    if args.hmmsearch_args is None:
        hmmsearch_args = ",".join(["None" for _ in input_hmms])
    hmmsearch_args = list(map(lambda x: x.strip(), hmmsearch_args.split(",")))
    hmmsearch_args = list(map(lambda x: None if x == 'None' else x, hmmsearch_args))
    hmmer_output_dir = os.path.join(args.outdir, 'hmmer_outputs/')
    synteny_table = args.outdir / f"{args.prefix}synteny_matched.tsv"
    
    logger.info('Searching database by synteny structure')
    synteny_hits = filterFASTAbySyntenyStructure(
        synteny_structure=args.synteny_struc,
        unordered=args.unordered,
        input_fasta=args.data,
        input_hmms=input_hmms,
        hmm_meta=args.hmm_meta,
        hmmer_output_dir=hmmer_output_dir,
        reuse_hmmer_results=True,
        method='hmmsearch',
        processes=args.processes,
        additional_args=hmmsearch_args
    )
    synteny_hits.writeToTSV(synteny_table)
    logger.info("Writing matching sequences to FASTA files")
    synteny_hits.writeHitSequencesToFASTAfiles(
        sequence_database=args.data,
        output_dir=args.outdir,
        output_prefix=args.prefix
    )
    if temp_hmm_dir.exists():
        shutil.rmtree(temp_hmm_dir)
    logger.info('Finished!')
    return synteny_hits

def build_database(args) -> None:
    """
    Build annotated peptide database from input assembly
    or GenBank data.
    """
    if args.logfile is None:
        args.logfile = Path(os.devnull)
    elif not Path(args.logfile.parent).exists():
        Path(args.logfile.parent).mkdir(parents=True)
    logging.basicConfig(format='%(asctime)s | %(levelname)s: %(message)s',
                        handlers=[
                            logging.FileHandler(args.logfile),
                            logging.StreamHandler(sys.stdout)
                            ],
                        level=logging.NOTSET)
    logger = logging.getLogger(__name__)
    if args.processes is None:
        args.processes = os.cpu_count() - 1
    logger.info("Building annotated peptide database")
    database = Database(args.data)
    database.build(
        output_file=args.outfile,
    )
    logger.info("Database built successfully!")

def parse_gene_ids(args) -> str:
    """
    Convert gene symbols to hmm names.
    """
    if args.logfile is None:
        args.logfile = Path(os.devnull)
    elif not Path(args.logfile.parent).exists():
        Path(args.logfile.parent).mkdir(parents=True)
    logging.basicConfig(format='%(asctime)s | %(levelname)s: %(message)s',
                        handlers=[
                            logging.FileHandler(args.logfile),
                            logging.StreamHandler(sys.stdout)
                            ],
                        level=logging.NOTSET)
    logger = logging.getLogger(__name__)
    config = ConfigParser.get_default_config()
    if args.hmm_meta is None:
        if not config.get_field("data_downloaded"):
            logger.error("Please download hmm database meta file or provide path to existing one first.")
            sys.exit(1)
        else:
            args.hmm_meta = Path(config.get_field("PGAP_meta_file"))
    gene_synteny_struc, gene_to_hmm_group = SyntenyParser.parseGenesInSyntenyStructure(
        synteny_structure=args.synteny_struc,
        hmm_meta=args.hmm_meta
        )
    logger.info(f'Translated \n "{args.synteny_struc}" \n to \n "{gene_synteny_struc}" \n according to provided HMM database metadata')
    return gene_synteny_struc

def download_hmms(args) -> None:
    """
    Download HMM (PGAP) database from NCBI.
    """
    if args.logfile is None:
        args.logfile = Path(os.devnull)
    elif not Path(args.logfile.parent).exists():
        Path(args.logfile.parent).mkdir(parents=True)

    logging.basicConfig(format='%(asctime)s | %(levelname)s: %(message)s',
                        handlers=[
                            logging.FileHandler(args.logfile),
                            logging.StreamHandler(sys.stdout)
                            ],
                        level=logging.NOTSET)
    logger = logging.getLogger(__name__)
    module_dir = Path(__file__).parent
    config = ConfigParser.get_default_config()
    if config.get_field("data_downloaded"):
        logger.info("PGAP database already downloaded. Skipping download")
        sys.exit(1)
    if args.outdir is None:
        download_dir = Path(module_dir.parent) / "data"
    else:
        download_dir = Path(args.outdir).absolute()
    if not download_dir.exists():
        os.makedirs(download_dir, exist_ok=True)

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
        logger.info("Database dowloaded successfully\n")
        config.update_config("data_downloaded", True)
        config.update_config("PGAP_database", PGAP_file.as_posix())
        config.update_config("PGAP_meta_file",  meta_file.as_posix())
    except Exception as e:
        logger.exception("Failed to download PGAP database. Please check your internet connection.")
        sys.exit(1)
    logger.info("Removing missing entries from PGAP metadata file")
    PGAP(meta_file).removeMissingHMMsFromMetadata(
        PGAP_file,
        meta_file
    )
    if args.unpack:
        logger.info("Unpacking PGAP database")
        unpacked_PGAP_dir = download_dir / "hmm_PGAP"
        PGAP.extractPGAPtoDirectory(
            PGAP_file,
            output_dir=unpacked_PGAP_dir
        )
        os.remove(PGAP_file)
        config.update_config("PGAP_database", unpacked_PGAP_dir.as_posix())
        logger.info("PGAP database unpacked successfully")

def run_app() -> None:
    """
    Run Pynteny app through streamlit
    """
    config = ConfigParser.get_default_config()
    config_path = config.get_config_path()
    logfile = str(Path(config_path.parent) / "streamlit.log")
    config.update_config("streamlit_log", logfile)
    app_path = Path(Path(__file__).parent) / "app" / "main_page.py"
    log_str = f"--logger.level=info 2> {logfile}"
    cmd_str = (
        f"streamlit run {app_path} --browser.gatherUsageStats False --server.fileWatcherType none "
        f"{log_str}"
        )
    terminalExecute(cmd_str)

def run_tests() -> None:
    """
    Run unit and integration tests
    """
    tests_path = Path(Path(Path(__file__).parent).parent) / "tests"
    cmd_str = f"python -m unittest discover {tests_path}"
    terminalExecute(cmd_str, work_dir=tests_path)

def get_citation(args, silent: bool = False) -> str:
    """
    Get Pynteny citation string
    """
    citation = (
        "Semidán Robaina Estévez (2022). Pynteny: synteny-aware hmm searches made easy"
        f"(Version {args.version}). Zenodo. https://doi.org/10.5281/zenodo.7048685"
    )
    if not silent:
        print("If you use this software, please cite it as below: ")
        print(citation)
    return citation
