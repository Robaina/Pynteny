#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to parse Hmmer output and PGAP (HMM) database
"""

from __future__ import annotations

import logging
import os
import sys
from collections import defaultdict
from pathlib import Path
import tempfile

import pandas as pd
from Bio import SearchIO

import pynteny.wrappers as wrappers
from pynteny.utils import (
    extract_tar_file,
    flatten_directory,
    is_tar_file,
    list_tar_dir,
    download_file,
    extract_gz_file,
    split_hmms,
)

logger = logging.getLogger(__name__)


class HMMER:
    """Run Hmmer on multiple hmms and parse output"""

    def __init__(
        self,
        input_hmms: list[Path],
        hmm_output_dir: Path,
        input_data: Path,
        additional_args: list[str] = None,
        processes: int = None,
    ):
        """Initialize class HMMER

        Args:
            input_hmms (list[Path]): list of paths to input HMM files.
            hmm_output_dir (Path): path to output directory to HMMER output files.
            input_data (Path): path to input fasta file with sequence database.
            additional_args (list[str]): additional arguments to hmmsearch or hmmscan. Each
                element in the list is a string with additional arguments for each
                input hmm (arranged in the same order), an element can also take a
                value of None to avoid passing additional arguments for a specific
                input hmm. A single string may also be passed, in which case the
                same additional argument is passed to hmmsearch for all input hmms.
                Defaults to None.
            processes (int, optional): maximum number of threads to be employed.
                Defaults to all minus one.
        """
        if additional_args is None:
            additional_args = [None for _ in range(len(input_hmms))]
        self._hmmer_output_dir = Path(hmm_output_dir)
        self._input_hmms = [Path(hmm) for hmm in input_hmms]
        self._input_fasta = Path(input_data)
        self._additional_args = additional_args
        self._processes = processes

    @property
    def hmm_names(self) -> list[str]:
        """Get file names of input HMMs

        Returns:
            list[str]: list of file names.
        """
        return [hmm_path.stem for hmm_path in self._input_hmms]

    @staticmethod
    def parse_HMM_search_output(hmmer_output: Path) -> pd.DataFrame:
        """Parse hmmsearch or hmmscan summary table output file.

        Args:
            hmmer_output (str): path to HMMER output file.

        Returns:
            pd.DataFrame: a dataframe containing parsed HMMER output.
        """
        hmmer_output = Path(hmmer_output)
        attribs = ["id", "bias", "bitscore", "description"]
        hits = defaultdict(list)
        with open(hmmer_output, encoding="UTF-8") as handle:
            for queryresult in SearchIO.parse(handle, "hmmer3-tab"):
                for hit in queryresult.hits:
                    for attrib in attribs:
                        hits[attrib].append(getattr(hit, attrib))
        return pd.DataFrame.from_dict(hits)

    def get_HMMER_tables(
        self, reuse_hmmer_results: bool = True, method: str = "hmmsearch"
    ) -> dict[pd.DataFrame]:
        """Run hmmer for given hmm list

        Args:
            reuse_hmmer_results (bool, optional): if True then HMMER3 won't be run
                again for HMMs already searched in the same output directory. Defaults to True.
            method (str, optional): select between 'hmmsearch' or 'hmmscan'. Defaults to 'hmmsearch'.

        Returns:
            dict[pd.DataFrame]: dict of HMMER hits as pandas dataframes.
        """
        hmm_hits = {}
        for hmm_model, add_args in zip(self._input_hmms, self._additional_args):
            hmm_name = hmm_model.stem
            hmmer_output = Path(self._hmmer_output_dir) / f"hmmer_output_{hmm_name}.txt"

            if not (reuse_hmmer_results and os.path.isfile(hmmer_output)):
                wrappers.run_HMM_search(
                    hmm_model=hmm_model,
                    input_fasta=self._input_fasta,
                    output_file=hmmer_output,
                    method=method,
                    processes=self._processes,
                    additional_args=add_args,
                )
            elif reuse_hmmer_results and os.path.isfile(hmmer_output):
                logger.info(f"Reusing Hmmer results for HMM: {hmm_name}")
            hmm_hits[hmm_name] = HMMER.parse_HMM_search_output(hmmer_output)
        return hmm_hits


class HMMDatabase:
    """Base class for HMM databases"""

    def __init__(
        self,
        database_directory: Path,
        metadata_file: Path,
        metadata_columns: list[str] = None,
    ) -> None:
        """Initialize class HMMDatabase

        Args:
            database_directory (Path): path to directory containing HMM database
                                       with individual files for each HMM.
            metadata_file (Path): path to metadata file.
            metadata_columns (list[str], optional): list of metadata columns to be
                                                    used. Defaults to None (all columns).
        """
        self._data_dir = Path(database_directory)
        self._metadata_file = Path(metadata_file)
        self._meta = pd.read_csv(
            self._metadata_file, sep="\t", usecols=metadata_columns
        )

    @property
    def meta(self) -> pd.DataFrame:
        """Return PFAM metadata

        Returns:
            pd.DataFrame: PFAM metadata
        """
        return self._meta

    @property
    def meta_file(self) -> Path:
        """Return path to  metadata file

        Returns:
            Path: metadata file
        """
        return self._metadata_file

    @property
    def data_dir(self) -> Path:
        """Return data directory

        Returns:
            Path: data directory
        """
        return self._data_dir

    def get_HMM_names_by_gene_symbol(self, gene_symbol: str) -> list[str]:
        """Try to retrieve HMM by its gene symbol, more
           than one HMM may map to a single gene symbol

        Args:
            gene_symbol (str): gene symbol to be searched for HMM.

        Returns:
            list[str]: list of HMM names matching gene symbol.
        """
        meta = self._meta  # .dropna(subset=["gene_symbol", "label"], axis=0)
        return meta[((meta.gene_symbol == gene_symbol) | (meta.label == gene_symbol))][
            "accession"
        ].values.tolist()

    def get_HMM_group_for_gene_symbol(self, gene_symbol: str) -> str:
        """Get HMMs corresponding to gene symbol in PGAP metadata.
           If more than one HMM matching gene symbol, return a HMM group

        Args:
            gene_symbol (str): gene symbol to be searched for HMM.

        Returns:
            str: string of HMM names (group separated by |)
        """
        hmms = self.get_HMM_names_by_gene_symbol(gene_symbol)
        if not hmms:
            logger.error(f"No HMM found for gene {gene_symbol}")
            sys.exit(1)
        if len(hmms) == 1:
            return hmms.pop()
        else:
            return "|".join(hmms)

    def get_HMM_gene_ID(self, hmm_name: str) -> list[str]:
        """Get gene symbols matching given hmm.

        Args:
            hmm_name (str): query HMM name.

        Returns:
            list[str]: list of gene symbols matching given HMM.
        """
        meta = self._meta.dropna(subset=["accession"], axis=0)
        return meta[meta["accession"] == hmm_name]["gene_symbol"].values.tolist()

    def get_meta_info_for_HMM(self, hmm_name: str) -> dict:
        """Get meta info for given hmm.

        Args:
            hmm_name (str): query HMM name.

        Returns:
            dict: metadata of provided HMM.
        """
        meta = self._meta.dropna(subset=["accession"], axis=0).applymap(
            lambda x: x if not pd.isna(x) else ""
        )
        metadata = {
            k: list(v.values())[0] if list(v.values())[0] else "undef"
            for k, v in meta[meta["accession"] == hmm_name].to_dict().items()
        }
        if not metadata:
            logger.warning(f"No metadata for HMM: {hmm_name}")
        return metadata


class PGAP(HMMDatabase):
    """Tools to parse PGAP hmm database metadata"""

    def __init__(self, *args, **kwargs):
        """Initialize class PGAP"""
        super().__init__(*args, **kwargs)
        self._meta = self._meta.rename(columns={"#ncbi_accession": "accession"})
        # self._meta = self.remove_missing_HMMs_from_metadata(meta_outfile=None)

    @staticmethod
    def remove_missing_HMMs_from_metadata(self, meta_outfile: Path = None) -> None:
        """Remove HMMs from metadata that are not in HMM directory

        Args:
            hmm_dir (Path): path to directory containing PGAP database.
            meta_outfile (Path, optional): path to output file. Defaults to None.
        """
        logger.info("Removing missing HMMs from PGAP metadata")
        hmm_dir = self._database_dir
        if meta_outfile is None:
            meta_outfile = (
                self._metadata_file.parent
                / f"{self._metadata_file.stem}_preprocessed.tsv"
            )
        else:
            meta_outfile = Path(meta_outfile)
        if is_tar_file(hmm_dir):
            hmm_file_names = [
                Path(hmm_file).stem.strip() for hmm_file in list_tar_dir(hmm_dir)
            ]
        else:
            hmm_file_names = [hmm_file.stem.strip() for hmm_file in hmm_dir.iterdir()]
        not_found = set()
        for i, row in self._meta.iterrows():
            if row["accession"].strip() not in hmm_file_names:
                not_found.add(i)
        self._meta = self._meta.drop(not_found)
        self._metadata_file = meta_outfile
        self._meta.to_csv(meta_outfile, sep="\t", index=False)


class PFAM(HMMDatabase):
    """Tools to preprocess the PFAM-A hmm database"""

    @classmethod
    def from_gz_file(
        cls, hmm_gz_file: Path, hmm_outdir: Path = None, meta_outfile: Path = None
    ) -> PFAM:
        """Initialize PFAM class from hmm file

        Args:
            hmm_gz_file (Path): path to hmm gz file.
            hmm_outdir (Path): path to output directory to store individual HMM files.
            meta_outfile (Path, optional): path to metadata output file. Defaults to None.
        """
        hmm_gz_file = Path(hmm_gz_file)
        if hmm_outdir is not None:
            hmm_outdir = Path(hmm_outdir)
        else:
            hmm_outdir = hmm_gz_file.parent / "pfam_hmms"
        if meta_outfile is not None:
            meta_outfile = Path(meta_outfile)
        else:
            meta_outfile = hmm_outdir / f"{hmm_gz_file.stem}_meta.tsv"
        with tempfile.NamedTemporaryFile() as temp:
            extract_gz_file(hmm_gz_file, temp.name)
            split_hmms(temp.name, hmm_outdir)
            cls.construct_meta_file(hmm_outdir, meta_outfile)
        return cls(hmm_outdir, meta_outfile)

    def construct_meta_file(self, meta_outfile: Path = None) -> None:
        """Construct metadata file from individual HMM files.

        Args:
            meta_outfile (Path): path to metadata file.
        """
        logger.info("Constructing metadata file for PFAM-A database")
        if meta_outfile is None:
            meta_outfile = self._database_dir / "PFAM_meta.tsv"
        else:
            meta_outfile = Path(meta_outfile)
        hmm_meta_lines = ["accession\tgene_symbol\tdescription\tlength\n"]
        for hmm_file in self._database_dir.glob("*.hmm"):
            with open(hmm_file, "r") as f:
                hmm_text = f.read()

            acc_code = [
                entry for entry in hmm_text.split("\n") if entry.startswith("ACC")
            ]
            name = [entry for entry in hmm_text.split("\n") if entry.startswith("NAME")]
            description = [
                entry for entry in hmm_text.split("\n") if entry.startswith("DESC")
            ]
            length = [
                entry for entry in hmm_text.split("\n") if entry.startswith("LENG")
            ]
            name = name[0].split()[1] if name else "Unspecified"
            description = description[0].split()[1] if description else "Unspecified"
            length = length[0].split()[1] if length else "Unspecified"

            if acc_code:
                acc_name = acc_code[0].split()[1]
                hmm_meta_lines.append(f"{acc_name}\t{name}\t{description}\t{length}\n")
        with open(meta_outfile, "w+") as f:
            f.writelines(hmm_meta_lines)
        self._metadata_file = meta_outfile
        self._meta = pd.read_csv(meta_outfile, sep="\t")


def download_pgap(download_dir: Path, unpack: bool = False) -> tuple[Path, Path]:
    """Download PGAP database

    Args:
        download_dir (Path): path to output directory.
        unpack (bool, optional): if True then PGAP database will be extracted
    """
    if download_dir.exists():
        logger.warning(
            f"{download_dir} already exists. Downloader may overwrite files."
        )

    data_url = "https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz"
    meta_url = "https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv"
    PGAP_file = download_dir / "hmm_PGAP.HMM.tgz"
    meta_file = download_dir / "hmm_PGAP.tsv"
    download_file(data_url, PGAP_file)
    download_file(meta_url, meta_file)
    if unpack:
        destination_path = download_dir / "pgap_hmms"
        extract_pgap_to_directory(PGAP_file, destination_dir=destination_path)
        return destination_path, meta_file
    else:
        return PGAP_file, meta_file


def download_pfam(download_dir: Path, unpack: bool = False) -> Path:
    """Download PFAM database

    Args:
        unpack (bool, optional): if True then PFAM database will be extracted
    """
    if download_dir.exists():
        logger.warning(
            f"{download_dir} already exists. Downloader may overwrite files."
        )
    PFAM_file = download_dir / "Pfam-A.gz"
    logger.info("Downloading PFAM-A hmm database")
    url = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
    download_file(url, PFAM_file)
    if unpack:
        destination_path = download_dir / "pfam_hmms"
        extract_pfam_to_directory(PFAM_file, destination_dir=destination_path)
        return destination_path
    else:
        return PFAM_file


def extract_pgap_to_directory(pgap_tar: Path, destination_dir: Path) -> None:
    """Extract PGAP hmm database (tar.gz) to downlaod directory

    Args:
        pgap_tar (Path): path to compressed PGAP database.
    """
    pgap_tar = Path(pgap_tar)
    if not is_tar_file(pgap_tar):
        logger.warning(f"{pgap_tar} is not a tar file. Skipping extraction")
        return
    logger.info("Extracting hmm files to target directory")
    extract_tar_file(pgap_tar, destination_dir)
    flatten_directory(destination_dir)
    os.remove(pgap_tar)
    logger.info("PGAP database unpacked successfully")


def extract_pfam_to_directory(pfam_gz: Path, destination_dir: Path) -> None:
    """Extract PFAM hmm database (gz) to downlaod directory

    Args:
        pfam_gz (Path): path to compressed PFAM database.
    """
    pfam_gz = Path(pfam_gz)
    if not pfam_gz.is_file():
        logger.warning(f"{pfam_gz} is not a file. Skipping extraction")
        return
    logger.info("Extracting hmm files to target directory")
    extract_gz_file(pfam_gz, destination_dir)
    flatten_directory(destination_dir)
    os.remove(pfam_gz)
    logger.info("PGAP database unpacked successfully")
