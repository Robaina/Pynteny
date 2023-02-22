#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to parse Hmmer output and PGAP (HMM) database
"""

from __future__ import annotations

import logging
import os
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SearchIO

import pynteny.wrappers as wrappers
from pynteny.utils import extract_tar_file, flatten_directory, is_tar_file, list_tar_dir

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
                    n_processes=self._processes,
                    additional_args=add_args,
                )
            elif reuse_hmmer_results and os.path.isfile(hmmer_output):
                logger.info(f"Reusing Hmmer results for HMM: {hmm_name}")
            hmm_hits[hmm_name] = HMMER.parse_HMM_search_output(hmmer_output)
        return hmm_hits


class PGAP:
    """Tools to parse PGAP hmm database metadata"""

    def __init__(self, meta_file: Path):
        """Initialize class PGAP

        Args:
            meta_file (Path): path to PGAP's metadata file.
        """
        meta_file = Path(meta_file)
        meta = pd.read_csv(
            str(meta_file),
            sep="\t",
            usecols=[
                "#ncbi_accession",
                "gene_symbol",
                "label",
                "product_name",
                "ec_numbers",
            ],
        )
        # meta = meta[
        #     ["#ncbi_accession", "gene_symbol", "label", "product_name", "ec_numbers"]
        # ]
        self._meta = meta
        self._meta_file = meta_file

    @staticmethod
    def extract_PGAP_to_directory(pgap_tar: Path, output_dir: Path) -> None:
        """Extract PGAP hmm database (tar.gz) to given directory

        Args:
            pgap_tar (Path): path to compressed PGAP database.
            output_dir (Path): path to output directory.
        """
        pgap_tar = Path(pgap_tar)
        if not is_tar_file(pgap_tar):
            logger.warning(f"{pgap_tar} is not a tar file. Skipping extraction")
            sys.exit(1)
        logger.info("Extracting hmm files to target directory")
        output_dir = Path(output_dir)
        if output_dir.exists():
            shutil.rmtree(output_dir)
        extract_tar_file(tar_file=pgap_tar, dest_dir=output_dir)
        flatten_directory(output_dir)

    def remove_missing_HMMs_from_metadata(
        self, hmm_dir: Path, output_file: Path = None
    ) -> None:
        """Remove HMMs from metadata that are not in HMM directory

        Args:
            hmm_dir (Path): path to directory containing PGAP database.
            output_file (Path, optional): path to output file. Defaults to None.
        """
        hmm_dir = Path(hmm_dir)
        if output_file is None:
            output_file = self._meta.parent / f"{self._meta.stem}_missing_hmms.tsv"
        else:
            output_file = Path(output_file)
        if is_tar_file(hmm_dir):
            hmm_file_names = [
                Path(hmm_file).stem.strip() for hmm_file in list_tar_dir(hmm_dir)
            ]
        else:
            hmm_file_names = [hmm_file.stem.strip() for hmm_file in hmm_dir.iterdir()]
        not_found = set()
        for i, row in self._meta.iterrows():
            if row["#ncbi_accession"].strip() not in hmm_file_names:
                not_found.add(i)
        self._meta = self._meta.drop(not_found)
        self._meta.to_csv(output_file, sep="\t", index=False)

    def get_HMM_names_by_gene_symbol(self, gene_symbol: str) -> list[str]:
        """Try to retrieve HMM by its gene symbol, more
           than one HMM may map to a single gene symbol

        Args:
            gene_symbol (str): gene symbol to be searched for HMM.

        Returns:
            list[str]: list of HMM names matching gene symbol.
        """
        meta = self._meta  # .dropna(subset=["gene_symbol", "label"], axis=0)
        return meta[
            (
                (meta.gene_symbol == gene_symbol)
                |
                # (meta.label.str.contains(gene_id))
                (meta.label == gene_symbol)
            )
        ]["#ncbi_accession"].values.tolist()

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
        meta = self._meta.dropna(subset=["#ncbi_accession"], axis=0)
        return meta[meta["#ncbi_accession"] == hmm_name]["gene_symbol"].values.tolist()

    def get_meta_info_for_HMM(self, hmm_name: str) -> dict:
        """Get meta info for given hmm.

        Args:
            hmm_name (str): query HMM name.

        Returns:
            dict: metadata of provided HMM.
        """
        meta = self._meta.dropna(subset=["#ncbi_accession"], axis=0).applymap(
            lambda x: x if not pd.isna(x) else ""
        )
        metadata = {
            k: list(v.values())[0] if list(v.values())[0] else "undef"
            for k, v in meta[meta["#ncbi_accession"] == hmm_name].to_dict().items()
        }
        if not metadata:
            logger.warning(f"No metadata for HMM: {hmm_name}")
        return metadata
