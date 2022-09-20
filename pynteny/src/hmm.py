#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to parse Hmmer output and PGAP (HMM) database
"""

from __future__ import annotations
import os
import sys
import shutil
import logging
from pathlib import Path
from collections import defaultdict

import pandas as pd
from Bio import SearchIO

import pynteny.src.wrappers as wrappers
from pynteny.src.utils import isTarFile, extractTarFile, flattenDirectory, listTarDir


logger = logging.getLogger(__name__)


class HMMER:
    def __init__(self, input_hmms: list[Path],
                 hmm_output_dir: Path,
                 input_data: Path,
                 additional_args: list[str],
                 processes: int = None) -> None:
        """
        Run Hmmer on multiple hmms and parse output
        """
        self._hmmer_output_dir = hmm_output_dir
        self._input_hmms = input_hmms
        self._input_fasta = input_data
        self._additional_args = additional_args
        self._processes = processes

    @property
    def hmm_names(self) -> list[str]:
        return [hmm_path.stem for hmm_path in self._input_hmms]

    @staticmethod
    def parseHMMsearchOutput(hmmer_output: str) -> pd.DataFrame:
        """
        Parse hmmsearch or hmmscan summary table output file
        """
        attribs = ['id', 'bias', 'bitscore', 'description']
        hits = defaultdict(list)
        with open(hmmer_output) as handle:
            for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
                for hit in queryresult.hits:
                    for attrib in attribs:
                        hits[attrib].append(getattr(hit, attrib))
        return pd.DataFrame.from_dict(hits)
    
    def getHMMERtables(self,
                       reuse_hmmer_results: bool = True,
                       method: str = None) -> dict[pd.DataFrame]:
        """
        Run hmmer for given hmm list
        """
        hmm_hits = {}
        for hmm_model, add_args in zip(self._input_hmms, self._additional_args):
            hmm_name = hmm_model.stem
            hmmer_output = Path(os.path.join(self._hmmer_output_dir, f'hmmer_output_{hmm_name}.txt'))
            if not (reuse_hmmer_results and os.path.isfile(hmmer_output)):
                wrappers.runHMMsearch(
                    hmm_model=hmm_model,
                    input_fasta=self._input_fasta,
                    output_file=hmmer_output,
                    method=method,
                    n_processes=self._processes,
                    additional_args=add_args
                    )
            elif reuse_hmmer_results and os.path.isfile(hmmer_output):
                logger.info(f"Reusing Hmmer results for HMM: {hmm_name}")
            hmm_hits[hmm_name] = HMMER.parseHMMsearchOutput(hmmer_output)
        return hmm_hits

    
class PGAP:
    def __init__(self, meta_file: Path) -> None:
        """
        Tools to parse PGAP hmm database metadata
        """
        meta = pd.read_csv(str(meta_file), sep="\t")
        meta = meta[["#ncbi_accession", "gene_symbol", "label", "product_name", "ec_numbers"]]
        self._meta = meta
        self._meta_file = meta_file

    @staticmethod
    def extractPGAPtoDirectory(pgap_tar: Path, output_dir: Path) -> None:
        """
        Extract PGAP hmm database (tar.gz) to given directory
        """
        if not isTarFile(pgap_tar):
            logger.warning(f"{pgap_tar} is not a tar file. Skipping extraction")
            sys.exit(1)
        logger.info("Extracting hmm files to target directory")
        if output_dir.exists():
            shutil.rmtree(output_dir)
        extractTarFile(
            tar_file=pgap_tar,
            dest_dir=output_dir
        )
        flattenDirectory(
            output_dir
        )

    def removeMissingHMMsFromMetadata(self, hmm_dir: Path,
                                      outfile: Path = None) -> None:
        """
        Remove HMMs from metadata that are not in HMM directory
        """
        if outfile is None:
            outfile = self._meta.parent / f"{self._meta.stem}_missing_hmms.tsv"
        if isTarFile(hmm_dir):
            hmm_file_names = [
                Path(hmm_file).stem.strip()
                for hmm_file in listTarDir(hmm_dir)
                ]
        else:
            hmm_file_names = [
                hmm_file.stem.strip()
                for hmm_file in hmm_dir.iterdir()
                ]
        not_found = set() 
        for i, row in self._meta.iterrows():
            if row["#ncbi_accession"].strip() not in hmm_file_names:
                not_found.add(i)
        self._meta = self._meta.drop(not_found)
        self._meta.to_csv(outfile, sep="\t", index=False)

    def getHMMnamesByGeneSymbol(self, gene_id: str) -> list[str]:
        """
        Try to retrieve HMM by its gene symbol, more
        than one HMM may map to a single gene symbol
        """
        meta = self._meta#.dropna(subset=["gene_symbol", "label"], axis=0)
        try:
            return meta[
                (
                    (meta.gene_symbol == gene_id) |
                    # (meta.label.str.contains(gene_id))
                    (meta.label == gene_id)
                    )
                ]["#ncbi_accession"].values.tolist()
        except:
            return list()

    def getHMMgroupForGeneSymbol(self, gene_symbol: str) -> str:
        """
        Get HMMs corresponding to gene symbol in PGAP metadata.
        If more than one HMM matching gene symbol, return a HMM group
        """
        hmms = self.getHMMnamesByGeneSymbol(gene_symbol)
        if not hmms:
            logger.error(f"No HMM found for gene {gene_symbol}")
            sys.exit(1)
        if len(hmms) == 1:
            return hmms.pop()
        else:
            return "|".join(hmms)
    
    def getHMMgeneID(self, hmm_name: str) -> list[str]: 
        """
        Get gene symbol of given hmm
        """
        meta = self._meta.dropna(subset=["#ncbi_accession"], axis=0)
        try:
            return meta[
                meta["#ncbi_accession"] == hmm_name
            ]["gene_symbol"].values.tolist()
        except:
            return None

    def getMetaInfoForHMM(self, hmm_name: str) -> dict:
        """
        Get meta info for given hmm
        """
        meta = self._meta.dropna(subset=["#ncbi_accession"], axis=0).applymap(
            lambda x: x if not pd.isna(x) else ""
            )
        try: 
            return {
                k: list(v.values())[0] if list(v.values())[0] else "undef"
                for k, v in meta[meta["#ncbi_accession"] == hmm_name].to_dict().items()
                }
        except:
            logger.warning(f"No metadata for HMM: {hmm_name}")
            return dict()