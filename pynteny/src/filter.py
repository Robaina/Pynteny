#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to filter HMM hits by synteny structure
"""

from __future__ import annotations
import time
import os
import sys
import logging
from pathlib import Path

import pandas as pd

from pynteny.src.preprocessing import FASTA
from pynteny.src.hmm import HMMER, PGAP
from pynteny.src.parser import SyntenyParser, LabelParser

logger = logging.getLogger(__name__)


class SyntenyPatternFilters():
    def __init__(self, synteny_structure: str, unordered: bool = False) -> None:
        """
        Methods to filter hmm hits in the  same contig by synteny structure
        or collinearity. To be applied in pandas.Series.rolling.apply method.

        @param: synteny_structure: string coontaining the synteny structure
                to be searched for
        @param: unordered, boolean, whether the HMMs should be arranged in the
                exact same order displayed in the synteny_structure or in 
                any order If ordered, the filters would filter collinear rather
                than syntenic structures.
        """
        parsed_structure = SyntenyParser.parseSyntenyStructure(synteny_structure)
        hmm_order_dict = dict(
            zip(parsed_structure["hmm_groups"], range(len(parsed_structure["hmm_groups"])))
            )
        hmm_codes = list(hmm_order_dict.values())
        self.hmm_code_order_pattern = hmm_codes

        if unordered:
            max_distance = max(parsed_structure["distances"])
            distances = [float("inf")] + [max_distance for _ in parsed_structure["distances"]]
            self.distance_order_pattern = [dist + 1 for dist in distances]
        else:
            distances = [float("inf")] + parsed_structure["distances"]
            self.distance_order_pattern = [dist + 1 for dist in distances]

        self.strand_order_pattern = list(map(
            lambda strand : -1 if strand == "neg" else (1 if strand ==  "pos" else 0),
            parsed_structure["strands"]
            ))
        self._unordered = unordered

    def contains_hmm_pattern(self, data: pd.Series) -> int:
        if self._unordered:
            return 1 if set(data.values) == set(self.hmm_code_order_pattern) else 0
        else:
            return 1 if data.values.tolist() == self.hmm_code_order_pattern else 0

    def contains_distance_pattern(self, data: pd.Series) -> int:
        return 1 if all(
            [
                (data_dist <= pattern_dist) and (data_dist > 0)
                for data_dist, pattern_dist in zip(data.values.tolist(), self.distance_order_pattern)
            ]
            ) else 0

    def contains_strand_pattern(self, data: pd.Series) -> int:
        strand_comparisons = []
        for data_strand, pattern_strand in zip(data.values, self.strand_order_pattern):
            if pattern_strand != 0:
                strand_comparisons.append(data_strand == pattern_strand)
            else:
                strand_comparisons.append(True)
        return 1 if all(strand_comparisons) == True else 0


class SyntenyHMMfilter():
    """
    Tools to search for synteny structures among sets of hmm models
    """
    def __init__(self, hmm_hits: dict, synteny_structure: str, unordered: bool = True) -> None:
        """
        Search for contigs that satisfy the given gene synteny structure

        @param: hmm_hits, a dict of pandas DataFrames, as output by
                parseHMMsearchOutput with keys corresponding to hmm names

        @param: synteny_structure, a str describing the desired synteny structure,
                structured as follows:

                'hmm_a N_ab hmm_b'

                where N_ab corresponds to the maximum number of genes separating 
                gene found by hmm_a and gene found by hmm_b, and hmm_ corresponds 
                to the name of the hmm as provided in the keys of hmm_hits.
                More than two hmms can be concatenated.

        @param: uordered, boolean, whether the HMMs should be arranged in the
                exact same order displayed in the synteny_structure or in 
                any order If ordered, the filters would filter collinear rather
                than syntenic structures.
        """
        self._unordered = unordered
        self._hmm_hits = hmm_hits
        self._hmms = list(hmm_hits.keys())
        self._synteny_structure = synteny_structure
        self._parsed_structure = SyntenyParser.parseSyntenyStructure(self._synteny_structure)
        if self._unordered:
            self._parsed_structure["strands"] = ["" for _ in self._parsed_structure["strands"]]
            logger.warning("Disregarding strand location as unordered structured selected")
        self._hmm_order_dict = dict(
            zip(self._parsed_structure["hmm_groups"], range(len(self._parsed_structure["hmm_groups"])))
            )
        self._n_hmm_groups = len(self._hmm_order_dict)
        self._unordered = unordered

    def _assignCodeToHMM(self, hmm_group_name: str) -> int:
        hmm_group_name = str(hmm_group_name)
        code = [
            code for hmm_group, code in self._hmm_order_dict.items()
            if set(hmm_group_name.split("|")).issubset(set(hmm_group.split("|")))
            ]
        if len(code) > 1:
            logger.error(
            f"HMM: {hmm_group_name} found in more than one hmm group in synteny structure"
            )
            sys.exit(1)
        if not code:
            logger.error(
            f"No code found for HMM: {hmm_group_name}"
            )
            sys.exit(1)
        return code[0]

    def _addMetaCodesToHMMhits(self, all_hit_labels: pd.Dataframe) -> pd.Dataframe:
        """
        Add numeric codes for each hmm and strand, compute distance between genes
        """
        all_hit_labels["gene_pos_diff"] = all_hit_labels.gene_pos.diff()
        all_hit_labels.loc[0, "gene_pos_diff"] = 1 # required for rolling (skips first nan)
        all_hit_labels["hmm_code"] = all_hit_labels.hmm.apply(self._assignCodeToHMM)
        all_hit_labels["strand"] = all_hit_labels.strand.apply(
            lambda strand : -1 if strand == "neg" else (1 if strand ==  "pos" else 0)
            )
        return all_hit_labels

    def mergeHitsByHMMgroup(self, hits: pd.DataFrame):
        groups = hits.groupby(["full_label"]).groups
        for group_idxs in groups.values():
            hmm_names = set(hits.loc[group_idxs, "hmm"].values)
            candidate_hmm_group = [
                hmm_group_str for hmm_group_str in self._parsed_structure["hmm_groups"]
                if hmm_names.issubset(set(hmm_group_str.split("|")))
                ]
            if candidate_hmm_group:
                hmm_group = "|".join(hmm_names)
            else:
                hmm_group = "discard"
            hits.loc[group_idxs, "hmm"] = hmm_group
        hits = hits[hits.hmm != "discard"].drop_duplicates()
        hits.hmm = hits.hmm.apply(lambda x: str(x))
        return hits

    def getAllHMMhits(self) -> pd.DataFrame:
        """
        Group and preprocess all hit labels into a single dataframe
        """
        hit_labels = {}
        labelparser = LabelParser()
        for hmm, hits in self._hmm_hits.items():
            labels = hits.id.values.tolist()
            if not labels:
                logger.error(
                    f'No records found in database matching HMM: {hmm}'
                    )
                sys.exit(1)
            hit_labels[hmm] = labelparser.parse_from_list(labels)
            hit_labels[hmm]["hmm"] = hmm
        # Create single dataframe with new column corresponding to HMM and all hits
        # Remove contigs with less hits than the number of hmms in synteny structure (drop if enabling partial searching)
        all_hit_labels = pd.concat(
            hit_labels.values()
            ).groupby("contig").filter(
                lambda x : len(x) >= self._n_hmm_groups
                ).sort_values(["contig", "gene_pos"], ascending=True)
        all_hit_labels.reset_index(drop=True, inplace=True)
        all_hit_labels_merged = self.mergeHitsByHMMgroup(all_hit_labels)
        all_hit_labels_merged.reset_index(drop=True, inplace=True)
        return self._addMetaCodesToHMMhits(all_hit_labels_merged)

    def filterHitsBySyntenyStructure(self) -> dict:
        """
        Search for contigs that satisfy the given gene synteny structure
        """
        all_matched_hits = {}
        filters = SyntenyPatternFilters(self._synteny_structure, unordered=self._unordered)
        all_hit_labels = self.getAllHMMhits()
        contig_names = all_hit_labels.contig.unique()

        for contig in contig_names:
            contig_hits = all_hit_labels[all_hit_labels.contig == contig].reset_index(drop=True)
            matched_hit_labels = {hmm_group: [] for hmm_group in contig_hits.hmm.unique()}

            if len(contig_hits.hmm.unique()) >= self._n_hmm_groups:
                
                hmm_match = contig_hits.hmm_code.rolling(window=self._n_hmm_groups).apply(
                    filters.contains_hmm_pattern
                    )
                strand_match = contig_hits.strand.rolling(window=self._n_hmm_groups).apply(
                    filters.contains_strand_pattern
                    )
                if self._n_hmm_groups > 1:
                    distance_match = contig_hits.gene_pos_diff.rolling(
                        window=self._n_hmm_groups).apply(filters.contains_distance_pattern)
                    matched_rows = contig_hits[
                        (hmm_match == 1) &
                        (strand_match == 1) &
                        (distance_match == 1)
                    ]
                else:
                    matched_rows = contig_hits[
                        (hmm_match == 1) &
                        (strand_match == 1)
                    ]
                
                for i, _ in matched_rows.iterrows():
                    matched_hits = contig_hits.iloc[i - (self._n_hmm_groups - 1): i + 1, :]
                    for label, hmm in zip((matched_hits.full_label.values), matched_hits.hmm):
                        matched_hit_labels[hmm].append(label)

                all_matched_hits[contig] = matched_hit_labels
        return all_matched_hits


class SyntenyHits():
    def __init__(self, synteny_hits: pd.DataFrame) -> None:
        """
        Class to store and manipulate synteny hits by contig
        """
        self._synteny_hits = synteny_hits

    @staticmethod
    def _hitsToDataframe(hits_by_contig) -> pd.DataFrame:
        """
        Return synteny hits as a dataframe
        """
        data = []
        columns = ["contig", "gene_id", "gene_number", "locus", "strand", "full_label", "hmm"]
        for contig, matched_hits in hits_by_contig.items():
            for hmm, labels in matched_hits.items():
                for label in labels:
                    parsed_label = LabelParser.parse(label)
                    data.append(
                        [
                            parsed_label['contig'], parsed_label['gene_id'],
                            parsed_label['gene_pos'], parsed_label['locus_pos'],
                            parsed_label['strand'], parsed_label["full_label"], hmm
                            ]
                    )
        return pd.DataFrame(data, columns=columns).sort_values(["contig", "gene_number"])

    @classmethod
    def fromHitsDict(cls, hits_by_contig: dict) -> SyntenyHits:
        """
        Return SyntenyHits object from hits_by_contig dictionary
        """
        return cls(cls._hitsToDataframe(hits_by_contig))

    def getSyntenyHits(self) -> pd.DataFrame:
        """
        Return synteny hits
        """
        return self._synteny_hits

    def addHMMmetaInfoToHits(self, hmm_meta: Path) -> SyntenyHits:
        """
        Add molecular metadata to synteny hits
        """
        fields = ["gene_symbol", "label", "product", "ec_number"]
        if all([f in self._synteny_hits.columns for f in fields]):
            return self._synteny_hits
        pgap = PGAP(hmm_meta)
        self._synteny_hits[fields] = ""
        for i, row in self._synteny_hits.iterrows():
            meta_values = [
                [
                    str(v).replace("nan", "") 
                    for k, v in pgap.getMetaInfoForHMM(hmm).items()
                    if k != "#ncbi_accession"
                ]
                for hmm in row.hmm.split("|")
            ]
            self._synteny_hits.at[
                i, fields
                ] = ["|".join(v) for v in zip(*meta_values)]
        return SyntenyHits(self._synteny_hits)

    def writeToTSV(self, output_tsv: Path) -> None:
        """
        Write synteny hits to a TSV file
        """
        self._synteny_hits.to_csv(output_tsv, sep="\t", index=False)

    def writeHitSequencesToFASTAfiles(self, sequence_database: Path,
                                      output_dir: Path = None,
                                      output_prefix: str = None) -> None:
        """
        Write matching sequences to FASTA files
        """
        fasta = FASTA(sequence_database)
        hmm_groups = self._synteny_hits.hmm.unique().tolist()

        for hmm_group in hmm_groups:
            record_ids = self._synteny_hits[
                self._synteny_hits.hmm == hmm_group
                ].full_label.values.tolist()

            if ("gene_symbol" in self._synteny_hits.columns and 
                "label" in self._synteny_hits.columns):
                gene_symbol = self._synteny_hits[self._synteny_hits.hmm == hmm_group].gene_symbol.unique()[0]
                gene_label = self._synteny_hits[self._synteny_hits.hmm == hmm_group].label.unique()[0]
                gene_id = (
                    gene_symbol if not pd.isna(gene_symbol)
                    else (gene_label if not pd.isna(gene_label) else "")
                    )
            else:
                gene_id = ""
            gene_id_str = gene_id + "_" if gene_id else ""

            output_fasta = output_dir / f"{output_prefix}{gene_id_str.replace('|', '_')}{hmm_group.replace('|', '_')}_hits.fasta"
            if record_ids:
                fasta.filterByIDs(
                    record_ids=record_ids,
                    output_file=output_fasta,
                    point_to_new_file=False

                )
            else:
                logger.warning(f"No record matches found in synteny structure for HMM: {hmm_group}")



def filterFASTAbySyntenyStructure(synteny_structure: str,
                                  input_fasta: Path,
                                  input_hmms: list[Path],
                                  unordered: bool = False,
                                  hmm_meta: Path = None,
                                  hmmer_output_dir: Path = None,
                                  reuse_hmmer_results: bool = True,
                                  method: str = 'hmmsearch',
                                  processes: int = None,
                                  additional_args: list[str] = None) -> SyntenyHits:
    """
    Generate protein-specific database by filtering sequence database
    to only contain sequences which satisfy the provided (gene/hmm)
    structure
    
    @Param
    unordered: bool, whether HMM hits should follow the exact order
    displayed in the synteny structure string or not, i.e., whether
    to search for only synteny (colocation) or collinearity as well
    (same order).
    @Param
    additional_args: additional arguments to hmmsearch or hmmscan. Each
    element in the list is a string with additional arguments for each 
    input hmm (arranged in the same order), an element can also take a 
    value of None to avoid passing additional arguments for a specific 
    input hmm. A single string may also be passed, in which case the 
    same additional argument is passed to hmmsearch for all input hmms
    """
    if hmmer_output_dir is None:
        hmmer_output_dir = Path(input_fasta.parent) / "hmmer_outputs"

    if additional_args is None:
        additional_args = [None for _ in input_hmms]
    
    if type(additional_args) == str:
        additional_args = [additional_args for _ in input_hmms]

    elif type(additional_args) == list:
        if len(additional_args) == 1:
            additional_args = [additional_args[0] for _ in input_hmms]

        if (len(additional_args) > 1) and (len(additional_args) < len(input_hmms)):
            logger.error("Provided additional argument strings are less than the number of input hmms.")
    else:
        logger.error("Additional arguments must be: 1) a list[str], 2) a str, or 3) None")
        sys.exit(1)
    if not os.path.isdir(hmmer_output_dir):
        os.mkdir(hmmer_output_dir)

    logger.info('Running Hmmer')
    hmmer = HMMER(
        input_hmms=input_hmms,
        input_data=input_fasta,
        hmm_output_dir=hmmer_output_dir,
        processes=processes,
        additional_args=additional_args
    )
    hmm_hits = hmmer.getHMMERtables(
        reuse_hmmer_results=reuse_hmmer_results,
        method=method
    )
    logger.info('Filtering results by synteny structure')
    syntenyfilter = SyntenyHMMfilter(hmm_hits, synteny_structure, unordered=unordered)
    hits_by_contig = syntenyfilter.filterHitsBySyntenyStructure()
    if hmm_meta is not None:
        return SyntenyHits.fromHitsDict(hits_by_contig).addHMMmetaInfoToHits(hmm_meta)
    else:
        return SyntenyHits.fromHitsDict(hits_by_contig)