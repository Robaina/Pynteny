#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to filter HMM hits by synteny structure
"""

from __future__ import annotations

import logging
import os
import sys
from pathlib import Path

import pandas as pd

import pynteny.parsers.labelparser as labelparser
import pynteny.parsers.syntenyparser as syntenyparser
from pynteny.hmm import HMMER, PGAP
from pynteny.preprocessing import FASTA

logger = logging.getLogger(__name__)


class SyntenyPatternFilters:
    """Methods to filter hmm hits in the  same contig by synteny structure
    or collinearity. These filters are inputs to pandas.Dataframe.rolling method.
    """

    def __init__(self, synteny_structure: str, unordered: bool = False) -> None:
        """Initialize filter class from synteny structure.

        Args:
            synteny_structure (str): a str describing the desired synteny structure,
                structured as follows:

                '>hmm_a N_ab hmm_b bc <hmm_c'

                where N_ab corresponds to the maximum number of genes separating
                gene found by hmm_a and gene found by hmm_b, and hmm_ corresponds
                to the name of the hmm as provided in the keys of hmm_hits.
                More than two hmms can be concatenated. Strand location may be
                specified by using '>' for sense and '<' for antisense.
            unordered (bool, optional): whether the HMMs should be arranged in the
                exact same order displayed in the synteny_structure or in
                any order. If ordered, the filters would filter collinear rather
                than syntenic structures. Defaults to False.
        """
        parsed_structure = syntenyparser.parse_synteny_structure(synteny_structure)
        hmm_codes = list(range(len(parsed_structure["hmm_groups"])))
        self.hmm_code_order_pattern = hmm_codes

        if unordered:
            max_distance = max(parsed_structure["distances"])
            max_distances = [float("inf")] + [
                max_distance for _ in parsed_structure["distances"]
            ]
        else:
            max_distances = [float("inf")] + parsed_structure["distances"]
        min_distances = [-float("inf")] + [0 for _ in parsed_structure["distances"]]
        self.distance_order_pattern = [
            (min_dist, max_dist + 1)
            for min_dist, max_dist in zip(min_distances, max_distances)
        ]

        self.strand_order_pattern = list(
            map(
                lambda strand: -1 if strand == "neg" else (1 if strand == "pos" else 0),
                parsed_structure["strands"],
            )
        )
        self._unordered = unordered

    def contains_hmm_pattern(self, data: pd.Series) -> int:
        """Check if series items contain a profile HMM

        Args:
            data (pd.Series): a series resulting from calling rolling on a pandas column.

        Returns:
            int: 1 for True 0 for False.
        """
        if self._unordered:
            return 1 if set(data.values) == set(self.hmm_code_order_pattern) else 0
        else:
            return 1 if data.values.tolist() == self.hmm_code_order_pattern else 0

    def contains_distance_pattern(self, data: pd.Series) -> int:
        """Check if series items satisfy the maximum distance
           between HMM hits.

        Args:
            data (pd.Series): a series resulting from calling rolling on a pandas column.

        Returns:
            int: 1 for True 0 for False.
        """
        return (
            1
            if all(
                [
                    (data_dist <= max_dist) and (data_dist > min_dist)
                    for data_dist, (min_dist, max_dist) in zip(
                        data.values.tolist(), self.distance_order_pattern
                    )
                ]
            )
            else 0
        )

    def contains_strand_pattern(self, data: pd.Series) -> int:
        """Check if series items satisfy the strand pattern
           between HMM hits.

        Args:
            data (pd.Series): a series resulting from calling rolling on a pandas column.

        Returns:
            int: 1 for True 0 for False.
        """
        strand_comparisons = []
        for data_strand, pattern_strand in zip(data.values, self.strand_order_pattern):
            if pattern_strand != 0:
                strand_comparisons.append(data_strand == pattern_strand)
            else:
                strand_comparisons.append(True)
        return 1 if all(strand_comparisons) else 0


class SyntenyHMMfilter:
    """Tools to search for synteny structures among sets of hmm models"""

    def __init__(
        self, hmm_hits: dict, synteny_structure: str, unordered: bool = True
    ) -> None:
        """Search for contigs that satisfy the given gene synteny structure.

        Args:
            hmm_hits (dict): a dict of pandas DataFrames, as output by
                parseHMMsearchOutput with keys corresponding to hmm names
            synteny_structure (str): a str describing the desired synteny structure,
                structured as follows:

                '>hmm_a N_ab hmm_b bc <hmm_c'

                where N_ab corresponds to the maximum number of genes separating
                gene found by hmm_a and gene found by hmm_b, and hmm_ corresponds
                to the name of the hmm as provided in the keys of hmm_hits.
                More than two hmms can be concatenated. Strand location may be
                specificed by using '>' for sense and '<' for antisense.
            unordered (bool, optional): whether the HMMs should be arranged in the
                exact same order displayed in the synteny_structure or in
                any order. If ordered, the filters would filter collinear rather
                than syntenic structures. Defaults to True.
        """
        self._unordered = unordered
        self._hmm_hits = hmm_hits
        self._hmms = list(hmm_hits.keys())
        self._synteny_structure = synteny_structure
        self._contains_hmm_groups = syntenyparser.contains_HMM_groups(
            self._synteny_structure
        )
        self._parsed_structure = syntenyparser.parse_synteny_structure(
            self._synteny_structure
        )
        if self._unordered:
            self._parsed_structure["strands"] = [
                "" for _ in self._parsed_structure["strands"]
            ]
            logger.warning(
                "Disregarding strand location as unordered structured selected"
            )
        self._hmm_order_dict = dict(
            zip(
                self._parsed_structure["hmm_groups"],
                range(len(self._parsed_structure["hmm_groups"])),
            )
        )
        self._n_hmm_groups = len(self._hmm_order_dict)
        self._unordered = unordered

    def _assign_code_to_HMM(self, hmm_group_name: str) -> int:
        hmm_group_name = str(hmm_group_name)
        code = [
            code
            for hmm_group, code in self._hmm_order_dict.items()
            if set(hmm_group_name.split("|")).issubset(set(hmm_group.split("|")))
        ]
        if len(code) > 1:
            logger.error(
                f"HMM: {hmm_group_name} found in more than one hmm group in synteny structure"
            )
            sys.exit(1)
        if not code:
            logger.error(f"No code found for HMM: {hmm_group_name}")
            sys.exit(1)
        return code[0]

    def _add_meta_codes_to_HMM_hits(self, all_hit_labels: pd.Dataframe) -> pd.Dataframe:
        """Add numeric codes for each hmm and strand, compute distance between genes"""
        all_hit_labels["gene_pos_diff"] = all_hit_labels.gene_pos.diff()
        all_hit_labels.loc[
            0, "gene_pos_diff"
        ] = 1  # required for rolling (skips first nan)
        all_hit_labels["hmm_code"] = all_hit_labels.hmm.apply(self._assign_code_to_HMM)
        all_hit_labels["strand"] = all_hit_labels.strand.apply(
            lambda strand: -1 if strand == "neg" else (1 if strand == "pos" else 0)
        )
        return all_hit_labels

    def _merge_hits_by_HMM_group(self, hits: pd.DataFrame):
        """Merge hit sequences within an HMM group representing
        the same gene symbol. Drop duplicated hits within
        each HMM group.
        """
        groups = hits.groupby(["full_label"]).groups
        for group_idxs in groups.values():
            hmm_names = set(hits.loc[group_idxs, "hmm"].values)
            candidate_hmm_group = [
                hmm_group_str
                for hmm_group_str in self._parsed_structure["hmm_groups"]
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

    def get_all_HMM_hits(self) -> pd.DataFrame:
        """Group and preprocess all hit labels into a single dataframe.

        Returns:
            pd.DataFrame: HMMER3 hit labels matching provided HMMs.
        """
        hit_labels = {}
        for hmm, hits in self._hmm_hits.items():
            labels = hits.id.values.tolist()
            if not labels:
                logger.error(f"No records found in database matching HMM: {hmm}")
                sys.exit(1)
            hit_labels[hmm] = labelparser.parse_from_list(labels)
            hit_labels[hmm]["hmm"] = hmm
        # Create single dataframe with new column corresponding to HMM and all hits
        # Remove contigs with less hits than the number of hmms in synteny structure (drop if enabling partial searching)
        all_hit_labels = (
            pd.concat(hit_labels.values())
            .groupby("contig")
            .filter(lambda x: len(x) >= self._n_hmm_groups)
            .sort_values(["contig", "gene_pos"], ascending=True)
        )
        all_hit_labels = all_hit_labels.reset_index(drop=True)
        if self._contains_hmm_groups:
            all_hit_labels = self._merge_hits_by_HMM_group(all_hit_labels)
        all_hit_labels = all_hit_labels.reset_index(drop=True)
        return self._add_meta_codes_to_HMM_hits(all_hit_labels)

    def filter_hits_by_synteny_structure(self) -> dict:
        """Search for contigs that satisfy the given gene synteny structure.

        Returns:
            dict: HMMER3 hits separated by contig.
        """
        all_matched_hits = {}
        filters = SyntenyPatternFilters(
            self._synteny_structure, unordered=self._unordered
        )
        all_hit_labels = self.get_all_HMM_hits()
        if all_hit_labels.full_label.duplicated().any():
            logger.warning(
                "At least two different HMMs produced identical sequence hits"
            )
        if all_hit_labels.full_label.duplicated().sum() == all_hit_labels.shape[0] / 2:
            logger.error(
                (
                    "All input profile HMMs rendered identical sequence hits. "
                    "Inspect HMMER output tables to evaluate hits."
                )
            )
            sys.exit(1)
        contig_names = all_hit_labels.contig.unique()

        for contig in contig_names:
            contig_hits = all_hit_labels[all_hit_labels.contig == contig].reset_index(
                drop=True
            )
            matched_hit_labels = {
                hmm_group: [] for hmm_group in contig_hits.hmm.unique()
            }

            if contig_hits.hmm.nunique() >= self._n_hmm_groups:
                hmm_match = contig_hits.hmm_code.rolling(
                    window=self._n_hmm_groups
                ).apply(filters.contains_hmm_pattern)
                strand_match = contig_hits.strand.rolling(
                    window=self._n_hmm_groups
                ).apply(filters.contains_strand_pattern)
                if self._n_hmm_groups > 1:
                    distance_match = contig_hits.gene_pos_diff.rolling(
                        window=self._n_hmm_groups
                    ).apply(filters.contains_distance_pattern)
                    matched_rows = contig_hits[
                        (hmm_match == 1) & (strand_match == 1) & (distance_match == 1)
                    ]
                else:
                    matched_rows = contig_hits[(hmm_match == 1) & (strand_match == 1)]
                for i in matched_rows.index:
                    matched_hits = contig_hits.iloc[
                        i - (self._n_hmm_groups - 1) : i + 1, :
                    ]
                    for label, hmm in zip(
                        (matched_hits.full_label.values), matched_hits.hmm
                    ):
                        matched_hit_labels[hmm].append(label)

                all_matched_hits[contig] = matched_hit_labels
        return all_matched_hits


class SyntenyHits:
    """Store and manipulate synteny hits by contig"""

    def __init__(self, synteny_hits: pd.DataFrame) -> None:
        """
        Initialize from synteny hits object.
        """
        self._synteny_hits = synteny_hits.drop_duplicates()

    @staticmethod
    def _hits_to_dataframe(hits_by_contig: dict) -> pd.DataFrame:
        """Return synteny hits as a dataframe"""
        data = []
        columns = [
            "contig",
            "gene_id",
            "gene_number",
            "locus",
            "strand",
            "full_label",
            "hmm",
        ]
        for contig, matched_hits in hits_by_contig.items():
            for hmm, labels in matched_hits.items():
                for label in labels:
                    parsed_label = labelparser.parse(label)
                    data.append(
                        [
                            parsed_label["contig"],
                            parsed_label["gene_id"],
                            parsed_label["gene_pos"],
                            parsed_label["locus_pos"],
                            parsed_label["strand"],
                            parsed_label["full_label"],
                            hmm,
                        ]
                    )
        return pd.DataFrame(data, columns=columns).sort_values(
            ["contig", "gene_number"]
        )

    @classmethod
    def from_hits_dict(cls, hits_by_contig: dict) -> SyntenyHits:
        """Initialize SyntenyHits object from hits_by_contig dictionary.

        Args:
            hits_by_contig (dict): HMMER3 hit labels separated by contig name.

        Returns:
            SyntenyHits: initialized object of class SyntenyHits.
        """
        return cls(cls._hits_to_dataframe(hits_by_contig))

    @property
    def hits(self) -> pd.DataFrame:
        """Return synteny hits.

        Returns:
            pd.DataFrame: Synteny hits as dataframe.
        """
        return self._synteny_hits

    def add_HMM_meta_info_to_hits(self, hmm_meta: Path) -> SyntenyHits:
        """Add molecular metadata to synteny hits.

        Args:
            hmm_meta (Path): path to PGAP metadata file.

        Returns:
            SyntenyHits: and instance of class SyntenyHits.
        """
        hmm_meta = Path(hmm_meta)
        fields = ["gene_symbol", "label", "product", "ec_number"]
        if all([f in self._synteny_hits.columns for f in fields]):
            return self._synteny_hits
        pgap = PGAP(hmm_meta)
        self._synteny_hits[fields] = ""
        for row in self._synteny_hits.itertuples():
            meta_values = [
                [
                    str(v).replace("nan", "")
                    for k, v in pgap.get_meta_info_for_HMM(hmm).items()
                    if k != "#ncbi_accession"
                ]
                for hmm in row.hmm.split("|")
            ]
            self._synteny_hits.loc[row.Index, fields] = [
                "|".join(v) for v in zip(*meta_values)
            ]
        return SyntenyHits(self._synteny_hits)

    def write_to_TSV(self, output_tsv: Path) -> None:
        """Write synteny hits to a TSV file.

        Args:
            output_tsv (Path): path to output tsv file.
        """
        output_tsv = Path(output_tsv)
        self._synteny_hits.to_csv(output_tsv, sep="\t", index=False)

    def write_hit_sequences_to_FASTA_files(
        self,
        sequence_database: Path,
        output_dir: Path,
        output_prefix: str = None,
    ) -> None:
        """Write matching sequences to FASTA files.

        Args:
            sequence_database (Path): path to the peptide or nucleotide sequence database
                in which the synteny search was conducted.
            output_dir (Path): path to output directory.
            output_prefix (str, optional): prefix for output files. Defaults to None.
        """
        sequence_database = Path(sequence_database)
        output_dir = Path(output_dir)
        fasta = FASTA(sequence_database)
        hmm_groups = self._synteny_hits.hmm.unique().tolist()

        for hmm_group in hmm_groups:
            record_ids = self._synteny_hits[
                self._synteny_hits.hmm == hmm_group
            ].full_label.values.tolist()

            if (
                "gene_symbol" in self._synteny_hits.columns
                and "label" in self._synteny_hits.columns
            ):
                gene_symbol = self._synteny_hits[
                    self._synteny_hits.hmm == hmm_group
                ].gene_symbol.unique()[0]
                gene_label = self._synteny_hits[
                    self._synteny_hits.hmm == hmm_group
                ].label.unique()[0]
                gene_id = (
                    gene_symbol
                    if not pd.isna(gene_symbol)
                    else (gene_label if not pd.isna(gene_label) else "")
                )
            else:
                gene_id = ""
            gene_id_str = gene_id + "_" if gene_id else ""

            output_fasta = (
                output_dir
                / f"{output_prefix}{gene_id_str.replace('|', '_')}{hmm_group.replace('|', '_')}_hits.fasta"
            )
            if record_ids:
                fasta.filter_by_IDs(
                    record_ids=record_ids,
                    output_file=output_fasta,
                    point_to_new_file=False,
                )
            else:
                logger.warning(
                    f"No record matches found in synteny structure for HMM: {hmm_group}"
                )


def filter_FASTA_by_synteny_structure(
    synteny_structure: str,
    input_fasta: Path,
    input_hmms: list[Path],
    unordered: bool = False,
    hmm_meta: Path = None,
    hmmer_output_dir: Path = None,
    reuse_hmmer_results: bool = True,
    method: str = "hmmsearch",
    processes: int = None,
    additional_args: list[str] = None,
) -> SyntenyHits:
    """Generate protein-specific database by filtering sequence database
       to only contain sequences which satisfy the provided (gene/hmm) structure.

    Args:
        synteny_structure (str): a str describing the desired synteny structure,
            structured as follows:

            '>hmm_a N_ab hmm_b bc <hmm_c'

            where N_ab corresponds to the maximum number of genes separating
            gene found by hmm_a and gene found by hmm_b, and hmm_ corresponds
            to the name of the hmm as provided in the keys of hmm_hits.
            More than two hmms can be concatenated. Strand location may be
            specificed by using '>' for sense and '<' for antisense.
        input_fasta (Path): input fasta containing sequence database to be searched.
        input_hmms (list[Path]): list containing paths to hmms contained in synteny structure.
        unordered (bool, optional): whether HMM hits should follow the exact order
            displayed in the synteny structure string or not, i.e., whether
            to search for only synteny (colocation) or collinearity as well
            (same order). Defaults to False.
        hmm_meta (Path, optional): path to PGAP's metadata file. Defaults to None.
        hmmer_output_dir (Path, optional): output directory to store HMMER3 output files. Defaults to None.
        reuse_hmmer_results (bool, optional): if True then HMMER3 won't be run again for HMMs already
            searched in the same output directory. Defaults to True.
        method (str, optional): select between 'hmmsearch' or 'hmmscan'. Defaults to 'hmmsearch'.
        processes (int, optional): maximum number of threads to be employed. Defaults to all minus one.
        additional_args (list[str], optional): additional arguments to hmmsearch or hmmscan. Each
            element in the list is a string with additional arguments for each
            input hmm (arranged in the same order), an element can also take a
            value of None to avoid passing additional arguments for a specific
            input hmm. A single string may also be passed, in which case the
            same additional argument is passed to hmmsearch for all input hmms. Defaults to None.

    Returns:
        SyntenyHits: object of class SyntenyHits containing labels matching synteny structure.
    """
    input_fasta = Path(input_fasta)
    if hmmer_output_dir is None:
        hmmer_output_dir = Path(input_fasta.parent) / "hmmer_outputs"
    else:
        hmmer_output_dir = Path(hmmer_output_dir)

    if additional_args is None:
        additional_args = [None for _ in input_hmms]

    # if type(additional_args) == str:
    if isinstance(additional_args, str):
        logger.warning(f"Repeating hmmsearch arg: '{additional_args}' for all HMMs")
        additional_args = [additional_args for _ in input_hmms]

    # elif type(additional_args) == list:
    elif isinstance(additional_args, list):
        if len(additional_args) == 1:
            logger.warning(
                f"Repeating hmmsearch arg: '{additional_args[0]}' for all HMMs"
            )
            additional_args = [additional_args[0] for _ in input_hmms]

        if (len(additional_args) > 1) and (len(additional_args) < len(input_hmms)):
            logger.error(
                "Provided additional argument strings are less than the number of input hmms."
            )
    else:
        logger.error(
            "Additional arguments must be: 1) a list[str], 2) a str, or 3) None"
        )
        sys.exit(1)
    if not os.path.isdir(hmmer_output_dir):
        os.mkdir(hmmer_output_dir)

    logger.info("Running Hmmer")
    hmmer = HMMER(
        input_hmms=input_hmms,
        input_data=input_fasta,
        hmm_output_dir=hmmer_output_dir,
        processes=processes,
        additional_args=additional_args,
    )
    hmm_hits = hmmer.get_HMMER_tables(
        reuse_hmmer_results=reuse_hmmer_results, method=method
    )
    hmms_without_hits = []
    for hmm_name, hmmer_hits in hmm_hits.items():
        if hmmer_hits.empty:
            hmms_without_hits.append(hmm_name)
    if hmms_without_hits:
        for hmm_name in hmms_without_hits:
            logger.error(f"No records found in database matching HMM: {hmm_name}")
        sys.exit(1)

    logger.info("Filtering results by synteny structure")
    syntenyfilter = SyntenyHMMfilter(hmm_hits, synteny_structure, unordered=unordered)
    hits_by_contig = syntenyfilter.filter_hits_by_synteny_structure()
    if hmm_meta is not None:
        return SyntenyHits.from_hits_dict(hits_by_contig).add_HMM_meta_info_to_hits(
            hmm_meta
        )
    else:
        return SyntenyHits.from_hits_dict(hits_by_contig)
