#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to parse synteny structure and filter HMM hits
"""

from __future__ import annotations
import os
import sys
import logging
from pathlib import Path

import pandas as pd

from pynteny.src.preprocessing import FASTA
from pynteny.src.utils import setDefaultOutputPath
from pynteny.src.hmm import HMMER, PGAP

logger = logging.getLogger(__name__)

class LabelParser():
    """
    Parse entry label to extract coded info
    """
    @staticmethod
    def parse(label: str) -> dict:
        """
        Parse sequence labels to obtain contig and locus info
        """
        parsed_dict = {
            'full': label,
            'gene_id': '',
            'contig': '',
            'gene_pos': None,
            'locus_pos': None,
            'strand': ''
            } 
        try: 
            entry = label.split('__')[0]
            meta = label.split('__')[1]
            strand = meta.split('_')[-1]
            locus_pos = tuple([int(pos) for pos in meta.split('_')[-3:-1]])
            gene_pos = int(meta.split('_')[-4])
            contig = '_'.join(meta.split('_')[:-4])

            parsed_dict['gene_id'] = entry
            parsed_dict['contig'] = contig
            parsed_dict['gene_pos'] = gene_pos 
            parsed_dict['locus_pos'] = locus_pos
            parsed_dict['strand'] = strand
        except Exception:
            pass
        return parsed_dict
    
    @staticmethod
    def parse_from_list(labels=list) -> pd.DataFrame: 
        """
        Parse labels in list of labels and return DataFrame
        """
        return pd.DataFrame(
            [LabelParser.parse(label) for label in labels]
        )


class SyntenyParser():
    """
    Tools to parse synteny structure strings
    """
    @staticmethod
    def splitStrandFromLocus(locus_str: str,
                             parsed_symbol: bool = True) -> tuple[str]:
        locus_str = locus_str.strip()
        if locus_str[0] == '<' or locus_str[0] == '>':
            sense = locus_str[0]
            locus_str = locus_str[1:]
            if parsed_symbol:
                strand = 'pos' if sense == '>' else 'neg'
            else:
                strand = sense
        else:
            strand = ""
        return (strand, locus_str)

    @staticmethod
    def getHMMgroupsInStructure(synteny_structure: str) -> list[str]:
        """
        Get hmm names employed in synteny structure,
        if more than one hmm for the same gene, return
        a list with all of them.
        """
        links = synteny_structure.replace("(","").replace(")","").strip().split()
        if not links:
            logger.error("Invalid format for synteny structure")
            sys.exit(1)
        hmm_groups = [
            SyntenyParser.splitStrandFromLocus(h)[1] 
            for h in links if not h.isdigit()
            ]
        return hmm_groups

    @staticmethod
    def getGeneSymbolsInStructure(synteny_structure: str) -> list[str]:
        """
        Get gene symbols employed in synteny structure
        """
        links = synteny_structure.strip().split()
        if not links:
            logger.error("Invalid format for synteny structure")
            sys.exit(1)
        gene_symbols = [
            SyntenyParser.splitStrandFromLocus(h)[1] 
            for h in links if not h.isdigit()
            ]
        return gene_symbols

    @staticmethod
    def getAllHMMsInStructure(synteny_structure: str) -> list[str]:
        """
        Get hmm names employed in synteny structure,
        if more than one hmm for the same gene, return
        a list with all of them.
        """
        hmm_groups = SyntenyParser.getHMMgroupsInStructure(synteny_structure)
        hmm_names = [hmm for hmm_group in hmm_groups for hmm in hmm_group.split("|")]
        return hmm_names

    @staticmethod
    def getStrandsInStructure(synteny_structure: str,
                              parsed_symbol: bool = True) -> list[str]:
        """
        Get strand sense list in structure
        """
        links = synteny_structure.strip().split()
        if not links:
            logger.error("Invalid format for synteny structure")
            sys.exit(1)
        return [
            SyntenyParser.splitStrandFromLocus(h, parsed_symbol)[0] 
            for h in links if not h.isdigit()
            ]

    @staticmethod
    def getMaximumDistancesInStructure(synteny_structure: str) -> list[int]:
        """
        Get maximum gene distances in structure
        """
        links = synteny_structure.strip().split()
        if not links:
            logger.error("Invalid format for synteny structure")
            sys.exit(1)
        return [int(dist) for dist in links if dist.isdigit()]

    @staticmethod
    def parseSyntenyStructure(synteny_structure: str) -> dict:
        """
        Parse synteny structure string. A synteny structure
        is a string like the following:

        >hmm_a n_ab <hmm_b n_bc hmm_c

        where '>' indicates a hmm target located on the positive strand,
        '<' a target located on the negative strand, and n_ab cooresponds 
        to the maximum number of genes separating matched gene a and b. 
        Multiple hmms may be employed (limited by computational capabilities).
        No order symbol in a hmm indicates that results should be independent
        of strand location.
        """
        max_dists = SyntenyParser.getMaximumDistancesInStructure(synteny_structure)
        hmm_groups = SyntenyParser.getHMMgroupsInStructure(synteny_structure)
        strands = SyntenyParser.getStrandsInStructure(synteny_structure)
        return {"hmm_groups": hmm_groups, "strands": strands, "distances": max_dists}

    @staticmethod
    def parseGenesInSyntenyStructure(synteny_structure: str, hmm_meta: Path) -> tuple[str,dict]:
        """
        Convert gene-based synteny structure into a HMM-based one.
        If a gene symbol matches more than one HMM, return a HMM group
        like: (HMM1 | HMM2 | ...)
        """
        pgap = PGAP(hmm_meta)
        # links = synteny_structure.strip().split()
        # if not links:
        #     logger.error("Invalid format for synteny structure")
        #     sys.exit(1)
        # gene_symbols = [
        #     SyntenyParser.splitStrandFromLocus(h)[1] 
        #     for h in links if not h.isdigit()
        #     ]
        gene_symbols = SyntenyParser.getGeneSymbolsInStructure(
            synteny_structure
            )
        strand_locs = SyntenyParser.getStrandsInStructure(
            synteny_structure, parsed_symbol=False
            )
        gene_dists = SyntenyParser.getMaximumDistancesInStructure(
            synteny_structure
            )
        # hmm_names = {
        #     gene_symbol: pgap.getHMMnamesByGeneSymbol(gene_symbol)
        #     for gene_symbol in gene_symbols
        # }
        # unmatched_genes = [
        #     gene_id for gene_id, hmms in hmm_names.items()
        #     if not hmms
        # ]
        hmm_groups = {
            gene_symbol: pgap.getHMMgroupForGeneSymbol(gene_symbol)
            for gene_symbol in gene_symbols
        }
        unmatched_genes = [
            gene_id for gene_id, hmm_group in hmm_groups.items()
            if not hmm_group
        ]
        if unmatched_genes:
            logger.error(
                f"These genes did not get a HMM match in database: {unmatched_genes}"
                )
            sys.exit(1)
      
        hmm_synteny_struc = ""
        # for strand, dist, hmms in zip(
        #     strand_locs, [""] + gene_dists, hmm_names.values()
        #     ):
        #     if len(hmms) == 1:
        #         hmm_group = f"{dist} {strand}{hmms.pop()} "
        #     else:
        #         hmm_group = f"{dist} {strand}({'|'.join(hmms)}) "
        #     hmm_synteny_struc += hmm_group

        for strand, dist, hmm_group in zip(
            strand_locs, [""] + gene_dists, hmm_groups.values()
            ):
            hmm_synteny_struc += f"{dist} {strand}{hmm_group} "

        return hmm_synteny_struc.strip(), hmm_groups



class SyntenyPatternFilters():
    def __init__(self, synteny_structure: str) -> None:
        parsed_structure = SyntenyParser.parseSyntenyStructure(synteny_structure)
        hmm_order_dict = dict(
            zip(parsed_structure["hmm_groups"], range(len(parsed_structure["hmm_groups"])))
            )
        hmm_codes = list(hmm_order_dict.values())
        self.hmm_code_order_pattern = hmm_codes
        parsed_distances = [float("inf")] + parsed_structure["distances"]
        self.distance_order_pattern = [dist + 1 for dist in parsed_distances]
        self.strand_order_pattern = list(map(
            lambda strand : -1 if strand == "neg" else (1 if strand ==  "pos" else 0),
            parsed_structure["strands"]
            ))

    def contains_hmm_pattern(self, data: pd.Series) -> int:
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
    def __init__(self, hmm_hits: dict, synteny_structure: str) -> None:
        """
        Search for contigs that satisfy the given gene synteny structure

        @param: hmm_hits, a dict of pandas DataFrames, as output by
                parseHMMsearchOutput with keys corresponding to hmm names
        """
        self._hmm_hits = hmm_hits
        self._hmms = list(hmm_hits.keys())
        self._synteny_structure = synteny_structure
        self._parsed_structure = SyntenyParser.parseSyntenyStructure(self._synteny_structure)
        self._hmm_order_dict = dict(
            zip(self._parsed_structure["hmm_groups"], range(len(self._parsed_structure["hmm_groups"])))
            )
        self._n_hmm_groups = len(self._hmm_order_dict)

    def mergeHitsByHMMgroup(self, hits: pd.DataFrame):
        g = hits.groupby(["full"]).groups
        _, group_idxs = list(g.keys()), list(g.values())
        for group in group_idxs:
            hmm_names = set(hits.loc[group, "hmm"].values)
            candidate_hmm_group = [
                hmm_group_str for hmm_group_str in self._parsed_structure["hmm_groups"]
                if hmm_names == set(hmm_group_str.split("|"))
                ]
            if candidate_hmm_group:
                hmm_group = candidate_hmm_group[0]
            else:
                hmm_group = "discard"
            hits.loc[group, "hmm"] = hmm_group
        hits = hits[hits.hmm != "discard"].drop_duplicates()
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
        return self._addMetaInfoToHMMhits(all_hit_labels_merged)
    
    def _assignCodeToHMM(self, hmm_name: str):
        code = [
            code for hmm_group, code in self._hmm_order_dict.items()
            if hmm_name in hmm_group
            ]
        if len(code) > 1:
            logger.error(
            f"HMM: {hmm_name} found in more than one hmm group in synteny structure"
            )
            sys.exit(1)
        return code[0]

    def _addMetaInfoToHMMhits(self, all_hit_labels: pd.Dataframe) -> pd.Dataframe:
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

    def _writeAllHitsToTSV(self, all_matched_hits: dict,
                           output_file: str, hmm_meta: Path = None) -> None:
        """
        Write hits matching synteny structure to TSV file
        """
        if hmm_meta is not None:
            pgap = PGAP(hmm_meta)
            header = (
                "contig\tgene_id\tgene_number\tlocus\tstrand\tfull_label\t"
                "HMM\tgene_symbol\tlabel\tproduct\tec_number\n"
                )
        else:
            header = "contig\tgene_id\tgene_number\tlocus\tstrand\tfull_label\tHMM\n"
        output_lines = []
        for contig, matched_hits in all_matched_hits.items():
            for hmm, labels in matched_hits.items():
                if hmm_meta is not None:
                    hmm_meta_info = pgap.getMetaInfoForHMM(hmm)
                    meta_str = "\t".join(
                        [str(v) for k, v in hmm_meta_info.items() if k != "#ncbi_accession"]
                        ).replace("nan", "")
                else:
                    meta_str = ""
                for label in labels:
                    parsed_label = LabelParser.parse(label)
                    output_lines.append(
                        (
                            f"{parsed_label['contig']}\t{parsed_label['gene_id']}\t"
                            f"{parsed_label['gene_pos']}\t{parsed_label['locus_pos']}\t"
                            f"{parsed_label['strand']}\t{parsed_label['full']}\t{hmm}\t{meta_str}\n"
                            )
                    )
        with open(output_file, "w") as outfile:
            outfile.write(header)
            outfile.writelines(output_lines)
        # sort values by gene pos and contig
        pd.read_csv(output_file, sep="\t").sort_values(["contig", "gene_number"]).to_csv(output_file, sep="\t")

    def filterHitsBySyntenyStructure(self, output_tsv: str = None,
                                     hmm_meta: Path = None) -> dict:
        """
        Search for contigs that satisfy the given gene synteny structure
        @param: synteny_structure, a str describing the desired synteny structure,
                structured as follows:

                'hmm_a N_ab hmm_b'

                where N_ab corresponds to the maximum number of genes separating 
                gene found by hmm_a and gene found by hmm_b, and hmm_ corresponds 
                to the name of the hmm as provided in the keys of hmm_hits.
                More than two hmms can be concatenated.
        """
        all_matched_hits = {}
        filters = SyntenyPatternFilters(self._synteny_structure)
        all_hit_labels = self.getAllHMMhits()
        contig_names = all_hit_labels.contig.unique()

        for contig in contig_names:

            matched_hit_labels = {hmm_group: [] for hmm_group in self._parsed_structure["hmm_groups"]}
            contig_hits = all_hit_labels[all_hit_labels.contig == contig].reset_index(drop=True)

            if len(contig_hits.hmm.unique()) == self._n_hmm_groups:
                
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
                    for label, hmm in zip((matched_hits.full.values), matched_hits.hmm):
                        matched_hit_labels[hmm].append(label)

                all_matched_hits[contig] = matched_hit_labels

        if output_tsv is not None:
            self._writeAllHitsToTSV(all_matched_hits, output_file=output_tsv, hmm_meta=hmm_meta)
        return all_matched_hits



def filterFASTABySyntenyStructure(synteny_structure: str,
                                  input_fasta: Path,
                                  input_hmms: list[Path],
                                  hmm_meta: Path = None,
                                  output_dir: Path = None,
                                  output_prefix: str = None,
                                  hmmer_output_dir: Path = None,
                                  reuse_hmmer_results: bool = True,
                                  method: str = 'hmmsearch',
                                  additional_args: list[str] = None) -> None:
    """
    Generate protein-specific database by filtering sequence database
    to only contain sequences which satisfy the provided (gene/hmm)
    structure
    
    @Arguments:
    additional_args: additional arguments to hmmsearch or hmmscan. Each
    element in the list is a string with additional arguments for each 
    input hmm (arranged in the same order), an element can also take a 
    value of None to avoid passing additional arguments for a specific 
    input hmm. A single string may also be passed, in which case the 
    same additional argument is passed to hmmsearch for all input hmms
    """
    if hmmer_output_dir is None:
        hmmer_output_dir = os.path.join(
            setDefaultOutputPath(input_fasta, only_dirname=True), 'hmmer_outputs')
    
    if output_prefix is None:
        output_prefix = ""
        
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
    
    if output_dir is None:
        output_dir = setDefaultOutputPath(input_fasta, only_dirname=True)
    else:
        output_dir = output_dir

    results_table = output_dir / f"{output_prefix}synteny_matched.tsv"

    logger.info('* Running Hmmer')
    hmmer = HMMER(
        input_hmms=input_hmms,
        input_data=input_fasta,
        hmm_output_dir=hmmer_output_dir,
        additional_args=additional_args
    )
    hmm_hits = hmmer.getHMMERtables(
        reuse_hmmer_results=reuse_hmmer_results,
        method=method
    )

    logger.info('* Filtering results by synteny structure')
    syntenyfilter = SyntenyHMMfilter(hmm_hits, synteny_structure)
    hmm_groups = syntenyfilter._parsed_structure["hmm_groups"]
    matches = syntenyfilter.filterHitsBySyntenyStructure(
        output_tsv=results_table,
        hmm_meta=hmm_meta
        )
    # logger.info("* Writing matching sequences to FASTA files")
    # SyntenyParser.parseGenesInSyntenyStructure
    # fasta = FASTA(input_fasta)
    # df = pd.read_csv(results_table, sep="\t")
    # for hmm_group in hmm_groups:
    #     record_ids = df[df.HMM == hmm_group].full_label.values.tolist()
    #     output_fasta = output_dir / f"{output_prefix}{hmm_group}_hits.fasta"
    #     if record_ids:
    #         fasta.filterByIDs(
    #             record_ids=record_ids,
    #             output_file=output_fasta      
    #         )
    #     else:
    #         logger.warning(f"No record matches found in synteny structure for HMM: {hmm_group}")

def writeMatchingSequencesToFASTA(input_fasta: Path,
                                  synteny_results: Path,
                                  output_dir: Path,
                                  output_prefix: str,
                                  hmm_group_to_gene: dict = None
                                  ) -> None:
    """
    Write matching sequences to FASTA files
    """
    SyntenyParser.parseGenesInSyntenyStructure
    fasta = FASTA(input_fasta)
    df = pd.read_csv(synteny_results, sep="\t")
    hmm_groups = df.HMM.unique().tolist()

    for hmm_group in hmm_groups:
        record_ids = df[df.HMM == hmm_group].full_label.values.tolist()
        if hmm_group_to_gene is not None:
            gene_symbol = hmm_group_to_gene[hmm_group] if hmm_group in hmm_group_to_gene else "no_gene_symbol"
            output_fasta = output_dir / f"{output_prefix}{gene_symbol}_{hmm_group}_hits.fasta"
        else:
            output_fasta = output_dir / f"{output_prefix}{hmm_group}_hits.fasta"
        if record_ids:
            fasta.filterByIDs(
                record_ids=record_ids,
                output_file=output_fasta      
            )
        else:
            logger.warning(f"No record matches found in synteny structure for HMM: {hmm_group}")