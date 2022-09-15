#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to parse synteny structure
"""

from __future__ import annotations
import sys
import logging
from pathlib import Path

import pandas as pd

from pynteny.src.hmm import PGAP
from pynteny.src.utils import isRightListNestedType

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
            "full_label": label,
            "gene_id": "",
            "contig": "",
            "gene_pos": None,
            "locus_pos": None,
            "strand": ""
            } 
        try: 
            entry = label.split("__")[0]
            meta = label.split("__")[1]
            strand = meta.split("_")[-1]
            locus_pos = tuple([int(pos) for pos in meta.split("_")[-3:-1]])
            gene_pos = int(meta.split("_")[-4])
            contig = '_'.join(meta.split("_")[:-4])

            parsed_dict["gene_id"] = entry
            parsed_dict["contig"] = contig
            parsed_dict["gene_pos"] = gene_pos 
            parsed_dict["locus_pos"] = locus_pos
            parsed_dict["strand"] = strand
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
    def reformatSyntenyStructure(synteny_structure: str) -> str:
        """
        Reformat synteny structure
        """
        synteny_structure = synteny_structure.replace(
            " |", "|").replace("| ", "|").strip()
        return synteny_structure

    @staticmethod
    def isValidStructure(synteny_structure: str) -> bool:
        """
        Validate synteny structure format
        """
        synteny_structure = synteny_structure.replace(
            " |", "|").replace("| ", "|").strip()
        parsed_struc = SyntenyParser.parseSyntenyStructure(synteny_structure)
        right_type = all(
            (isRightListNestedType(parsed_struc["hmm_groups"], str),
            isRightListNestedType(parsed_struc["strands"], str),
            isRightListNestedType(parsed_struc["distances"], int))
        )
        right_format = (
            len(parsed_struc["hmm_groups"]) == len(parsed_struc["strands"]) and
            len(parsed_struc["hmm_groups"]) == (len(parsed_struc["distances"]) + 1)
        )
        return False if (not right_type or not right_format) else True

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
        gene_symbols = SyntenyParser.getGeneSymbolsInStructure(
            synteny_structure
            )
        strand_locs = SyntenyParser.getStrandsInStructure(
            synteny_structure, parsed_symbol=False
            )
        gene_dists = SyntenyParser.getMaximumDistancesInStructure(
            synteny_structure
            )
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

        for strand, dist, hmm_group in zip(
            strand_locs, [""] + gene_dists, hmm_groups.values()
            ):
            if "|" in hmm_group:
                hmm_group = f"({hmm_group})"
            hmm_synteny_struc += f"{dist} {strand}{hmm_group} "

        return hmm_synteny_struc.strip(), hmm_groups