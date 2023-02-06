#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to parse synteny structure strings
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from pynteny.hmm import PGAP
from pynteny.utils import is_right_list_nested_type

logger = logging.getLogger(__name__)


def reformat_synteny_structure(synteny_structure: str) -> str:
    """Remove illegal symbols and extra white space."""
    synteny_structure = synteny_structure.replace(" |", "|").replace("| ", "|").strip()
    return synteny_structure


def contains_HMM_groups(synteny_structure: str) -> bool:
    """Check whether structure contains groups of gene-equivalent HMMs."""
    return "|" in synteny_structure


def is_valid_structure(synteny_structure: str) -> bool:
    """Validate synteny structure format."""
    synteny_structure = synteny_structure.replace(" |", "|").replace("| ", "|").strip()
    parsed_struc = parse_synteny_structure(synteny_structure)
    right_type = all(
        (
            is_right_list_nested_type(parsed_struc["hmm_groups"], str),
            is_right_list_nested_type(parsed_struc["strands"], str),
            is_right_list_nested_type(parsed_struc["distances"], int),
        )
    )
    right_format = len(parsed_struc["hmm_groups"]) == len(
        parsed_struc["strands"]
    ) and len(parsed_struc["hmm_groups"]) == (len(parsed_struc["distances"]) + 1)
    return False if (not right_type or not right_format) else True


def split_strand_from_locus(locus_str: str, parsed_symbol: bool = True) -> tuple[str]:
    """Split strand info from locus tag / HMM model.

    Args:
        locus_str (str): a substring of a synteny structure containing
            a gene symbol / HMM name and strand info.
        parsed_symbol (bool, optional): if True, strand info '>' is parsed
            as 'pos' and '<' as 'neg'. Defaults to True.

    Returns:
        tuple[str]: tuple with parsed strand info and gene symbol / HMM name.
    """
    locus_str = locus_str.strip()
    if locus_str[0] == "<" or locus_str[0] == ">":
        sense = locus_str[0]
        locus_str = locus_str[1:]
        if parsed_symbol:
            strand = "pos" if sense == ">" else "neg"
        else:
            strand = sense
    else:
        strand = ""
    return (strand, locus_str)


def get_HMM_groups_in_structure(synteny_structure: str) -> list[str]:
    """Get hmm names employed in synteny structure,
    if more than one hmm for the same gene, return
    a list with all of them.
    """
    links = synteny_structure.replace("(", "").replace(")", "").strip().split()
    if not links:
        logger.error("Invalid format for synteny structure")
        sys.exit(1)
    hmm_groups = [split_strand_from_locus(h)[1] for h in links if not h.isdigit()]
    return hmm_groups


def get_gene_symbols_in_structure(synteny_structure: str) -> list[str]:
    """Retrieve gene symbols contained in synteny structure."""
    links = synteny_structure.strip().split()
    if not links:
        logger.error("Invalid format for synteny structure")
        sys.exit(1)
    gene_symbols = [split_strand_from_locus(h)[1] for h in links if not h.isdigit()]
    return gene_symbols


def get_all_HMMs_in_structure(synteny_structure: str) -> list[str]:
    """Get hmm names employed in synteny structure,
    if more than one hmm for the same gene, return
    a list with all of them.
    """
    hmm_groups = get_HMM_groups_in_structure(synteny_structure)
    hmm_names = [hmm for hmm_group in hmm_groups for hmm in hmm_group.split("|")]
    return hmm_names


def get_strands_in_structure(
    synteny_structure: str, parsed_symbol: bool = True
) -> list[str]:
    """Get strand sense list in structure.

    Args:
        synteny_structure (str): a synteny structure.
        parsed_symbol (bool, optional): if True, strand info '>' is parsed
            as 'pos' and '<' as 'neg'. Defaults to True.

    Returns:
        list[str]: parsed synteny structure as a list of tuples containing
            HMM name and strand info for each HMM group.
    """
    links = synteny_structure.strip().split()
    if not links:
        logger.error("Invalid format for synteny structure")
        sys.exit(1)
    return [
        split_strand_from_locus(h, parsed_symbol)[0] for h in links if not h.isdigit()
    ]


def get_maximum_distances_in_structure(synteny_structure: str) -> list[int]:
    """Get maximum gene distances in synteny structure."""
    links = synteny_structure.strip().split()
    if not links:
        logger.error("Invalid format for synteny structure")
        sys.exit(1)
    return [int(dist) for dist in links if dist.isdigit()]


def parse_synteny_structure(synteny_structure: str) -> dict:
    """Parse synteny structure string.

    Args:
        synteny_structure (str): a string like the following:
            >hmm_a n_ab <hmm_b n_bc hmm_c

            where '>' indicates a hmm target located on the positive strand,
            '<' a target located on the negative strand, and n_ab cooresponds
            to the maximum number of genes separating matched gene a and b.
            Multiple hmms may be employed (limited by computational capabilities).
            No order symbol in a hmm indicates that results should be independent
            of strand location.

    Returns:
        dict: parsed synteny structure.
    """
    max_dists = get_maximum_distances_in_structure(synteny_structure)
    hmm_groups = get_HMM_groups_in_structure(synteny_structure)
    strands = get_strands_in_structure(synteny_structure)
    return {"hmm_groups": hmm_groups, "strands": strands, "distances": max_dists}


def parse_genes_in_synteny_structure(
    synteny_structure: str, hmm_meta: Path
) -> tuple[str, dict]:
    """Convert gene-based synteny structure into a HMM-based one.
    If a gene symbol matches more than one HMM, return a HMM group
    like: (HMM1 | HMM2 | ...).

    Args:
        synteny_structure (str): a string like the following:
            >hmm_a n_ab <hmm_b n_bc hmm_c

            where '>' indicates a hmm target located on the positive strand,
            '<' a target located on the negative strand, and n_ab cooresponds
            to the maximum number of genes separating matched gene a and b.
            Multiple hmms may be employed (limited by computational capabilities).
            No order symbol in a hmm indicates that results should be independent
            of strand location.
        hmm_meta (Path): path to PGAP's metadata file.

    Returns:
        tuple[str,dict]: parsed synteny structure where gene symbols are replaced
            by HMM names.
    """
    pgap = PGAP(hmm_meta)
    gene_symbols = get_gene_symbols_in_structure(synteny_structure)
    strand_locs = get_strands_in_structure(synteny_structure, parsed_symbol=False)
    gene_dists = get_maximum_distances_in_structure(synteny_structure)
    hmm_groups = {
        gene_symbol: pgap.get_HMM_group_for_gene_symbol(gene_symbol)
        for gene_symbol in gene_symbols
    }
    unmatched_genes = [
        gene_id for gene_id, hmm_group in hmm_groups.items() if not hmm_group
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
