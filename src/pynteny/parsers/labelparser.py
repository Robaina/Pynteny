#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to parse record labels to extract coded info
"""

from __future__ import annotations

import logging

import pandas as pd

logger = logging.getLogger(__name__)


def parse(label: str) -> dict:
    """Parse sequence labels to obtain contig and locus info

    Args:
        label (str): sequence label

    Returns:
        dict: dictionary with parsed label info.
    """
    parsed_dict = {
        "full_label": label,
        "gene_id": "",
        "contig": "",
        "gene_pos": None,
        "locus_pos": None,
        "strand": "",
    }
    try:
        entry = label.split("__")[0]
        meta = label.split("__")[1]
        strand = meta.split("_")[-1]
        locus_pos = tuple([int(pos) for pos in meta.split("_")[-3:-1]])
        gene_pos = int(meta.split("_")[-4])
        contig = "_".join(meta.split("_")[:-4])

        parsed_dict["gene_id"] = entry
        parsed_dict["contig"] = contig
        parsed_dict["gene_pos"] = gene_pos
        parsed_dict["locus_pos"] = locus_pos
        parsed_dict["strand"] = strand
    except Exception:
        pass
    return parsed_dict


def parse_from_list(labels: list[str]) -> pd.DataFrame:
    """Parse labels in list of labels and return DataFrame.

    Args:
        labels (list, optional): list of labels as stringgs.

    Returns:
        pd.DataFrame: Dataframe containing parsed information from labels.
    """
    return pd.DataFrame([parse(label) for label in labels])
