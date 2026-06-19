#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests for the synteny filter, including the paralog cross-hit case
(see issue #96).
"""

import unittest

import pandas as pd

from pynteny.filter import SyntenyHMMfilter


def _paralog_crosshit_hits() -> dict:
    """Build a minimal hmm_hits dict reproducing a paralog cross-hit scenario.

    A textbook syntenic operon on a single contig, all on the positive strand,
    at gene positions 33 -> 34 -> 35. Each peptide has a single fixed label
    (independent of which HMM hits it). The paralogous nifD/nifK models also
    cross-hit each other's peptides, so peptides p34 and p35 each appear in two
    HMM hit tables. The highest-scoring HMM per peptide is exactly the canonical
    nifH -> nifD -> nifK operon.
    """
    p33 = "p33__ctg_33_1_100_pos"  # true nifH peptide
    p34 = "p34__ctg_34_1_100_pos"  # true nifD peptide
    p35 = "p35__ctg_35_1_100_pos"  # true nifK peptide
    nifH = pd.DataFrame({"id": [p33], "bitscore": [452.7]})
    nifD = pd.DataFrame(
        {"id": [p34, p35], "bitscore": [856.5, 95.0]}  # p34 real, p35 cross-hit
    )
    nifK = pd.DataFrame(
        {"id": [p34, p35], "bitscore": [88.0, 861.6]}  # p34 cross-hit, p35 real
    )
    return {"nifH": nifH, "nifD": nifD, "nifK": nifK}


SYNTENY_STRUC = ">nifH 0 >nifD 0 >nifK"


def _matched_labels(matched_hits: dict) -> set:
    """Flatten the {contig: {hmm: [labels]}} structure into a set of labels."""
    return {
        label
        for hmm_to_labels in matched_hits.values()
        for labels in hmm_to_labels.values()
        for label in labels
    }


class TestParalogCrossHits(unittest.TestCase):
    def test_default_drops_operon_on_crosshit(self):
        """Without best_hmm_wins, cross-hits silently drop the operon (issue #96)."""
        synteny_filter = SyntenyHMMfilter(
            _paralog_crosshit_hits(), SYNTENY_STRUC, unordered=False
        )
        labels = _matched_labels(synteny_filter.filter_hits_by_synteny_structure())
        self.assertNotIn(
            "p33__ctg_33_1_100_pos",
            labels,
            "Expected the operon to be dropped under cross-hits with default behaviour",
        )

    def test_best_hmm_wins_recovers_operon(self):
        """With best_hmm_wins, each peptide is assigned its best HMM and the
        canonical operon is recovered."""
        synteny_filter = SyntenyHMMfilter(
            _paralog_crosshit_hits(),
            SYNTENY_STRUC,
            unordered=False,
            best_hmm_wins=True,
        )
        labels = _matched_labels(synteny_filter.filter_hits_by_synteny_structure())
        self.assertTrue(
            {
                "p33__ctg_33_1_100_pos",
                "p34__ctg_34_1_100_pos",
                "p35__ctg_35_1_100_pos",
            }.issubset(labels),
            "Expected best_hmm_wins to recover the full syntenic operon",
        )


if __name__ == "__main__":
    unittest.main()
