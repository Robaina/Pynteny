#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path

import pandas as pd
import streamlit as st
from st_aggrid import AgGrid, GridOptionsBuilder


parent_dir = Path(__file__).parent


def plot_dataframe(data: pd.DataFrame) -> AgGrid:
    """
    Plot dataframe in webpage
    themes: streamlit, balham, alpine, material
    """
    gb = GridOptionsBuilder.from_dataframe(data)
    gb.configure_pagination(paginationAutoPageSize=True)
    gb.configure_side_bar()
    gridOptions = gb.build()
    grid_response = AgGrid(
        data,
        gridOptions=gridOptions,
        data_return_mode="AS_INPUT",
        update_mode="MODEL_CHANGED",
        fit_columns_on_grid_load=False,
        theme="alpine",
        enable_enterprise_modules=True,
        height=350,
        width="100%",
        reload_data=True,
    )
    return grid_response


def set_example():
    example_data_dir = Path(Path(parent_dir.parent).parent) / "tests"
    st.session_state.sequence_data_uploaded = True
    st.session_state.search_state.prefix = "example_"
    st.session_state.search_state.data = example_data_dir / "test_data/MG1655.fasta"
    st.session_state.search_state.hmm_dir = example_data_dir / "test_data/hmms"
    st.session_state.search_state.hmm_meta = example_data_dir / "test_data/hmm_meta.tsv"
    search_outdir = Path(st.session_state.outdir) / "pynteny_example"
    if not st.session_state.outdir.exists():
        search_outdir.mkdir(parents=True, exist_ok=True)
    st.session_state.outdir = search_outdir
