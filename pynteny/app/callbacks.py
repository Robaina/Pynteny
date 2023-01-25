#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Callback definitions to be used by streamlit app"""

import shutil
from pathlib import Path
import streamlit as st

import pynteny.app.filemanager as filemanager
from pynteny.subcommands import synteny_search, build_database, download_hmms
from pynteny.utils import ConfigParser


parent_dir = Path(__file__).parent


def select_log_path():
    if st.session_state.log == "Yes":
        logfile = Path(st.session_state.outdir) / "pynteny.log"
    else:
        logfile = None
    st.session_state.pynteny_log = logfile


def update_log():
    """
    Update Pynteny log with Streamlit log info
    """
    config = ConfigParser.get_default_config()
    streamlit_log = config.get_field("streamlit_log")
    select_log_path()
    if st.session_state.pynteny_log is not None:
        shutil.copy(streamlit_log, st.session_state.pynteny_log)


def search():
    config = ConfigParser.get_default_config()
    if (st.session_state.search_state.hmm_dir is None) and (
        not config.get_field("data_downloaded")
    ):
        with st.spinner("Downloading HMM database, please wait..."):
            download_hmms(st.session_state.download_state)
        st.success("HMM database downloaded!")
    st.session_state.search_state.outdir = st.session_state.outdir
    if (
        st.session_state.search_state.data is not None
        and st.session_state.search_state.data.exists()
    ):
        synhits = synteny_search(st.session_state.search_state).getSyntenyHits()
        st.session_state.search_state.synteny_hits = synhits[
            [c for c in synhits.columns if c != "full_label"]
        ]
        update_log()
    else:
        st.warning("Please, first upload a sequence database file")


def build():
    if not st.session_state.sequence_data_uploaded:
        st.warning("Please, first upload assembly data file")
    else:
        st.session_state.build_state.data = st.session_state.search_state.data
        st.session_state.build_state.outdir = st.session_state.search_state.outdir
        st.session_state.build_state.outfile = (
            Path(st.session_state.search_state.data.parent)
            / f"{st.session_state.search_state.data.stem}_labelled.faa"
        )
        st.session_state.search_state.data = st.session_state.build_state.outfile
        build_database(st.session_state.build_state)
        update_log()


def upload_data():
    selected_path = filemanager.open_file_explorer()
    st.session_state.search_state.data = selected_path
    st.session_state.sequence_data_uploaded = True


def update_output_dir():
    selected_outdir = filemanager.open_directory_explorer()
    if selected_outdir is not None:
        st.session_state.outdir = selected_outdir


def update_output_prefix():
    st.session_state.search_state.prefix = st.session_state.prefix


def update_output_subdirectory():
    subdir = Path(st.session_state.outdir) / st.session_state.subdirectory
    if not subdir.exists():
        subdir.mkdir(parents=True, exist_ok=False)
    st.session_state.outdir = subdir


def select_HMM_dir():
    selected_dir = filemanager.open_directory_explorer()
    st.session_state.search_state.hmm_dir = selected_dir
    st.success(f"Selected HMM database: {selected_dir}")


def select_HMM_meta():
    selected_file = filemanager.open_file_explorer()
    st.session_state.search_state.hmm_meta = selected_file
    st.success(f"Selected HMM metadata: {selected_file}")


def set_number_of_processes():
    st.session_state.build_state.processes = st.session_state.processes
    st.session_state.search_state.processes = st.session_state.processes


def close_session():
    st.text(" ")
    st.text(" ")
    st.markdown("# Thanks for using")
    st.text(" ")
    st.image(st.session_state.sidebar_icon)
    st.text(" ")
    st.markdown("### Please stop the server by pressing control + c in terminal")
    st.stop()
