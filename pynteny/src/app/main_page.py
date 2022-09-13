from distutils.command.build import build
from pathlib import Path

import streamlit as st
from PIL import Image
import tkinter as tk
from tkinter import filedialog

from pynteny.src.utils import CommandArgs
from pynteny.src.subcommands import synteny_search, build_database
import pynteny.src.app.app_helpers as helpers



def close_session():
    st.markdown("Thanks for using Pynteny!")
    st.markdown("Please stop server by pressing control + c in terminal")
    st.stop()


APP_TITLE = "Pynteny — Synteny-aware HMM searches made easy"
icon = Image.open("assets/logo.png")
favicon = Image.open("assets/favicon.png")

st.set_page_config(
    page_title=APP_TITLE,
    page_icon=favicon,
    layout="centered",
    initial_sidebar_state="auto",
)


# Side bar
st.sidebar.image(
    icon, use_column_width=True, caption="Pynteny v0.0.1"
)

st.sidebar.success("Select a demo above.")

st.sidebar.button("Close session", key="close", on_click=close_session)

st.sidebar.header("Search database")
st.sidebar.markdown("Search...")

st.sidebar.info(
    """
    Sequence data can be either:
    - nucleotide assembly data in FASTA format or
    - a GenBank file containing sequence annotations.

    **Note: This Pynteny instance is run locally, thus files are always kept in your machine.
"""
)

st.sidebar.info(
    """
    Synteny blocks are specified by strings of ordered HMM names or gene IDs with the following format:\n
    $$\lt HMM_a \space n_{ab} \space \lt HMM_b \space n_{bc} \space \lt(HMM_{c1}|HMM_{c2}|HMM_{c3}),$$\n
    where $n_{ab}$ corresponds to the maximum number of genes between $HMM_a$ and $HMM_b$. Results can be 
    strand-specific, in that case $>$ preceding a HMM name indicates that the corresponding ORF must be
    located in the positive (or sense) strand. Likewise, a $<$ symbol indicates that the ORF must be located
    in the negative (antisense) strand. Searches can be made strand-insensitive by omitting the $>$ or $<$ symbol. 
    Several HMMs can be assigned to the same ORF, in which case the search is performed for all of them.
    In this case, HMM names must be separated by "|" and grouped within parentheses, as shown above.
    """
)


# Main page
st.title("Pynteny — Synteny-aware HMM searches made easy")
st.markdown("Welcome! This is a web app for the Pynteny package.")
st.markdown("Pynteny is a Python package for synteny-aware HMM searches.")
st.markdown("Please select a command from the sidebar to get started.")


st.markdown("# Search")
st.markdown(
    """Search database by synteny-aware HMMs."""
)

# State variables
search_state = CommandArgs(
    data=None,
    synteny_struc=None,
    hmm_dir=Path("/home/robaina/Documents/Pynteny/hmm_data/hmm_PGAP"),
    hmm_meta=Path("/home/robaina/Documents/Pynteny/hmm_data/hmm_PGAP_no_missing.tsv"),
    outdir=Path("/home/robaina/Documents/Pynteny/test_results"),
    prefix="",
    hmmsearch_args=None,
    gene_ids=True,
    logfile=None,
    synteny_hits=None
    )

build_state = CommandArgs(
    data=None,
    outfile=None,
    processes=None,
    logfile=None
)

download_state = CommandArgs(
    outdir=None,
    unpack=True,
    logfile=None
)

if "search_state" not in st.session_state:
    st.session_state.search_state = search_state
if "build_state" not in st.session_state:
    st.session_state.build_state = build_state

def search():
    if st.session_state.search_state.data is not None and st.session_state.search_state.data.exists():
        synhits = synteny_search(st.session_state.search_state).getSyntenyHits()
        st.session_state.search_state.synteny_hits = synhits[[c for c in synhits.columns if c !="full_label"]]
        st.success("Search completed!")
    else:
        st.warning("Please, first upload a sequence database file")


with st.expander("Select sequence data", expanded=True):
    col1, col2 = st.columns([1, 1])
    with col1:
        file_uploaded = st.button("Upload file")
    with col2:
        database_builded = st.button("Build database")

if file_uploaded:
    root = tk.Tk()
    root.withdraw()
    st.session_state.search_state.data = Path(filedialog.askopenfilename())
    st.info(f"Uploaded file: {st.session_state.search_state.data.name}")
    root.destroy()

if database_builded:
    if not file_uploaded:
        st.warning("Please, first upload assembly data file")
    else:
        build_database(st.session_state.build_state)


with st.expander("Enter synteny structure", expanded=True):
    st.session_state.search_state.synteny_struc = st.text_input("", "<leuD 0 <leuC 1 leuA")


st.button("Search!", on_click=search)

plot_div = st.empty()
if st.session_state.search_state.synteny_hits is not None:
    with plot_div.container():
        helpers.plot_dataframe(st.session_state.search_state.synteny_hits)