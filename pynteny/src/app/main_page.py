import os
from importlib import metadata
from pathlib import Path

import streamlit as st
from PIL import Image

from pynteny.src.utils import CommandArgs
import pynteny.src.app.app_helpers as helpers
from pynteny.src.app.app_helpers import Callbacks, ExampleSearch


parent_dir = Path(Path(__file__).parent)
meta = metadata.metadata("pynteny")
__version__ = meta["Version"]
__author__ = meta["Author"]


APP_TITLE = "Pynteny — Synteny-aware HMM searches made easy"
icon = Image.open("assets/logo.png")
favicon = Image.open("assets/favicon.png")

st.set_page_config(
    page_title=APP_TITLE,
    page_icon=favicon,
    layout="centered",
    initial_sidebar_state="auto",
)

with open(parent_dir / "styles.css", "r") as file:
    css_text = file.read()
st.markdown(f"<style>{css_text}</style>", unsafe_allow_html=True)

# State variables
search_state = CommandArgs(
    data=None,
    synteny_struc=None,
    hmm_dir=None, #Path("/home/robaina/Documents/Pynteny/hmm_data/hmm_PGAP"),
    hmm_meta=None, #Path("/home/robaina/Documents/Pynteny/hmm_data/hmm_PGAP_no_missing.tsv"),
    outdir=None,
    prefix="",
    processes=None,
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
if "download_state" not in st.session_state:
    st.session_state.download_state = download_state
if "outdir" not in st.session_state:
    st.session_state.outdir = Path.cwd()
if "sequence_data_uploaded" not in st.session_state:
    st.session_state.sequence_data_uploaded = False

ExampleSearch.setExample()

# Side bar
st.sidebar.image(
    icon, use_column_width=True, caption=f"Pynteny v{__version__}"
)


st.sidebar.header("Search database")

st.sidebar.button("Select output directory", on_click=Callbacks.updateOutputDir)
col1, col2 = st.sidebar.columns([.5, .5])
with col1:
    st.text_input("",
                value="", max_chars=None,
                key="subdirectory", on_change=Callbacks.updateOutputSubdirectory,
                placeholder="Enter subdirectory")
with col2:
    st.text_input("",
                  value="", max_chars=None,
                  key="prefix", on_change=Callbacks.updateOutputPrefix,
                  placeholder="Enter output prefix")

if st.session_state.outdir is not None:
    files_div = st.empty()
    with files_div.container():
        helpers.show_files_in_dir(st.session_state.outdir, sidebar=True)

st.sidebar.text(" ")

with st.sidebar.expander("Advanced parameters:", expanded=False):

    st.markdown("Select custom HMM database")
    col1, col2 = st.columns([1, 1])
    with col1:
        st.button("HMM directory", on_click=Callbacks.selectHMMdir)
    with col2:
        st.button("HMM metadata", on_click=Callbacks.selectHMMmeta)
    
    col1, col2 = st.columns([0.4, 0.6])
    with col1:
        st.slider("Processes",
                  min_value=1, max_value=os.cpu_count(),
                  value=os.cpu_count() - 1, step=1,
                  on_change=Callbacks.setNumberOfProcesses, key="processes"
                  )



st.sidebar.button("Close session", key="close", on_click=Callbacks.close_session)


# Main page
st.title("Pynteny — Synteny-aware HMM searches made easy")
st.markdown("Welcome! This is a web app for the Pynteny package.")
st.markdown("Pynteny is a Python package for synteny-aware HMM searches.")


with st.expander("Search database by synteny-aware HMMs.", expanded=False):
    st.info(
        """
        Synteny blocks are specified by strings of ordered HMM names or gene IDs with the following format:\n
        $$\lt HMM_a \space n_{ab} \space \lt HMM_b \space n_{bc} \space \lt(HMM_{c1}|HMM_{c2}|HMM_{c3}),$$\n
        where $n_{ab}$ corresponds to the maximum number of genes between $HMM_a$ and $HMM_b$. Results can be 
        strand-specific, in that case $>$ preceding a HMM name indicates that the corresponding ORF must be
        located in the positive (or sense) strand. Likewise, a $<$ symbol indicates that the ORF must be located
        in the negative (antisense) strand. Searches can be made strand-insensitive by omitting the $>$ or $<$ symbol. 
        Several HMMs can be assigned to the same ORF, in which case the search is performed for all of them.
        In this case, HMM names must be separated by "|" and grouped within parentheses, as shown above.

        Sequence data can be either:
        - nucleotide assembly data in FASTA format or
        - a GenBank file containing sequence annotations.

        **Note: This Pynteny instance is run locally, thus files are always kept in your machine.
        """
    )


with st.expander("Select sequence data", expanded=True):
    col1, col2 = st.columns([1, 1])
    with col1:
        file_uploaded = st.button("Upload file", on_click=Callbacks.uploadData)
    with col2:
        database_builded = st.button("Build database", on_click=Callbacks.build)

if st.session_state.sequence_data_uploaded:
    st.info(f"Uploaded file: {st.session_state.search_state.data.name}")


with st.expander("Enter synteny structure", expanded=True):
    col1, col2 = st.columns([0.8, 0.2])
    with col1:
        st.session_state.search_state.synteny_struc = st.text_input("", "<leuD 0 <leuC 1 leuA")
    with col2:
        locus_repr = st.select_slider("", options=["Gene ID", "HMM"], value="Gene ID")
    gene_ids = True if locus_repr == "Gene ID" else False
    st.session_state.search_state.gene_ids = gene_ids


st.button("Search!", on_click=Callbacks.search)

plot_div = st.empty()
if st.session_state.search_state.synteny_hits is not None:
    with plot_div.container():
        helpers.plot_dataframe(st.session_state.search_state.synteny_hits)