from importlib import metadata
from pathlib import Path

import streamlit as st
from PIL import Image

from pynteny.src.utils import CommandArgs
from pynteny.src.app.helpers import ExampleSearch
from pynteny.src.app.components import Sidebar, Mainpage


parent_dir = Path(Path(__file__).parent)
meta = metadata.metadata("pynteny")
__version__ = meta["Version"]
__author__ = meta["Author"]


APP_TITLE = "Pynteny — Synteny-aware HMM searches made easy"
icon = Image.open(parent_dir / "assets/img/logo.png")
favicon = Image.open(parent_dir / "assets/img/favicon.png")
st.session_state.sidebar_icon = icon

st.set_page_config(
    page_title=APP_TITLE,
    page_icon=favicon,
    layout="centered",
    initial_sidebar_state="auto",
    menu_items={
        'Get Help': None,
        'Report a bug': None,
        'About': (
            "## Pynteny — Synteny-aware HMM searches made easy\n"
            "A python tool to query sequence databases based on hmms "
            "arranged in a synteny block. For more info visit Pynteny's "
            "[repo](https://github.com/Robaina/Pynteny).\n"
            "### Citation\n"
            "If you use this software, please cite it as below:  \n"
            " Semidán Robaina Estévez. (2022). Pynteny: synteny-aware hmm searches made easy (Version 0.0.2). Zenodo. https://doi.org/10.5281/zenodo.7048685 \n"
            "### Using:\n"
            )
    }
)

with open(parent_dir / "assets/styles.css", "r") as file:
    css_text = file.read()
st.markdown(f"<style>{css_text}</style>", unsafe_allow_html=True)

with open(parent_dir / "assets/script.js", "r") as file:
    js_text = file.read()
st.components.v1.html(f"<script>{js_text}</script>")

# State variables
search_state = CommandArgs(
    data=None,
    synteny_struc=None,
    hmm_dir=None,
    hmm_meta=None,
    outdir=None,
    prefix="",
    processes=None,
    hmmsearch_args=None,
    gene_ids=True,
    logfile=None,
    synteny_hits=None,
    unordered=False
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
if "pynteny_log" not in st.session_state:
    st.session_state.pynteny_log = None


Sidebar.show()
Mainpage.show()
ExampleSearch.setExample()