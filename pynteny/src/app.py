from pathlib import Path

import streamlit as st
from st_aggrid import AgGrid, GridOptionsBuilder
from PIL import Image

import tkinter as tk  
from tkinter import filedialog  
import pandas as pd

from pynteny.src.utils import CommandArgs
from pynteny.src.subcommands import synteny_search



def plot_dataframe(data: pd.DataFrame, search_state: CommandArgs) -> None:
    # if search_state.synteny_hits is not None:
    gb = GridOptionsBuilder.from_dataframe(data)
    gb.configure_pagination(paginationAutoPageSize=True)
    gb.configure_side_bar()
    gridOptions = gb.build()
    grid_response = AgGrid(
        data,
        gridOptions=gridOptions,
        data_return_mode='AS_INPUT', 
        update_mode='MODEL_CHANGED', 
        fit_columns_on_grid_load=False,
        theme='blue',
        enable_enterprise_modules=True,
        height=350, 
        width='100%',
        reload_data=True
    )
    return grid_response


# Set the configs
APP_TITLE = "Pynteny — Synteny-aware HMM searches made easy"
icon = Image.open("assets/pynteny_logo_2.png")

st.set_page_config(
    page_title=APP_TITLE,
    page_icon=icon,
    layout="centered",
    initial_sidebar_state="auto",
)

st.title(APP_TITLE)


# Set the sidebar
st.sidebar.image(
    icon, use_column_width=True, caption="Pynteny v0.0.1"
)

st.sidebar.markdown(
    "# [Options](https://github.com/Robaina/Pynteny)"
)

st.sidebar.slider(
    "Set threads", min_value=1, max_value=16, value=1
)

st.sidebar.number_input(
    "Number of seqs", min_value=1, max_value=100, value=1
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


# def data_uploader(state):
#     """
#     Upload input data and return paths
#     """
#     return 

# Set main content
with st.expander("Select sequence data and synteny structure", expanded=False):
    st.info(
        """
        Sequence data can be either:
        - nucleotide assembly data in FASTA format or
        - a GenBank file containing sequence annotations.

        **Note: This Pynteny instance is run locally, thus files are always kept in your machine.
    """
    )

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
        """
    )

# file_buffer = st.file_uploader(
#     "Upload your dataset below", type=["fasta", "gbk", "gb"]
# )

file_uploaded = st.button("Upload file")
if file_uploaded:
    root = tk.Tk()  
    root.withdraw()  
    search_state.data = Path(filedialog.askopenfilename())
    st.info(f"Uploaded file: {search_state.data.name}")

search_state.synteny_struc = st.text_input("Enter synteny structure", "<leuD 0 <leuC 1 leuA")


def search():
    if search_state.data is not None and search_state.data.exists():
        synhits = synteny_search(search_state).getSyntenyHits()
        search_state.synteny_hits = synhits[[c for c in synhits.columns if c !="full_label"]]
        st.success("Search completed!")
    else:
        st.warning("Please upload a file")

    if search_state.synteny_hits is not None:
        plot_dataframe(search_state.synteny_hits)



st.button("Search!", on_click=search)




# grid_response = plot_dataframe(search_state.synteny_hits, search_state)

# if search_state.synteny_hits is not None:
#     # AgGrid(search_state.synteny_hits)
#     df = pd.read_csv("/home/robaina/Documents/Pynteny/test_results/synteny_matched.tsv", sep="\t")
#     plot_dataframe(df) #search_state.synteny_hits)

# synteny_hits = Path("/home/robaina/Documents/Pynteny/test_results") / "synteny_matched.tsv"
# if synteny_hits.exists():
#     df = pd.read_csv(synteny_hits, sep="\t")
#     AgGrid(df)

# if search_state.synteny_hits is not None:
#     AgGrid(search_state.synteny_hits)

# if file_buffer is not None:
#     st.info(os.path.abspath(file_buffer.name))