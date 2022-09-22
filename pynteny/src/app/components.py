import os
from importlib import metadata

import streamlit as st
from PIL import Image

from pynteny.src.app.helpers import Callbacks, FileManager, Plotter


meta = metadata.metadata("pynteny")
__version__ = meta["Version"]
__author__ = meta["Author"]
icon = Image.open("assets/logo.png")


class Sidebar:
    @staticmethod
    def show():

        st.sidebar.image(
            icon, use_column_width=True, caption=f"Synteny-aware HMM searches made easy (v{__version__})"
        )

        st.sidebar.text(" ")
        st.sidebar.text(" ")

        with st.sidebar.expander("Current output directory:", expanded=True):
            if st.session_state.outdir is not None:
                files_div = st.empty()
                with files_div.container():
                    FileManager.show_files_in_dir(st.session_state.outdir, sidebar=False)
            
            st.text(" ")
            st.button("Select directory", on_click=Callbacks.updateOutputDir)
            col1, col2 = st.columns([.5, .5])
            with col1:
                st.text_input("",
                            value="", max_chars=None,
                            key="subdirectory", on_change=Callbacks.updateOutputSubdirectory,
                            placeholder="Create subdirectory")
            with col2:
                st.text_input("",
                            value="", max_chars=None,
                            key="prefix", on_change=Callbacks.updateOutputPrefix,
                            placeholder="Enter output prefix")

        with st.sidebar.expander("Advanced parameters:", expanded=False):

            st.markdown("Select custom HMM database")
            col1, col2 = st.columns([1, 1])
            with col1:
                st.button("HMM directory", on_click=Callbacks.selectHMMdir)
            with col2:
                st.button("HMM metadata", on_click=Callbacks.selectHMMmeta)
            
            col1, col2 = st.columns([0.4, 0.6])
            with col1:
                st.markdown(("Processes:"))
                st.slider("Processes",
                        min_value=1, max_value=os.cpu_count(),
                        value=os.cpu_count() - 1, step=1,
                        on_change=Callbacks.setNumberOfProcesses, key="processes"
                        )
        
        col1, col2 = st.sidebar.columns([1, 1])
        with col2:
            st.button("Close session", key="close", on_click=Callbacks.close_session)


class Mainpage:
    @staticmethod
    def show():

        st.title("Pynteny — Synteny-aware HMM searches made easy")
        # with st.expander("Welcome!", expanded=False):
        #     st.info(
        #         """
        #         __Note__: This Pynteny instance is run locally, thus files are always kept in your machine.
        #         """
        #     )

        with st.expander("Select sequence data", expanded=False):
            st.info(
                """
                Sequence data can be either:

                - nucleotide assembly data in FASTA format or
                - a GenBank file containing sequence annotations.

                __Note__: This Pynteny instance is run locally, thus files are always kept in your machine.
                """
            )

        with st.expander("", expanded=True):
            col1, col2 = st.columns([1, 1])
            with col1:
                st.button("Upload file", on_click=Callbacks.uploadData)
            with col2:
                st.button("Build database", on_click=Callbacks.build)

        if st.session_state.sequence_data_uploaded:
            st.info(f"Uploaded file: {st.session_state.search_state.data.name}")

        with st.expander("Enter synteny structure:", expanded=False):
            st.info(
                """
                Synteny blocks are specified by strings of ordered HMM names or gene IDs with the following format: 

                $$\lt HMM_a \space n_{ab} \space \lt HMM_b \space n_{bc} \space \lt(HMM_{c1}|HMM_{c2}|HMM_{c3}),$$ 

                where $n_{ab}$ corresponds to the maximum number of genes between $HMM_a$ and $HMM_b$. Additionally:

                - Results can be strand-specific, in that case $>$ preceding a HMM name indicates that the corresponding ORF must be located in the positive (or sense) strand. Likewise, a $<$ symbol indicates that the ORF must be located in the negative (antisense) strand. 

                - Searches can be made strand-insensitive by omitting the $>$ or $<$ symbol. 

                - Several HMMs can be assigned to the same ORF, in which case the search is performed for all of them. In this case, HMM names must be separated by "|" and grouped within parentheses, as shown above.
                """
            )

        with st.expander("", expanded=True):
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
                Plotter.plot_dataframe(st.session_state.search_state.synteny_hits)