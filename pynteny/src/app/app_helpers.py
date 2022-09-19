from pathlib import Path

import pandas as pd
import tkinter as tk
from tkinter import filedialog
import streamlit as st
from st_aggrid import AgGrid, GridOptionsBuilder

from pynteny.src.subcommands import synteny_search, build_database


def set_welcome_page():
    """
    Main landing page
    """
    with st.session_state.welcome_page.container():
        st.title("Pynteny — Synteny-aware HMM searches made easy")
        st.markdown("Welcome! This is a web app for the Pynteny package.")
        st.markdown("Pynteny is a Python package for synteny-aware HMM searches.")
        st.markdown("Please select a command from the sidebar to get started.")

def empty_welcome_page():
    """
    Remove welcome page
    """
    st.session_state.welcome_page.empty()

def close_session():
    st.markdown("Thanks for using Pynteny!")
    st.markdown("Please stop server by pressing control + c in terminal")
    st.stop()

def show_files_in_dir(directory: Path, sidebar: bool = True) -> None:
    """
    Show a list of files in directory

    :open_file_folder:
    :page_facing_up:
    """
    filelist = [file.name for file in directory.iterdir() if not file.name.startswith(".")]
    filelist.sort(key = lambda x: x.lower())
    fileicons = []
    for file in filelist:
        if Path(file).suffix:
            fileicons.append(":page_facing_up:")
        else:
            fileicons.append(":open_file_folder:")

    markdown_table = f"\n| :avocado: | {directory} |\n| --- | --- | "
    for icon, object in zip(fileicons, filelist):
        markdown_table += f"\n| {icon} | {object} | "
    markdown_table += "\n\n"

    if sidebar:
        st.sidebar.markdown(markdown_table)
    else:
        st.markdown(markdown_table)

def open_file_explorer() -> Path:
    """
    Open a file explorer and return selected path
    """
    root = tk.Tk()
    root.withdraw()
    try:
        selected_path = Path(filedialog.askopenfilename())
    except:
        selected_path = None
    root.destroy()
    return selected_path

def open_directory_explorer() -> Path:
    """
    Open a explorer to select directory and return path
    """
    root = tk.Tk()
    root.withdraw()
    try:
        selected_dir = Path(filedialog.askdirectory(master=root))
    except:
        selected_dir = None
    root.destroy()
    return selected_dir

def plot_dataframe(data: pd.DataFrame) -> AgGrid:
    """
    Plot dataframe in webpage
    """
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



class Callbacks:
    @staticmethod
    def search():
        st.session_state.search_state.outdir = st.session_state.outdir
        if st.session_state.search_state.data is not None and st.session_state.search_state.data.exists():
            synhits = synteny_search(st.session_state.search_state).getSyntenyHits()
            st.session_state.search_state.synteny_hits = synhits[[c for c in synhits.columns if c !="full_label"]]
            st.success("Search completed!")
        else:
            st.warning("Please, first upload a sequence database file")

    @staticmethod
    def build():
        if not st.session_state.sequence_data_uploaded:
            st.warning("Please, first upload assembly data file")
        else:
            st.session_state.build_state.data = st.session_state.search_state.data
            st.session_state.build_state.outdir = st.session_state.search_state.outdir
            st.session_state.build_state.outfile = Path(st.session_state.search_state.data.parent) / f"{st.session_state.search_state.data.stem}_labelled.faa"
            st.session_state.search_state.data = st.session_state.build_state.outfile
            build_database(st.session_state.build_state)

    @staticmethod
    def uploadData():
        selected_path = open_file_explorer()
        st.session_state.search_state.data = selected_path
        st.session_state.sequence_data_uploaded = True

    @staticmethod
    def updateOutdir():
        selected_outdir = open_directory_explorer()
        st.session_state.outdir = selected_outdir

    @staticmethod
    def close_session():
        st.markdown("Thanks for using Pynteny!")
        st.markdown("Please stop server by pressing control + c in terminal")
        st.stop()
