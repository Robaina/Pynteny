#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tools to handle files in Streamlit App"""

import tkinter as tk
from pathlib import Path
from tkinter import filedialog

import streamlit as st

parent_dir = Path(__file__).parent


def show_files_in_dir(directory: Path, sidebar: bool = False) -> None:
    """
    Show a list of files in directory
    """
    filelist = [
        file.name
        for file in directory.iterdir()
        if (not file.name.startswith(".") and file.is_file())
    ]
    filelist.sort(key=lambda x: x.lower())
    dirlist = [
        file.name
        for file in directory.iterdir()
        if (not file.name.startswith(".") and file.is_dir())
    ]
    dirlist.sort(key=lambda x: x.lower())
    itemlist = filelist + dirlist
    iconlist = [":page_facing_up:" for _ in filelist] + [
        ":open_file_folder:" for _ in dirlist
    ]

    markdown_table = (
        f"\n| :diamond_shape_with_a_dot_inside: | {directory} |\n| --- | --- | "
    )
    for icon, object in zip(iconlist, itemlist):
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
    root.geometry("700x350")
    root.withdraw()
    try:
        selected_path = Path(filedialog.askopenfilename())
    except Exception:
        selected_path = None
    root.destroy()
    return selected_path


def open_directory_explorer() -> Path:
    """
    Open a explorer to select directory and return path
    """
    root = tk.Tk()
    root.geometry("700x350")
    root.withdraw()
    try:
        selected_dir = Path(filedialog.askdirectory(master=root))
    except Exception:
        selected_dir = None
    root.destroy()
    return selected_dir
