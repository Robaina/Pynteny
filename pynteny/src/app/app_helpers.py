import pandas as pd
import streamlit as st
from st_aggrid import AgGrid, GridOptionsBuilder


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