import streamlit as st
from PIL import Image


def close_session():
    st.markdown("Thanks for using Pynteny!")
    st.markdown("Please stop server by pressing control + c in terminal")
    st.stop()


APP_TITLE = "Pynteny — Synteny-aware HMM searches made easy"
icon = Image.open("assets/pynteny_logo_2.png")

st.set_page_config(
    page_title=APP_TITLE,
    page_icon=icon,
    layout="centered",
    initial_sidebar_state="auto",
)



st.sidebar.image(
    icon, use_column_width=True, caption="Pynteny v0.0.1"
)

st.sidebar.success("Select a demo above.")

st.sidebar.button("Close session", key="close", on_click=close_session)


st.title("Pynteny — Synteny-aware HMM searches made easy")
st.markdown("Welcome! This is a web app for the Pynteny package.")
st.markdown("Pynteny is a Python package for synteny-aware HMM searches.")
st.markdown("Please select a command from the sidebar to get started.")