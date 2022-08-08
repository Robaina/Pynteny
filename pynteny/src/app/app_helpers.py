import streamlit as st


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


def search_gui():
    """
    GUI for search command
    """
    with st.expander("Select sequence data and synteny structure", expanded=True):
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


def build_gui():
    """
    GUI for build command
    """
    pass

def download_gui():
    """
    GUI for download command
    """
    pass

def parse_gui():
    """
    GUI for parse command
    """
    pass