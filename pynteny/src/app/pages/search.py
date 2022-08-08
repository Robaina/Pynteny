import streamlit as st
from PIL import Image






APP_TITLE = "Pynteny — Synteny-aware HMM searches made easy"
icon = Image.open("assets/pynteny_logo_2.png")


st.set_page_config(page_title="Search", page_icon="📈")

st.sidebar.header("Search database")
st.sidebar.markdown("Search...")

st.markdown("# Search")

st.markdown(
    """Search database by synteny-aware HMMs."""
)

# st.sidebar.image(
#     icon, use_column_width=True, caption="Pynteny v0.0.1"
# )


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