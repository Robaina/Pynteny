// Add custom JavaScript functionality to Streamlit

function setFooter() {
    let footer  = window.parent.document.getElementsByTagName("footer")[0];
    footer.innerHTML = "< > by <a href='https://github.com/Robaina'>Semidán Robaina</a>, 2022. " + "Made with <a href='https://streamlit.io/'>Streamlit</a>.";
}

setFooter();