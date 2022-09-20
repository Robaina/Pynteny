// Add custom JavaScript functionality to Streamlit

let footer  = window.parent.document.getElementsByTagName("footer")[0];
old_text = footer.innerHTML;
footer.innerHTML = "< > by <a href='https://github.com/Robaina'>Semidán Robaina</a>, 2022. " + old_text;
console.log(footer.innerHTML)