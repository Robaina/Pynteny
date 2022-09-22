/* 
Add custom JavaScript functionality to Streamlit
 Note: JS is injected within an iframe, thus outside elements must be accessed with "window.parent.document"
*/

function setFooter() {
    let footer  = window.parent.document.getElementsByTagName("footer")[0];
    footer.innerHTML = "< > by <a href='https://github.com/Robaina'>Semidán Robaina</a>, 2022. " + "Made with <a href='https://streamlit.io/'>Streamlit</a>.";
}

function addInfoClassToExpanders() {
    /*
    Was a good idea, but streamlit seems to remake the page upon any event,
    and assigned css classes are "forgotten"
    */
    let expanders = window.parent.document.getElementsByClassName("streamlit-expanderHeader");
    console.log(expanders);
    for (let i = 0; i < expanders.length; i++){
        let expander = expanders[i];
        console.log(i, expander);
        if (expander.ariaExpanded === "false") {
            expander.classList.add("pynteny-expander-info");
        }
    }
  }


setFooter();
// setTimeout(addInfoClassToExpanders, 1000);