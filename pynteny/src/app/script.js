/* 
Add custom JavaScript functionality to Streamlit
 Note: JS is injected within an iframe, thus outside elements must be accessed with "window.parent.document"
*/

function setFooter() {
    let footer  = window.parent.document.getElementsByTagName("footer")[0];
    footer.innerHTML = "< > by <a href='https://github.com/Robaina'>Semidán Robaina</a>, 2022. " + "Made with <a href='https://streamlit.io/'>Streamlit</a>.";
}

function addEventToRestartButton() {
    let buttons = window.parent.document.getElementsByTagName("button");
    for (let i = 0; i < buttons.length; i++){
        let button = buttons[i];
        if (button.innerHTML.includes("Restart")) {
            button.addEventListener("click", function() {
                window.parent.location.reload();}
                )
        }
    }
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
setTimeout(addEventToRestartButton, 1000);