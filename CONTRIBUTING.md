# Contributing to Pynteny

First of all, thanks for taking the time to contribute! :tada::+1:

Here you will find a set of guidelines for contributing to Pynteny. Feel free to propose changes to this document in a pull request.

## Code of conduct

This project and everyone participating in it is governed by the [Contributor Covenant, v2.0](CODE_OF_CONDUCT.md) code of conduct. By participating, you are expected to uphold this code.

## I have a question!

If you only have a question about all things related to Pynteny, the best course of actions for you is to open a new [discussion](https://github.com/Robaina/Pynteny/discussions).

## How can I contribute?

### 1. Reporting bugs

We all make mistakes, and the developers behind Pynteny are no exception... So, if you find a bug in the source code, please open an [issue](https://github.com/Robaina/Pynteny/issues) and report it. Please, first search for similar issues that are currrently open.

### 2. Suggesting enhancements

Are you missing some feature that would like Pynteny to have? No problem! You can contribute by suggesting an enhancement, just open a new issue and tag it with the [```enhancement```](https://github.com/Robaina/Pynteny/labels/enhancement) label. Please, first search for similar issues that are currrently open.

### 3. Improving the documentation

Help is always needed at improving the [documentation](https://robaina.github.io/Pynteny/). Either adding more detailed docstrings, usage explanations or new examples.

## First contribution

Unsure where to begin contributing to Pynteny? You can start by looking for issues with the label [```good first issue```](https://github.com/Robaina/Pynteny/labels/good%20first%20issue). If you are unsure about how to set a developer environment for Pynteny, do take a look at the section below. Thanks!

## Setting up a local developer environment

Pynteny depends on packages that are not available in pip, namely [HMMER](https://github.com/EddyRivasLab/hmmer) and [Prodigal](https://github.com/hyattpd/Prodigal). These can be installed from the bioconda channel. Hence, to setup up a developer environment for Pynteny:

1. Fork and download repo, cd to downloaded directory. You should create a new branch to work on your issue.

2. Create conda environment with required dependencies:

The file `envs/pynteny-dev.yml` contains all dependencies required to use Pynteny. Conda is very slow solving the environment. It is recommended to use [mamba](https://github.com/mamba-org/mamba) instead:

```bash
mamba env create -n pynteny-dev -f envs/pynteny-dev.yml
conda activate pynteny-dev
```

3. Build package

```bash
(pynteny-dev) poetry build
```

4. Install Pynteny

```bash
(pynteny-dev) pip install dist/pynteny*.whl
```

5. Run tests

```bash
(pynteny-dev) python -m unittest discover tests
```

## Using GitHub codespaces

Alternatively, you can directly work on a developer environment in the browser within GitHub's codespaces. To set this environment up:

1. Fork repository and create new branch for your issue.

2. Start a new Codespace for the new branch. GitHub will build an environment with all required dependencies as well as pynteny the first time (it will take a couple of minutes). Pynteny will be installed in a conda environment named "pynteny-dev".

__Known issues__:
Port forwarding within the codespace doesn't work very well for Streamlit (pynteny app). This is a [known issue](https://discuss.streamlit.io/t/how-to-make-streamlit-run-on-codespaces/24526) regarding Streamlit and Codespaces.

## Building the documentation

The documentation is formed by a series of markdown files located in directory [docs](https://github.com/Robaina/Pynteny/tree/main/docs). This repo uses [mkdocs](https://www.mkdocs.org/) to automatically generate documentation pages from markdown files. Also, [MathJax](https://github.com/mathjax/MathJax) syntax is allowed!

This means that, to modify the [API reference](https://robaina.github.io/Pynteny/references/api/), all you need to do is to modify the docstring directly in the source file where the definion/class is located. And, to update the documentation pages, you just have to update the corresponding markdown file in the [docs](https://github.com/Robaina/Pynteny/tree/main/docs) directory. Note that, if you need to change the documentation structure (e.g., add or new pages),you would need to tell mkdocs about this change through its [configuration file](https://github.com/Robaina/Pynteny/blob/main/mkdocs.yml). Or just open an issue and ask for help!

When all the changes are ready to deploy, just open a pull request. After reviewing and merging the changes, the documentation will be automatically deployed.

Run the documentation locally with:

> mkdocs serve

## Tests on push and pull request to main

Pynteny's repo contains a [GitHub Action](https://github.com/features/actions) to perform build and integration tests which is triggered automatically on push and pull request events to the main brach. Currently the tests include building and installing Pynteny in Ubuntu and MacOS and running the [test](tests) suit.
