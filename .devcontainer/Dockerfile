FROM mcr.microsoft.com/devcontainers/miniconda:0-3

WORKDIR /tmp/conda-tmp/
# Copy repo to docker container
ADD ../ .
# Make conda environment and activate
RUN conda install mamba -n base -c conda-forge
RUN mamba env create -f envs/pynteny-dev.yml
SHELL ["conda", "run", "-n", "pynteny-dev", "/bin/bash", "-c"]
# Build and install Pynteny
RUN poetry build && pip install dist/pynteny*.whl && pynteny --help
# Give read/write permissions to install directory (needed to make config.json)
RUN chmod ugo+rw /opt/conda/envs/pynteny-dev/lib/python3.10/site-packages

# Initialize conda for default user: vscode
RUN runuser -l vscode -c 'conda init'
# Activate pynteny environment by default
RUN echo "conda activate pynteny-dev" >> /home/vscode/.bashrc
RUN source /home/vscode/.bashrc