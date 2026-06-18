FROM python:3.11-slim

WORKDIR /Pynteny
# Copy repo to docker container
COPY src/pynteny src/pynteny/
COPY tests tests/
COPY README.md .
COPY pyproject.toml .
COPY LICENSE .

# Build and install Pynteny (pure pip, no conda required).
# pyhmmer and pyrodigal ship manylinux wheels, so no compiler toolchain is needed.
RUN pip install --no-cache-dir poetry \
    && poetry build \
    && pip install --no-cache-dir dist/pynteny*.whl \
    && pynteny --help
