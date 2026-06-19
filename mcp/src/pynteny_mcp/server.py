"""
Pynteny MCP server.

Exposes the Pynteny API as Model Context Protocol tools so that any
MCP-compatible agent (Claude Desktop, Claude Code, or a custom client driving
Claude / DeepSeek / etc.) can run synteny-aware HMM searches over genomic
sequence data in natural language.

Pynteny (https://github.com/Robaina/Pynteny) finds synteny blocks in
(prokaryotic) sequence data: it searches a labelled peptide database with
profile HMMs (via HMMER / pyhmmer) and keeps only the hits whose genomic
arrangement — order, strand and gene spacing — matches a user-supplied synteny
structure such as ``>leuD 0 >leuC 1 <leuA``.

Run it (stdio transport, the default):

    pynteny-mcp
    # or, without installing this package:
    PYTHONPATH=src python -m pynteny_mcp.server

Pynteny is a pure-Python pip package (HMMER and Prodigal come bundled via
pyhmmer / pyrodigal), so no conda environment or external binaries are needed —
only an environment in which `pip install pynteny` has been run.
"""

from __future__ import annotations

from typing import Optional

from mcp.server.fastmcp import FastMCP

from . import service

INSTRUCTIONS = """\
Tools to run Pynteny synteny-aware HMM searches over genomic sequence data.

A *synteny structure* describes a target gene arrangement, e.g.
  ">leuD 0 >leuC 1 <leuA"
where each token is an HMM name (or gene symbol, with gene_ids=True), the
leading ">"/"<" is the strand (sense / antisense), and the integers are the
maximum number of (untargeted) genes allowed between neighbours. Groups of
interchangeable HMMs for one gene are written "(HMM_A|HMM_B)".

Typical workflow:
  1. `get_pynteny_info` — confirm the install and see which HMM databases
     (PGAP / PFAM) are already downloaded.
  2. `validate_synteny_structure` — cheaply check/parse a structure before a run.
  3. (once) `download_hmm_databases` — fetch the PGAP HMMs + metadata, unless you
     already have an `hmm_dir`.
  4. `build_peptide_database` — turn a nucleotide assembly / GenBank file into the
     labelled peptide FASTA that search consumes (skip if you already have one).
  5. `run_synteny_search` — the core tool: returns the matched hits table and the
     paths of the FASTA / TSV files Pynteny writes.

Use `parse_gene_symbols` to translate a gene-symbol structure into HMM names
against a metadata file. Searching by gene symbol directly is done by passing
gene_ids=True to `run_synteny_search`.
"""

mcp = FastMCP("pynteny", instructions=INSTRUCTIONS)


@mcp.tool()
def get_pynteny_info() -> dict:
    """Return the installed Pynteny version, author, citation, and which HMM
    databases (PGAP / PFAM) are currently downloaded and configured. A good
    first call to confirm the server is wired up and to decide whether you need
    to download HMMs before searching."""
    return service.pynteny_info()


@mcp.tool()
def validate_synteny_structure(synteny_structure: str) -> dict:
    """Validate and decompose a synteny structure string without running any
    search or touching any file. Returns whether it is well-formed plus its
    parsed parts: the HMM groups, the flat list of HMM names, the strand of each
    gene, and the max-distance constraints between neighbours. Use this to catch
    formatting mistakes before an expensive `run_synteny_search`.

    Args:
        synteny_structure: e.g. ">leuD 0 >leuC 1 <leuA". Tokens are HMM names (or
            gene symbols), ">"/"<" give the strand, integers are the maximum
            number of genes allowed between neighbours, and "(A|B)" groups
            interchangeable HMMs for one gene.
    """
    return service.validate_structure(synteny_structure)


@mcp.tool()
def parse_gene_symbols(synteny_structure: str, hmm_meta: Optional[str] = None) -> dict:
    """Translate a synteny structure written with *gene symbols* (e.g. "leuD")
    into one written with the corresponding *HMM names* (e.g. "TIGR00171.1"),
    using a PGAP/PFAM metadata table.

    Args:
        synteny_structure: gene-symbol structure, e.g. ">leuD 0 >leuC 1 <leuA".
        hmm_meta: path to the HMM metadata TSV. If omitted, the PGAP metadata
            file recorded in Pynteny's config (from a previous download) is used.
    """
    return service.parse_gene_symbols(synteny_structure, hmm_meta)


@mcp.tool()
def build_peptide_database(
    data: str,
    outfile: Optional[str] = None,
    prepend: bool = False,
    processes: Optional[int] = None,
    tempdir: Optional[str] = None,
    logfile: Optional[str] = None,
) -> dict:
    """Build a labelled peptide database from a nucleotide assembly (or a GenBank
    file / directory). ORFs are predicted (Prodigal/pyrodigal) and each peptide
    is labelled with its contig and gene-location info so that synteny search can
    reason about gene order and strand. Returns the output file path and the
    number of peptides written.

    Args:
        data: path to the input assembly FASTA, GenBank file, or a directory of
            such files.
        outfile: path for the labelled peptide FASTA. Defaults to a name derived
            from `data`.
        prepend: when `data` is a directory, prepend each file's name to its
            sequence IDs (keeps IDs unique across genomes).
        processes: max worker processes (defaults to all CPUs minus one).
        tempdir: directory for temporary files.
        logfile: path to a log file (logs otherwise go to stderr/devnull).
    """
    return service.build_database(
        data,
        outfile=outfile,
        prepend=prepend,
        processes=processes,
        tempdir=tempdir,
        logfile=logfile,
    )


@mcp.tool()
def run_synteny_search(
    data: str,
    synteny_structure: str,
    hmm_dir: Optional[str] = None,
    hmm_meta: Optional[str] = None,
    gene_ids: bool = False,
    unordered: bool = False,
    best_hmm_wins: bool = False,
    reuse: bool = False,
    outdir: Optional[str] = None,
    prefix: str = "",
    hmmsearch_args: Optional[str] = None,
    processes: Optional[int] = None,
    logfile: Optional[str] = None,
    max_hits: int = 200,
) -> dict:
    """Search a labelled peptide database for genes arranged in the given synteny
    structure. This is the core Pynteny operation. It runs HMMER for each HMM,
    then keeps only hits whose genomic context (order, strand, spacing) matches
    the structure. Returns the matched hits as records, plus the paths of the
    per-gene FASTA files and the `synteny_matched.tsv` table Pynteny writes.

    Args:
        data: path to the labelled peptide database (output of
            `build_peptide_database`).
        synteny_structure: target arrangement, e.g. ">leuD 0 >leuC 1 <leuA"
            (HMM names), or gene symbols if `gene_ids=True`.
        hmm_dir: directory containing the HMM files. If omitted, the downloaded
            PGAP directory from Pynteny's config is used; if none is available
            the call returns an error rather than triggering a large download.
        hmm_meta: HMM metadata TSV; needed for `gene_ids=True` and used to
            annotate hits with gene symbol / product / EC number when available.
        gene_ids: treat tokens in `synteny_structure` as gene symbols instead of
            HMM names (translated via `hmm_meta`).
        unordered: match a syntenic *set* in any order instead of the exact
            collinear order given.
        best_hmm_wins: when one peptide is hit by several HMMs (paralog
            cross-hits), keep only the highest-scoring HMM for that peptide.
        reuse: reuse existing HMMER outputs in `outdir` instead of re-running.
        outdir: output directory (defaults to the directory of `data`).
        prefix: prefix for output file names.
        hmmsearch_args: extra hmmsearch arguments (a single string applied to all
            HMMs, or comma-separated per-HMM with "None" to skip one).
        processes: max worker processes (defaults to all CPUs minus one).
        logfile: path to a log file.
        max_hits: cap on the number of hit records returned inline (the full set
            is always written to the TSV).
    """
    return service.run_search(
        data,
        synteny_structure,
        gene_ids=gene_ids,
        unordered=unordered,
        best_hmm_wins=best_hmm_wins,
        reuse=reuse,
        hmm_dir=hmm_dir,
        hmm_meta=hmm_meta,
        outdir=outdir,
        prefix=prefix,
        hmmsearch_args=hmmsearch_args,
        processes=processes,
        logfile=logfile,
        max_hits=max_hits,
    )


@mcp.tool()
def download_hmm_databases(
    outdir: str,
    pgap: bool = True,
    pfam: bool = False,
    unpack: bool = True,
    force: bool = False,
    logfile: Optional[str] = None,
) -> dict:
    """Download profile-HMM databases (PGAP from NCBI and/or PFAM-A) and register
    them in Pynteny's config so later searches can find them automatically. This
    downloads large files over the network and can take a while.

    Args:
        outdir: directory to download the HMM databases into.
        pgap: download the PGAP database (HMMs + metadata). Default True.
        pfam: download the PFAM-A database. Default False.
        unpack: unpack the archives now (otherwise Pynteny unpacks per session).
        force: re-download even if already present.
        logfile: path to a log file.
    """
    return service.download_databases(
        outdir,
        pgap=pgap,
        pfam=pfam,
        unpack=unpack,
        force=force,
        logfile=logfile,
    )


def main() -> None:
    """Console-script entry point: run the server over stdio."""
    import logging
    import os

    # Quiet the per-request INFO chatter from the MCP runtime unless debugging.
    if not os.environ.get("PYNTENY_MCP_DEBUG"):
        logging.getLogger("mcp").setLevel(logging.WARNING)
    mcp.run()


if __name__ == "__main__":
    main()
