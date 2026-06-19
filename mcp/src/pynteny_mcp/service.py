"""
Service layer for the Pynteny MCP server.

Wraps the :mod:`pynteny` public API (``pynteny.api.Search / Build / Download``
and the synteny-structure parsers) with helpers that:

* keep the MCP **stdio channel clean** — Pynteny logs to ``stdout`` and several
  helpers print there, which would corrupt the JSON-RPC stream, so every call
  into Pynteny is wrapped in :func:`_pynteny_run` which redirects ``stdout`` to
  ``stderr`` and isolates Pynteny's root-logger handlers, and
* turn the rich return values (``SyntenyHits`` / pandas ``DataFrame`` / written
  output files) into compact, JSON-serialisable summaries.

Everything here is transport-agnostic; ``server.py`` only adds the MCP tool
definitions on top.
"""

from __future__ import annotations

import contextlib
import logging
import os
import sys
from pathlib import Path
from typing import Any, Optional


def _log(msg: str) -> None:
    # NEVER write to stdout: it is the MCP stdio channel. Logs go to stderr.
    print(f"[pynteny-mcp] {msg}", file=sys.stderr, flush=True)


@contextlib.contextmanager
def _pynteny_run():
    """
    Run a block of Pynteny code without polluting the MCP stdio channel.

    Pynteny's ``init_logger`` calls ``logging.basicConfig`` with a
    ``StreamHandler(sys.stdout)`` and several code paths ``print`` to stdout. The
    MCP stdio transport uses stdout for protocol traffic, so any of that would
    break the connection. We:

    * redirect ``sys.stdout`` to ``sys.stderr`` for the duration of the call
      (so the StreamHandler Pynteny installs binds to stderr, and stray prints
      go to stderr too), and
    * detach the root logger's existing handlers first and restore them after,
      so repeated calls don't accumulate handlers or reuse ones Pynteny closed
      via ``logging.shutdown()``.
    """
    root = logging.getLogger()
    saved_handlers = root.handlers[:]
    for h in saved_handlers:
        root.removeHandler(h)
    try:
        with contextlib.redirect_stdout(sys.stderr):
            yield
    finally:
        # Drop whatever handlers Pynteny installed (it may have closed them in
        # logging.shutdown()), then restore the ones we started with.
        for h in root.handlers[:]:
            root.removeHandler(h)
        for h in saved_handlers:
            root.addHandler(h)


# --------------------------------------------------------------------------- #
# Config / metadata                                                           #
# --------------------------------------------------------------------------- #
def pynteny_info() -> dict[str, Any]:
    """Version, author, citation and which HMM databases are configured."""
    from pynteny.api import Command, __author__, __version__
    from pynteny.utils import CommandArgs, ConfigParser
    from pynteny.subcommands import get_citation

    citation = get_citation(CommandArgs(version=__version__, author=__author__), silent=True)

    config = ConfigParser.get_default_config()
    db: dict[str, Any] = {}
    for field in (
        "PGAP_data_downloaded",
        "PFAM_data_downloaded",
        "PGAP_database",
        "PGAP_meta_file",
        "PFAM_database",
        "PFAM_meta_file",
        "database_dir",
    ):
        try:
            db[field] = config.get_field(field)
        except Exception:
            db[field] = None
    return {
        "version": __version__,
        "author": __author__,
        "citation": citation,
        "config_file": str(config.get_config_path()),
        "hmm_databases": db,
    }


def _resolve_hmm_meta(hmm_meta: Optional[str]) -> Optional[str]:
    """Fall back to the PGAP metadata file recorded in Pynteny's config."""
    if hmm_meta:
        return hmm_meta
    from pynteny.utils import ConfigParser

    config = ConfigParser.get_default_config()
    for field in ("PGAP_meta_file", "PFAM_meta_file"):
        try:
            value = config.get_field(field)
        except Exception:
            value = None
        if value:
            return value
    return None


def _resolve_hmm_dir(hmm_dir: Optional[str]) -> Optional[str]:
    """Fall back to the PGAP HMM directory recorded in Pynteny's config."""
    if hmm_dir:
        return hmm_dir
    from pynteny.utils import ConfigParser

    config = ConfigParser.get_default_config()
    for field in ("PGAP_database", "PFAM_database"):
        try:
            value = config.get_field(field)
        except Exception:
            value = None
        if value:
            return value
    return None


# --------------------------------------------------------------------------- #
# Synteny-structure parsing (no I/O)                                          #
# --------------------------------------------------------------------------- #
def validate_structure(synteny_structure: str) -> dict[str, Any]:
    """Reformat, validate and decompose a synteny structure string.

    Pure string work — no HMMs, files or searches involved. Several Pynteny
    parser helpers call ``sys.exit`` on malformed input, so we trap
    ``SystemExit`` and report ``valid=False`` instead of killing the server.
    """
    import pynteny.parsers.syntenyparser as sp

    reformatted = sp.reformat_synteny_structure(synteny_structure)
    result: dict[str, Any] = {
        "input": synteny_structure,
        "reformatted": reformatted,
    }
    try:
        with _pynteny_run():
            valid = sp.is_valid_structure(reformatted)
            hmm_groups = sp.get_HMM_groups_in_structure(reformatted)
            hmm_names = sp.get_all_HMMs_in_structure(reformatted)
            strands = sp.get_strands_in_structure(reformatted)
            distances = sp.get_maximum_distances_in_structure(reformatted)
    except SystemExit:
        result["valid"] = False
        result["error"] = (
            "Could not parse the structure. Expected something like "
            "'>hmm_a 0 >hmm_b 3 <hmm_c' (strand '>'/'<', integers = max genes "
            "between neighbours)."
        )
        return result
    except Exception as exc:  # noqa: BLE001
        result["valid"] = False
        result["error"] = f"{type(exc).__name__}: {exc}"
        return result

    result.update(
        {
            "valid": bool(valid),
            "n_genes": len(hmm_groups),
            "hmm_groups": hmm_groups,
            "hmm_names": hmm_names,
            "strands": strands,
            "max_distances": distances,
            "contains_hmm_groups": sp.contains_HMM_groups(reformatted),
        }
    )
    return result


def parse_gene_symbols(synteny_structure: str, hmm_meta: Optional[str]) -> dict[str, Any]:
    """Translate a *gene-symbol* synteny structure into one based on HMM names,
    using a PGAP/PFAM metadata table."""
    from pynteny.api import Search

    resolved_meta = _resolve_hmm_meta(hmm_meta)
    if not resolved_meta:
        return {
            "error": (
                "No HMM metadata file available. Provide `hmm_meta`, or download "
                "the PGAP database first with `download_hmm_databases`."
            )
        }
    if not Path(resolved_meta).exists():
        return {"error": f"HMM metadata file not found: {resolved_meta}"}

    search = Search(data=".", synteny_struc=synteny_structure, hmm_meta=resolved_meta)
    try:
        with _pynteny_run():
            translated = search.parse_genes(synteny_structure)
    except SystemExit:
        return {
            "input": synteny_structure,
            "hmm_meta": resolved_meta,
            "error": (
                "One or more gene symbols did not match an HMM in the metadata "
                "table. Check the gene symbols and the metadata file."
            ),
        }
    return {
        "input": synteny_structure,
        "hmm_meta": resolved_meta,
        "translated_structure": translated,
    }


# --------------------------------------------------------------------------- #
# Build                                                                       #
# --------------------------------------------------------------------------- #
def _count_fasta_records(path: Path) -> int:
    n = 0
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith(">"):
                    n += 1
    except OSError:
        return -1
    return n


def build_database(
    data: str,
    *,
    outfile: Optional[str],
    prepend: bool,
    processes: Optional[int],
    tempdir: Optional[str],
    logfile: Optional[str],
) -> dict[str, Any]:
    """Translate a nucleotide assembly (or GenBank file/dir) into a labelled
    peptide database that synteny search can run on."""
    from pynteny.api import Build

    data_path = Path(data)
    if not data_path.exists():
        return {"error": f"Input data not found: {data}"}

    build = Build(
        data=data_path,
        prepend=prepend,
        outfile=outfile,
        logfile=logfile,
        processes=processes,
        tempdir=tempdir,
    )
    with _pynteny_run():
        build.run()

    out = Path(build._args.outfile) if build._args.outfile else None
    result: dict[str, Any] = {"input": str(data_path)}
    if out is not None and out.exists():
        result["output_file"] = str(out)
        result["output_size_bytes"] = out.stat().st_size
        result["n_peptides"] = _count_fasta_records(out)
    else:
        result["output_file"] = str(out) if out else None
        result["warning"] = (
            "Build completed but the output file could not be located; "
            "check the logfile."
        )
    return result


# --------------------------------------------------------------------------- #
# Search                                                                      #
# --------------------------------------------------------------------------- #
def run_search(
    data: str,
    synteny_structure: str,
    *,
    gene_ids: bool,
    unordered: bool,
    best_hmm_wins: bool,
    reuse: bool,
    hmm_dir: Optional[str],
    hmm_meta: Optional[str],
    outdir: Optional[str],
    prefix: str,
    hmmsearch_args: Optional[str],
    processes: Optional[int],
    logfile: Optional[str],
    max_hits: int = 200,
) -> dict[str, Any]:
    """Run a synteny-aware HMM search over a labelled peptide database and return
    the matched hits plus the paths of the files Pynteny wrote."""
    from pynteny.api import Search

    data_path = Path(data)
    if not data_path.exists():
        return {"error": f"Sequence database not found: {data}"}

    resolved_dir = _resolve_hmm_dir(hmm_dir)
    if resolved_dir and not Path(resolved_dir).exists():
        return {"error": f"HMM directory not found: {resolved_dir}"}
    if resolved_dir is None:
        return {
            "error": (
                "No HMM directory available. Provide `hmm_dir`, or download the "
                "PGAP database first with `download_hmm_databases` (the search "
                "would otherwise trigger a large download)."
            )
        }
    resolved_meta = _resolve_hmm_meta(hmm_meta)

    kwargs = dict(
        data=data_path,
        synteny_struc=synteny_structure,
        gene_ids=gene_ids,
        unordered=unordered,
        reuse=reuse,
        hmm_dir=resolved_dir,
        hmm_meta=resolved_meta,
        outdir=outdir,
        prefix=prefix,
        hmmsearch_args=hmmsearch_args,
        logfile=logfile,
        processes=processes,
    )
    # `best_hmm_wins` was added after the published 1.2.0 release; pass it only
    # when the installed Pynteny supports it so this works on both.
    import inspect

    best_hmm_wins_supported = "best_hmm_wins" in inspect.signature(Search.__init__).parameters
    if best_hmm_wins_supported:
        kwargs["best_hmm_wins"] = best_hmm_wins

    search = Search(**kwargs)

    with _pynteny_run():
        synteny_hits = search.run()
        if resolved_meta and Path(resolved_meta).exists():
            annotated = synteny_hits.add_HMM_meta_info_to_hits(resolved_meta)
            # Depending on the Pynteny version this returns either a SyntenyHits
            # or the underlying DataFrame directly.
            df = annotated.hits if hasattr(annotated, "hits") else annotated
        else:
            df = synteny_hits.hits

    out_dir = Path(search._args.outdir)
    records = df.where(df.notna(), None).to_dict(orient="records")
    fasta_files = sorted(str(p) for p in out_dir.glob(f"{prefix}*_hits.fasta"))
    synteny_table = out_dir / f"{prefix}synteny_matched.tsv"

    result: dict[str, Any] = {
        "synteny_structure": search._args.synteny_struc,
        "data": str(data_path),
        "n_hits": int(len(df)),
        "columns": list(df.columns),
        "hits": records[:max_hits],
        "hits_truncated": len(records) > max_hits,
        "output_dir": str(out_dir),
        "synteny_table": str(synteny_table) if synteny_table.exists() else None,
        "fasta_files": fasta_files,
    }
    if best_hmm_wins and not best_hmm_wins_supported:
        result["warning"] = (
            "best_hmm_wins was requested but the installed Pynteny does not "
            "support it; the option was ignored."
        )
    return result


# --------------------------------------------------------------------------- #
# Download                                                                    #
# --------------------------------------------------------------------------- #
def download_databases(
    outdir: str,
    *,
    pgap: bool,
    pfam: bool,
    unpack: bool,
    force: bool,
    logfile: Optional[str],
) -> dict[str, Any]:
    """Download PGAP and/or PFAM profile-HMM databases from NCBI / InterPro.

    The :class:`pynteny.api.Download` wrapper hard-wires PGAP only, so we drive
    the underlying ``download_hmms`` subcommand directly to expose the PFAM
    option too.
    """
    from pynteny.subcommands import download_hmms
    from pynteny.utils import CommandArgs, ConfigParser

    if not (pgap or pfam):
        return {"error": "Select at least one of pgap / pfam to download."}

    args = CommandArgs(
        outdir=Path(outdir),
        logfile=Path(logfile) if logfile else None,
        force=force,
        unpack=unpack,
        pgap=pgap,
        pfam=pfam,
    )
    try:
        with _pynteny_run():
            download_hmms(args)
    except SystemExit:
        # download_hmms exits(1) when the requested databases are already present.
        pass

    config = ConfigParser.get_default_config()
    return {
        "outdir": str(Path(outdir).absolute()),
        "requested": {"pgap": pgap, "pfam": pfam},
        "hmm_databases": {
            "PGAP_data_downloaded": config.get_field("PGAP_data_downloaded"),
            "PFAM_data_downloaded": config.get_field("PFAM_data_downloaded"),
            "PGAP_database": config.get_field("PGAP_database"),
            "PGAP_meta_file": config.get_field("PGAP_meta_file"),
            "PFAM_database": config.get_field("PFAM_database"),
            "PFAM_meta_file": config.get_field("PFAM_meta_file"),
        },
    }
