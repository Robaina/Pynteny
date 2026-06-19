#!/usr/bin/env python
"""
No-API-key smoke test for the Pynteny MCP server.

Launches the server over stdio, connects as an MCP client, and exercises the
tools end to end against Pynteny's own committed test data
(``tests/test_data/MG1655.fasta`` + the ``hmms/`` directory), asserting on the
known synteny hits for the *leu* operon. No LLM and no API keys are involved —
this validates the MCP plumbing and the tool logic, including that Pynteny's
stdout logging does not corrupt the stdio protocol.

    python tests/smoke_test.py
"""

from __future__ import annotations

import asyncio
import json
import os
import sys
import tempfile
from datetime import timedelta
from pathlib import Path

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

HERE = Path(__file__).resolve().parent
MCP_DIR = HERE.parent
SRC_DIR = MCP_DIR / "src"
REPO_ROOT = MCP_DIR.parent
TEST_DATA = REPO_ROOT / "tests" / "test_data"

# The leu-operon synteny structure used by Pynteny's own integration test.
SYNTENY_STRUC = (
    "<(TIGR00171.1|TIGR02084.1) 0 "
    "<(TIGR00170.1|TIGR02083.1) 1 "
    "<(TIGR00973.1|NF002084.0|TIGR00970.1)"
)
EXPECTED_LABELS = {
    "b0071__U00096_71_78847_79453_neg",
    "b0072__U00096_72_79463_80864_neg",
    "b0074__U00096_74_81957_83529_neg",
}

TIMEOUT = timedelta(seconds=300)


async def call(session: ClientSession, name: str, **arguments):
    result = await session.call_tool(name, arguments, read_timeout_seconds=TIMEOUT)
    assert not getattr(result, "isError", False), f"{name} returned an error: {result.content}"
    text = "\n".join(getattr(b, "text", "") for b in result.content if getattr(b, "text", ""))
    return json.loads(text)


def check(label: str, condition: bool, detail: str = "") -> None:
    status = "PASS" if condition else "FAIL"
    print(f"  [{status}] {label}" + (f" — {detail}" if detail else ""))
    if not condition:
        raise AssertionError(f"{label}: {detail}")


async def main() -> None:
    data = TEST_DATA / "MG1655.fasta"
    hmm_dir = TEST_DATA / "hmms"
    hmm_meta = TEST_DATA / "hmm_meta.tsv"
    for p in (data, hmm_dir, hmm_meta):
        if not p.exists():
            sys.exit(f"Test data not found: {p}")

    server_params = StdioServerParameters(
        command=sys.executable,
        args=["-m", "pynteny_mcp.server"],
        env={
            "PYTHONPATH": str(SRC_DIR),
            "PATH": os.environ.get("PATH", ""),
        },
    )

    with tempfile.TemporaryDirectory() as outdir:
        async with stdio_client(server_params) as (read, write):
            async with ClientSession(read, write) as session:
                await session.initialize()
                tools = (await session.list_tools()).tools
                print(f"Server exposes {len(tools)} tools: {', '.join(t.name for t in tools)}\n")
                check("expected tool count", len(tools) == 6, f"got {len(tools)}")

                info = await call(session, "get_pynteny_info")
                check("pynteny version reported", bool(info.get("version")), str(info.get("version")))
                check("citation present", "Pynteny" in info.get("citation", ""), "")

                good = await call(session, "validate_synteny_structure",
                                  synteny_structure=SYNTENY_STRUC)
                check("valid structure recognised", good["valid"] is True, str(good.get("valid")))
                check("structure has 3 genes", good["n_genes"] == 3, str(good.get("n_genes")))
                check("max distances parsed", good["max_distances"] == [0, 1],
                      str(good.get("max_distances")))
                check("strands parsed", good["strands"] == ["neg", "neg", "neg"],
                      str(good.get("strands")))

                # Two adjacent genes with no distance token between them: the
                # gene/distance counts don't line up, so this is malformed.
                bad = await call(session, "validate_synteny_structure",
                                 synteny_structure=">TIGR00171.1 >TIGR00170.1")
                check("invalid structure rejected", bad["valid"] is False, str(bad.get("valid")))

                hits = await call(
                    session, "run_synteny_search",
                    data=str(data),
                    synteny_structure=SYNTENY_STRUC,
                    hmm_dir=str(hmm_dir),
                    hmm_meta=str(hmm_meta),
                    reuse=True,
                    outdir=outdir,
                )
                check("search returned 3 hits", hits["n_hits"] == 3, str(hits.get("n_hits")))
                labels = {h.get("full_label") for h in hits["hits"]}
                check("hit labels match expected leu operon", labels == EXPECTED_LABELS,
                      str(labels))
                check("hits annotated with gene_symbol", "gene_symbol" in hits["columns"],
                      str(hits["columns"]))
                check("wrote synteny_matched.tsv", bool(hits.get("synteny_table")),
                      str(hits.get("synteny_table")))
                check("wrote per-gene FASTA files", len(hits.get("fasta_files", [])) >= 1,
                      str(hits.get("fasta_files")))

                # best_hmm_wins: must still return the leu hits. On a Pynteny old
                # enough to lack the option the service reports a `warning` and
                # carries on; on a new-enough Pynteny it is honored (no warning).
                bhw = await call(
                    session, "run_synteny_search",
                    data=str(data),
                    synteny_structure=SYNTENY_STRUC,
                    hmm_dir=str(hmm_dir),
                    hmm_meta=str(hmm_meta),
                    reuse=True,
                    best_hmm_wins=True,
                    outdir=outdir,
                )
                check("best_hmm_wins search still returns 3 hits", bhw["n_hits"] == 3,
                      str(bhw.get("n_hits")))
                if bhw.get("warning"):
                    print(f"  [NOTE] best_hmm_wins not honored by installed Pynteny: "
                          f"{bhw['warning']}")
                else:
                    print("  [NOTE] best_hmm_wins honored by installed Pynteny.")

    print("\nAll smoke-test checks passed.")


if __name__ == "__main__":
    asyncio.run(main())
