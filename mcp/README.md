# Pynteny MCP server

A [Model Context Protocol](https://modelcontextprotocol.io) server that exposes
the [Pynteny](../) API as tools, so any MCP-compatible agent — Claude Desktop,
Claude Code, or a custom client driving Claude / DeepSeek / etc. — can run
synteny-aware HMM searches over genomic sequence data in natural language.

It wraps `pynteny.api` (`Search` / `Build` / `Download`) and the synteny-structure
parsers, returning **compact JSON summaries** (the matched hits table and the
paths of the files Pynteny writes) instead of dumping large objects.

```
mcp/
├── pyproject.toml             # installable package: `pynteny-mcp`
├── requirements.txt           # server + example deps
├── .env.example               # copy to .env and fill in keys (gitignored)
├── src/pynteny_mcp/
│   ├── server.py              # FastMCP server + tool definitions
│   └── service.py             # API wrappers, stdout-safe logging, JSON summaries
├── examples/
│   └── synteny_search_agent.py    # LLM agent (Claude + DeepSeek) demo
└── tests/
    └── smoke_test.py          # no-API-key end-to-end check against the test data
```

## What is a synteny structure?

A synteny structure describes a target gene arrangement, e.g.

```
>leuD 0 >leuC 1 <leuA
```

where each token is an **HMM name** (or a **gene symbol**, with `gene_ids=true`),
the leading `>`/`<` gives the **strand** (sense / antisense), and the integers are
the **maximum number of (untargeted) genes** allowed between neighbours. Groups of
interchangeable HMMs for one gene are written `(HMM_A|HMM_B)`.

## Tools

| Tool | What it does |
|------|--------------|
| `get_pynteny_info` | Pynteny version, citation, and which HMM databases (PGAP/PFAM) are downloaded |
| `validate_synteny_structure` | Validate & decompose a structure (HMMs, strands, distances) — no I/O |
| `parse_gene_symbols` | Translate a gene-symbol structure into HMM names via a metadata table |
| `build_peptide_database` | Predict ORFs and label them from a nucleotide assembly / GenBank |
| `run_synteny_search` | **Core:** run the synteny-aware HMM search; return matched hits + output paths (incl. `best_hmm_wins` for paralog cross-hits) |
| `download_hmm_databases` | Download the PGAP and/or PFAM profile-HMM databases |

## Prerequisites

A Python ≥ 3.8 environment with **Pynteny** installed. Pynteny is now a
pure-Python pip package — HMMER and Prodigal come bundled via
[pyhmmer](https://github.com/althonos/pyhmmer) and
[pyrodigal](https://github.com/althonos/pyrodigal), so **no conda environment or
external binaries are required**:

```bash
pip install pynteny            # or: pip install git+https://github.com/Robaina/Pynteny.git
```

## Install

Install the server (and example) dependencies into that same environment:

```bash
cd mcp
pip install -r requirements.txt          # mcp, pynteny, anthropic, openai, python-dotenv
# optional — register the `pynteny-mcp` console script:
pip install -e .
```

## Configure

Only the example agent needs API keys; the server itself does not.

```bash
cp .env.example .env        # .env is gitignored
```

Edit `.env`:

```ini
ANTHROPIC_API_KEY=sk-ant-...
DEEPSEEK_API_KEY=sk-...
# optional defaults so the agent need not repeat paths:
# PYNTENY_DATA=/abs/path/to/labelled_peptides.faa
# PYNTENY_HMM_DIR=/abs/path/to/data/hmms
# PYNTENY_HMM_META=/abs/path/to/data/hmms/hmm_PGAP.tsv
```

## Verify (no API keys needed)

Runs the server against Pynteny's committed test genome and asserts the known
*leu*-operon synteny hits:

```bash
python tests/smoke_test.py
# → "All smoke-test checks passed."
```

## Run the example agent

An LLM decides how to validate the structure and run the search to answer a
question. By default it searches Pynteny's test genome
(`../tests/test_data/MG1655.fasta`), so it works with no extra data; set
`PYNTENY_DATA` / `PYNTENY_HMM_DIR` / `PYNTENY_HMM_META` to search your own.

```bash
python examples/synteny_search_agent.py --provider claude
python examples/synteny_search_agent.py --provider deepseek
python examples/synteny_search_agent.py --provider claude \
    --question "Is leuD-leuC-leuA syntenic in this genome, and on which strand?"
```

The example launches the MCP server itself (as a stdio subprocess using the same
Python interpreter), connects as an MCP client, converts the MCP tool schemas to
each provider's tool format, and runs a manual tool-use loop.

> **Models.** Claude defaults to `claude-opus-4-8` with adaptive thinking;
> DeepSeek defaults to `deepseek-v4-pro`. Override via `ANTHROPIC_MODEL` /
> `DEEPSEEK_MODEL` in `.env`.

## A typical end-to-end workflow

For real data (not the bundled test genome) the agent (or you) would:

1. `download_hmm_databases(outdir="data/hmms")` — once, to fetch PGAP HMMs + metadata.
2. `build_peptide_database(data="assembly.fa", outfile="labelled_peptides.faa")` — predict & label ORFs.
3. `validate_synteny_structure(">leuD 0 >leuC 1 <leuA")` — sanity-check the query.
4. `run_synteny_search(data="labelled_peptides.faa", synteny_structure=">leuD 0 >leuC 1 <leuA", gene_ids=true)` — search.

## Register with an MCP client

### Claude Desktop / Claude Code

Add to your MCP config (`claude_desktop_config.json`, or Claude Code's
`.mcp.json`). Use the absolute interpreter path of the environment where
`pynteny` + `mcp` are installed, and put `src/` on `PYTHONPATH`:

```json
{
  "mcpServers": {
    "pynteny": {
      "command": "/path/to/python",
      "args": ["-m", "pynteny_mcp.server"],
      "env": {
        "PYTHONPATH": "/home/robaina/Documents/Hapdera/Hapdera-Projects/Pynteny/mcp/src"
      }
    }
  }
}
```

If you ran `pip install -e .`, you can use the console script instead:
`"command": ".../bin/pynteny-mcp"` with no `args`/`PYTHONPATH`.

## Notes

- **stdout is sacred.** The stdio transport uses stdout for protocol traffic, but
  Pynteny logs and prints to stdout. The service layer wraps every Pynteny call
  so that its stdout (and logging) is redirected to stderr, keeping the protocol
  stream clean. Set `PYNTENY_MCP_DEBUG=1` to restore the MCP runtime's verbose
  per-request logging.
- **No big in-memory database.** Unlike a database server, each tool call is a
  Pynteny operation; `run_synteny_search` runs HMMER and writes per-gene FASTA
  files plus a `synteny_matched.tsv` table to `outdir`.
- **`best_hmm_wins`.** When one peptide is hit by several HMMs (paralog
  cross-hits), this keeps only the highest-scoring HMM for that peptide before
  matching the structure. It requires Pynteny ≥ 1.3.0; on an older build the
  service ignores it and returns a `warning` field rather than failing.
- **DeepSeek** is reached through its OpenAI-compatible endpoint, so the example
  uses the `openai` SDK pointed at `DEEPSEEK_BASE_URL`. The same MCP server and
  tool set serve both providers unchanged.
