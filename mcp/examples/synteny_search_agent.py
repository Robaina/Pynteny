#!/usr/bin/env python
"""
Example: an LLM agent that runs Pynteny synteny searches through the MCP server.

Instead of hand-writing the Pynteny calls, an LLM is given the MCP toolbox and
decides how to validate a synteny structure and run the search to answer a
question — e.g. "find the leucine-biosynthesis (leuD-leuC-leuA) synteny block in
this genome".

The same MCP toolbox is driven by two providers so you can compare them:

* Claude   — via the official `anthropic` SDK (`claude-opus-4-8` by default)
* DeepSeek — via the OpenAI-compatible `openai` SDK

Both connect to the *same* MCP server (``pynteny_mcp.server``) launched as a
subprocess over stdio; only the model and the tool-schema adapter differ.

Usage:
    cp .env.example .env      # then fill in your API keys (+ optional PYNTENY_* paths)
    python examples/synteny_search_agent.py --provider claude
    python examples/synteny_search_agent.py --provider deepseek
    python examples/synteny_search_agent.py --provider claude --question "..."

By default it runs against Pynteny's committed test genome
(../tests/test_data/MG1655.fasta) so it works with no extra data; point
PYNTENY_DATA / PYNTENY_HMM_DIR / PYNTENY_HMM_META at your own files to search
real data.
"""

from __future__ import annotations

import argparse
import asyncio
import json
import os
import sys
from datetime import timedelta
from pathlib import Path

from dotenv import load_dotenv
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

HERE = Path(__file__).resolve().parent
MCP_DIR = HERE.parent
SRC_DIR = MCP_DIR / "src"
REPO_ROOT = MCP_DIR.parent
ENV_PATH = MCP_DIR / ".env"

load_dotenv(ENV_PATH)

# Fall back to Pynteny's committed test data so the demo runs out of the box.
TEST_DATA = REPO_ROOT / "tests" / "test_data"
DEFAULT_DATA = os.environ.get("PYNTENY_DATA", str(TEST_DATA / "MG1655.fasta"))
DEFAULT_HMM_DIR = os.environ.get("PYNTENY_HMM_DIR", str(TEST_DATA / "hmms"))
DEFAULT_HMM_META = os.environ.get("PYNTENY_HMM_META", str(TEST_DATA / "hmm_meta.tsv"))

# Per-tool-call timeout: a search runs HMMER over the whole database.
TOOL_TIMEOUT = timedelta(seconds=600)
MAX_TURNS = 14  # guard against runaway tool loops

SYSTEM_PROMPT = (
    "You are a comparative-genomics assistant with access to Pynteny through "
    "tools. Pynteny finds synteny blocks: genes hit by profile HMMs whose "
    "genomic arrangement (order, strand, spacing) matches a synteny structure "
    "like '>leuD 0 >leuC 1 <leuA'. Use the tools rather than prior knowledge: "
    "validate a structure before searching, then run the search and report the "
    "matched genes (contig, gene id, strand, HMM, and gene symbol/product when "
    "available) and where the result files were written. Be concise.\n\n"
    f"Unless told otherwise, search this database: {DEFAULT_DATA}\n"
    f"Using this HMM directory: {DEFAULT_HMM_DIR}\n"
    f"And this HMM metadata file: {DEFAULT_HMM_META}\n"
    "Pass reuse=true so repeated runs don't recompute HMMER outputs."
)

DEFAULT_QUESTION = (
    "Find the leucine-biosynthesis synteny block in the genome. The genes "
    "leuD, leuC and leuA are expected adjacent on the antisense strand. Their "
    "HMMs are leuD=(TIGR00171.1|TIGR02084.1), leuC=(TIGR00170.1|TIGR02083.1), "
    "leuA=(TIGR00973.1|NF002084.0|TIGR00970.1). Allow 0 genes between leuD and "
    "leuC and up to 1 between leuC and leuA. Which genes match, and on which "
    "contig?"
)


# --------------------------------------------------------------------------- #
# MCP toolbox: connect to the server and adapt its tools to each provider     #
# --------------------------------------------------------------------------- #
class MCPToolbox:
    def __init__(self, session: ClientSession):
        self.session = session
        self.tools = []

    async def load(self) -> None:
        self.tools = (await self.session.list_tools()).tools

    async def call(self, name: str, arguments: dict) -> str:
        result = await self.session.call_tool(
            name, arguments, read_timeout_seconds=TOOL_TIMEOUT
        )
        parts = [
            getattr(b, "text", "") for b in result.content if getattr(b, "text", "")
        ]
        text = "\n".join(parts).strip()
        if not text and getattr(result, "structuredContent", None):
            text = json.dumps(result.structuredContent)
        if getattr(result, "isError", False):
            return f"ERROR: {text or '(unknown tool error)'}"
        return text or "(no content)"

    def anthropic_tools(self) -> list[dict]:
        return [
            {
                "name": t.name,
                "description": t.description or "",
                "input_schema": t.inputSchema,
            }
            for t in self.tools
        ]

    def openai_tools(self) -> list[dict]:
        return [
            {
                "type": "function",
                "function": {
                    "name": t.name,
                    "description": t.description or "",
                    "parameters": t.inputSchema,
                },
            }
            for t in self.tools
        ]


def _log_tool_call(name: str, arguments: dict) -> None:
    print(f"  → tool: {name}({json.dumps(arguments, ensure_ascii=False)})")


# --------------------------------------------------------------------------- #
# Claude backend (official anthropic SDK, manual agentic loop)                #
# --------------------------------------------------------------------------- #
async def run_claude(toolbox: MCPToolbox, question: str) -> str:
    from anthropic import AsyncAnthropic

    model = os.environ.get("ANTHROPIC_MODEL", "claude-opus-4-8")
    effort = os.environ.get("ANTHROPIC_EFFORT", "medium")
    client = AsyncAnthropic()  # reads ANTHROPIC_API_KEY from the environment
    tools = toolbox.anthropic_tools()
    messages: list[dict] = [{"role": "user", "content": question}]

    for _ in range(MAX_TURNS):
        resp = await client.messages.create(
            model=model,
            max_tokens=16000,
            system=SYSTEM_PROMPT,
            tools=tools,
            # Adaptive thinking lets Claude decide how much to reason between
            # tool calls; effort trades thoroughness against token cost.
            thinking={"type": "adaptive", "display": "summarized"},
            output_config={"effort": effort},
            messages=messages,
        )

        for block in resp.content:
            if block.type == "thinking" and getattr(block, "thinking", ""):
                print(f"  [thinking] {block.thinking.strip()[:300]}")
            elif block.type == "text" and block.text.strip():
                print(f"  [claude] {block.text.strip()}")

        if resp.stop_reason != "tool_use":
            return "".join(b.text for b in resp.content if b.type == "text").strip()

        messages.append({"role": "assistant", "content": resp.content})
        tool_results = []
        for block in resp.content:
            if block.type == "tool_use":
                _log_tool_call(block.name, block.input)
                output = await toolbox.call(block.name, block.input)
                tool_results.append(
                    {"type": "tool_result", "tool_use_id": block.id, "content": output}
                )
        messages.append({"role": "user", "content": tool_results})

    return "(stopped: reached the maximum number of tool-use turns)"


# --------------------------------------------------------------------------- #
# DeepSeek backend (OpenAI-compatible API, manual function-calling loop)       #
# --------------------------------------------------------------------------- #
async def run_deepseek(toolbox: MCPToolbox, question: str) -> str:
    from openai import AsyncOpenAI

    model = os.environ.get("DEEPSEEK_MODEL", "deepseek-v4-pro")
    client = AsyncOpenAI(
        api_key=os.environ["DEEPSEEK_API_KEY"],
        base_url=os.environ.get("DEEPSEEK_BASE_URL", "https://api.deepseek.com"),
    )
    tools = toolbox.openai_tools()
    messages: list[dict] = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {"role": "user", "content": question},
    ]

    for _ in range(MAX_TURNS):
        resp = await client.chat.completions.create(
            model=model, messages=messages, tools=tools
        )
        msg = resp.choices[0].message

        if msg.content and msg.content.strip():
            print(f"  [deepseek] {msg.content.strip()}")

        assistant: dict = {"role": "assistant", "content": msg.content or ""}
        if msg.tool_calls:
            assistant["tool_calls"] = [
                {
                    "id": tc.id,
                    "type": "function",
                    "function": {
                        "name": tc.function.name,
                        "arguments": tc.function.arguments,
                    },
                }
                for tc in msg.tool_calls
            ]
        messages.append(assistant)

        if not msg.tool_calls:
            return (msg.content or "").strip()

        for tc in msg.tool_calls:
            try:
                arguments = json.loads(tc.function.arguments or "{}")
            except json.JSONDecodeError:
                arguments = {}
            _log_tool_call(tc.function.name, arguments)
            output = await toolbox.call(tc.function.name, arguments)
            messages.append({"role": "tool", "tool_call_id": tc.id, "content": output})

    return "(stopped: reached the maximum number of tool-use turns)"


# --------------------------------------------------------------------------- #
# Driver                                                                       #
# --------------------------------------------------------------------------- #
async def main_async(provider: str, question: str) -> None:
    # Launch the MCP server as a subprocess, using the *same* interpreter (so it
    # shares this environment's `pynteny` + `mcp`), with src/ on PYTHONPATH.
    server_params = StdioServerParameters(
        command=sys.executable,
        args=["-m", "pynteny_mcp.server"],
        env={**os.environ, "PYTHONPATH": str(SRC_DIR)},
    )

    async with stdio_client(server_params) as (read, write):
        async with ClientSession(read, write) as session:
            await session.initialize()
            toolbox = MCPToolbox(session)
            await toolbox.load()

            print(
                f"Connected to Pynteny MCP server: {len(toolbox.tools)} tools available"
            )
            print(f"Provider: {provider}\n")
            print(f"Question:\n  {question}\n")
            print("--- agent trace ---")

            try:
                if provider == "claude":
                    answer = await run_claude(toolbox, question)
                else:
                    answer = await run_deepseek(toolbox, question)
            except (
                Exception
            ) as exc:  # noqa: BLE001 — surface a clean message, not a stack trace
                print(f"\nLLM call failed ({type(exc).__name__}): {exc}")
                return

            print("\n--- final answer ---")
            print(answer)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--provider",
        choices=["claude", "deepseek"],
        default="claude",
        help="Which LLM backend to use (default: claude).",
    )
    parser.add_argument(
        "--question",
        default=DEFAULT_QUESTION,
        help="Question to ask the agent (defaults to the leu-operon example).",
    )
    args = parser.parse_args()

    key_var = "ANTHROPIC_API_KEY" if args.provider == "claude" else "DEEPSEEK_API_KEY"
    if not os.environ.get(key_var):
        sys.exit(f"{key_var} is not set. Add it to {ENV_PATH} (see .env.example).")

    asyncio.run(main_async(args.provider, args.question))


if __name__ == "__main__":
    main()
