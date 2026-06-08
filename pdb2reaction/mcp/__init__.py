"""MCP (Model Context Protocol) server for pdb2reaction.

Exposes every `pdb2reaction <subcmd>` CLI as an MCP tool, callable by LLM
agents (Claude Desktop, Cursor, Codeium, custom MCP clients) via JSON-RPC
over stdio.

Each tool wraps the corresponding CLI invocation via `subprocess.run`,
parses the produced `summary.json`, and returns a structured dict. No
Python library API is required — the MCP layer is the sole programmatic
entry point besides the CLI itself.

Entry point: console_script `pdb2reaction-mcp` (alias `p2r-mcp`).
"""
from __future__ import annotations

__all__ = ["serve"]


def serve() -> None:
    """Console-script entry point.

    Lazily imports `pdb2reaction.mcp.server` so that `pdb2reaction --help`
    does not pull `mcp` into the import graph.
    """
    from pdb2reaction.mcp.server import serve as _serve

    _serve()
