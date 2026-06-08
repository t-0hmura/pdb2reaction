"""FastMCP server entry point for pdb2reaction.

Run as console_script:

    pdb2reaction-mcp        # equivalent to: python -m pdb2reaction.mcp.server
    p2r-mcp                 # short alias

The server speaks JSON-RPC over stdio (the default MCP transport). Register
in an MCP client (Claude Desktop, Cursor, etc.) via:

    {
      "mcpServers": {
        "pdb2reaction": {
          "command": "pdb2reaction-mcp",
          "args": []
        }
      }
    }

See `docs/mcp_server.md` for the full tool list, parameter schemas, and
return-value structure.
"""
from __future__ import annotations

from mcp.server.fastmcp import FastMCP

from pdb2reaction.mcp._tools import register_all

# One server instance, registered tools attached at import time.
mcp = FastMCP("pdb2reaction")
register_all(mcp)


def serve() -> None:
    """Run the MCP server on stdio (blocking call)."""
    mcp.run()


if __name__ == "__main__":
    serve()
