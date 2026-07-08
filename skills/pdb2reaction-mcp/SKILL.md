---
name: pdb2reaction-mcp
description: How to drive `pdb2reaction` from an MCP-speaking agent (Claude Desktop / Claude Code / Cursor / Codeium / custom Python or TypeScript MCP SDK clients) via the bundled `pdb2reaction-mcp` (alias `p2r-mcp`) server. Lists the 18 MCP tools (one per CLI subcommand) and the shared `SubcmdResult` schema. TRIGGER on questions about MCP setup, agent integration, tool names, return-value shape, or "how does an agent invoke pdb2reaction". SKIP for direct CLI usage — that's `pdb2reaction-cli`.
---

# pdb2reaction MCP server

`pdb2reaction-mcp` exposes every CLI subcommand as an MCP tool over stdio JSON-RPC. Any MCP client can drive the full reaction-path pipeline without spawning shell processes manually.

## Install + register

```bash
pip install "pdb2reaction[mcp]"   # adds the `mcp[cli]` extra
```

This registers two console scripts: `pdb2reaction-mcp` and `p2r-mcp` (alias).

Drop the snippet below into your client's MCP config file:

- Claude Desktop — `~/Library/Application Support/Claude/claude_desktop_config.json` (macOS) / `%APPDATA%\Claude\claude_desktop_config.json` (Windows)
- Cursor — `~/.cursor/mcp.json`
- Other MCP clients — consult the client's own MCP-server docs

```json
{
  "mcpServers": {
    "pdb2reaction": {
      "command": "pdb2reaction-mcp",
      "args": []
    }
  }
}
```

See [`examples/mcp_client_config.json`](../../examples/mcp_client_config.json) for a full sample with explicit `PATH` / `CUDA_VISIBLE_DEVICES` env-var overrides.

## 18 tools (1 per CLI subcommand)

### Stage runners

| MCP tool | wraps | purpose |
|---|---|---|
| `optimize_geometry` | `pdb2reaction opt` | Optimise one structure |
| `find_transition_state` | `pdb2reaction tsopt` | TS search (RS-P-RFO / Dimer / TRIM / RS-I-RFO) |
| `run_irc` | `pdb2reaction irc` | IRC integration from a TS geometry |
| `compute_frequencies` | `pdb2reaction freq` | Vibrational analysis + thermochemistry |
| `run_single_point` | `pdb2reaction sp` | Single-point MLIP energy + forces (+optional Hessian) |

### Scans / paths / pipeline

| MCP tool | wraps | purpose |
|---|---|---|
| `scan_1d` / `scan_2d` / `scan_3d` | `pdb2reaction scan{,2d,3d}` | Restraint-driven distance scans |
| `optimize_path` | `pdb2reaction path-opt` | Two-endpoint MEP optimisation |
| `search_paths` | `pdb2reaction path-search` | Recursive reaction-pathway search |
| `run_full_pipeline` | `pdb2reaction all` | End-to-end (extract → MEP → TS → IRC → freq → DFT) |
| `run_single_point_dft` | `pdb2reaction dft` | Single-point DFT via gpu4pyscf |

### Structure / I/O helpers

| MCP tool | wraps | purpose |
|---|---|---|
| `extract_active_site` | `pdb2reaction extract` | Cut a sphere around a ligand |
| `add_element_info` | `pdb2reaction add-elem-info` | Repair PDB element columns |
| `fix_altloc` | `pdb2reaction fix-altloc` | Resolve PDB alternate locations |
| `plot_trajectory` | `pdb2reaction trj2fig` | Energy profile PNG / SVG / PDF / HTML / CSV |
| `plot_energy_diagram` | `pdb2reaction energy-diagram` | Categorical energy diagram |
| `detect_bond_changes` | `pdb2reaction bond-summary` | Bond-change diff between two PDBs |

## `SubcmdResult` return schema

Every tool returns the same structured dict so the calling agent can dispatch on `status` without parsing stderr:

```python
{
    "schema_version": "1.0",         # pin to MCP_SUBCMD_RESULT_SCHEMA_VERSION
    "status": "ok" | "failed" | "summary_missing" | "summary_parse_error",
    "exit_code": int,                # subprocess exit code
    "out_dir": str | None,           # working directory the CLI wrote to
    "summary": dict,                 # parsed summary.json (CLI output schema)
    "stderr_tail": str,              # last ~60 lines of stderr
    "stdout_tail": str,              # last ~60 lines of stdout
    "hint": str | None,              # parsed `; recover: <hint>` suffix
    "argv": list[str],               # full argv that was executed (reproducible)
}
```

A failed subcommand additionally surfaces a structured exception envelope inside `summary`:

- `error_class_chain` — list of class names walking the MRO (e.g. `["CudaOutOfMemoryError", "RuntimeError", "Exception"]`)
- `error_module` — module path of the originating exception class
- `error_label` — high-level CLI stage label (e.g. `"IRC"`, `"optimization"`, `"frequency analysis"`, `"path optimization"`, fallback `"UnhandledError"`)

so an MCP client can pattern-match on the hierarchy instead of substring-matching `stderr_tail`.

## `find_transition_state` opt-mode kwarg

Beyond the default RS-P-RFO, `find_transition_state` accepts:

- `opt_mode="trim"` — Helgaker trust-region image-minimisation TS opt
- `opt_mode="rsprfo"` — Banerjee restricted-step P-RFO TS opt
- `opt_mode="dimer"` (alias `grad`) — Hessian-guided Dimer TS opt

## Sandbox / safety notes

- The MCP server inherits the calling environment's PATH, conda env, and CUDA setup. Long-running tools (opt / tsopt / irc) launch the `pdb2reaction` CLI in a subprocess — set `timeout_seconds` on each call to bound runaway computations.
- Output files land under the `out_dir` kwarg (defaults to a unique `tempfile.mkdtemp("p2r_mcp_<subcmd>_…")` so concurrent agent calls don't collide).
- The server does not modify `~/.bashrc` / login env, install software, or write outside `out_dir`. All MLIP weights / PDB inputs must already exist on disk.

## See also
- Full MCP server doc: [`docs/mcp_server.md`](../../docs/mcp_server.md)
- Sample MCP-client config: [`examples/mcp_client_config.json`](../../examples/mcp_client_config.json)
- Server / tool source: `pdb2reaction/mcp/{server,_tools,_runner}.py`
