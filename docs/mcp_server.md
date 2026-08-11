# pdb2reaction MCP server

`pdb2reaction-mcp` (alias `p2r-mcp`) is an [MCP](https://modelcontextprotocol.io/)
server that lets any MCP-speaking agent drive every `pdb2reaction` CLI
subcommand via JSON-RPC over stdio — including Claude Desktop / Claude Code /
Cursor / Codeium, and any custom agent built on the official Python or
TypeScript MCP SDKs.

## Install

```bash
pip install "pdb2reaction[mcp]"
```

This adds the `mcp[cli]` dependency and registers two console scripts:
`pdb2reaction-mcp` and `p2r-mcp` (alias).

## Tools

18 tools, one per CLI subcommand. Each tool returns a structured dict with:

- `schema_version`: envelope version. Live value: `pdb2reaction.mcp._runner.MCP_SUBCMD_RESULT_SCHEMA_VERSION`. A version bump signals a field-set / value-type change; pin against the constant rather than the literal in this doc.
- `status`: `ok` | `failed` | `summary_missing` | `summary_parse_error` | `summary_run_mismatch`
- `exit_code`: subprocess exit code
- `out_dir`: managed stage directory, or null for successful helper tools
- `summary`: parsed `summary.json`; successful helpers without a managed summary return an empty object
- `stderr_tail` / `stdout_tail`: last ~60 lines of process output
- `hint`: parsed `; recover: <hint>` suffix from CLI error messages, if any
- `argv`: the full argv that was executed (for reproducibility)
- `run_id`: UUID generated for this subprocess invocation. A parsed summary is
  returned only when it carries this exact ID; stale or mixed-generation
  output is reported as `summary_run_mismatch` with an empty `summary`.

Product CLI aliases are bound to the MCP server's current Python interpreter
as `python -m pdb2reaction`, and the imported package source root is first in
the child `PYTHONPATH`. The child keeps the caller's current working directory,
so relative scientific input paths retain their original meaning.

### Structured error envelope

When a subcommand fails, the parsed `summary` (or sibling `result.json`) carries an extended error envelope so agents can pattern-match the exception class hierarchy without parsing text:

- `error`: `str(exc)` of the original exception
- `error_type`: exception class name (e.g. `"OptimizationError"`)
- `error_class_chain`: MRO class names (e.g. `["OptimizationError", "RuntimeError", "Exception", "BaseException"]`)
- `error_module`: defining module of the exception class
- `error_label`: the high-level CLI stage label

### Stage runners

| MCP tool | CLI subcmd | Purpose |
|---|---|---|
| `optimize_geometry` | `pdb2reaction opt` | Optimize a single molecular geometry |
| `find_transition_state` | `pdb2reaction tsopt` | TS search (RS-P-RFO / Dimer / TRIM / RS-I-RFO) |
| `run_irc` | `pdb2reaction irc` | IRC integration from a TS geometry |
| `compute_frequencies` | `pdb2reaction freq` | Vibrational analysis + thermochemistry |
| `run_single_point` | `pdb2reaction sp` | Single-point MLIP energy + forces (+optional Hessian) |

### Scans / paths / pipeline

| MCP tool | CLI subcmd | Purpose |
|---|---|---|
| `scan_1d` / `scan_2d` / `scan_3d` | `pdb2reaction scan` / `pdb2reaction scan2d` / `pdb2reaction scan3d` | Restraint-driven distance scans |
| `optimize_path` | `pdb2reaction path-opt` | Two-endpoint MEP optimization |
| `search_paths` | `pdb2reaction path-search` | Recursive reaction-pathway search |
| `run_full_pipeline` | `pdb2reaction` (the `all` subcmd) | End-to-end: extract → MEP → TS → IRC → freq → DFT |
| `run_single_point_dft` | `pdb2reaction dft` | Single-point DFT via gpu4pyscf |

### Structure / I/O helpers

| MCP tool | CLI subcmd | Purpose |
|---|---|---|
| `extract_active_site` | `pdb2reaction extract` | Cut a sphere around a ligand |
| `add_element_info` | `pdb2reaction add-elem-info` | Repair PDB element columns |
| `fix_altloc` | `pdb2reaction fix-altloc` | Resolve PDB alternate locations |
| `plot_trajectory` | `pdb2reaction trj2fig` | Energy profile figure (PNG default; also SVG/PDF/HTML/CSV) |
| `plot_energy_diagram` | `pdb2reaction energy-diagram` | Categorical energy diagram |
| `detect_bond_changes` | `pdb2reaction bond-summary` | Bond-change diff between two PDBs |

### Charge and ordered inputs

For tools exposing `charge` and `ligand_charge`, omit `charge` and supply a
per-resname `ligand_charge` mapping for PDB inputs when you want the total to be
derived. XYZ input without PDB residue context needs an explicit total `charge`.
A valid GJF charge/multiplicity header supplies both values automatically unless
you deliberately override them. If `charge` and `ligand_charge` are both
supplied, explicit `charge` wins.

`search_paths` requires `input_pdb` (reactant) and `product_pdb`; pass any
ordered intermediates as `intermediate_pdbs`. For staged 1D/full-pipeline scans,
put the first literal in `scan_lists` and later literals in
`additional_scan_stages`; the adapter emits one `--scan-lists` flag followed by
all stage values.

## Opt-in IRC controls

`run_irc` (CLI: `pdb2reaction irc`) accepts `irc_pos_def: bool` — when set,
IRC convergence additionally requires a positive-definite mass-weighted
Hessian, blocking the IRC "shoulder" false-convergence where the rms-only
criterion calls success before reaching the local minimum. Defaults to
`None` (rms-only, legacy).

`run_irc` also accepts `step_size` and `never_stop`. When a branch stops after
only a few frames, reduce `step_size` (typically to `0.05`) first. Set
`never_stop=True` to ignore gradient and energy endpoint criteria and trace to
the max-cycle cap; numerical/integration failures still stop. `run_full_pipeline`
forwards the same controls as `irc_step_size` / `irc_never_stop` and exposes
`flatten` plus `refine_path` for TS recovery.

`find_transition_state` (CLI: `pdb2reaction tsopt`) also exposes the
alternative TS optimizers via `--opt-mode`:

- `opt_mode="trim"` — Helgaker (1991) trust-region image-minimization TS opt
- `opt_mode="rsprfo"` — Banerjee (1985) restricted-step P-RFO TS opt

## Client configuration

Client configuration schemas differ. The snippet below applies to clients that
accept a top-level `mcpServers` object; check the client's own MCP
documentation before choosing the config file and schema.

- Claude Desktop — `~/Library/Application Support/Claude/claude_desktop_config.json` (macOS) / `%APPDATA%\Claude\claude_desktop_config.json` (Windows)
- Cursor — `~/.cursor/mcp.json`
- Other clients — consult the client's own MCP-server docs

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

See [`examples/mcp_client_config.json`](../examples/mcp_client_config.json)
for a full example with explicit env-var overrides (PATH / CUDA_VISIBLE_DEVICES).

VS Code instead uses a top-level `servers` object in `.vscode/mcp.json`
([VS Code MCP configuration reference](https://code.visualstudio.com/docs/agents/reference/mcp-configuration)):

```json
{
  "servers": {
    "pdb2reaction": {
      "command": "pdb2reaction-mcp",
      "args": []
    }
  }
}
```

### Custom Python MCP client

```python
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

server_params = StdioServerParameters(command="pdb2reaction-mcp")
async with stdio_client(server_params) as (read, write):
    async with ClientSession(read, write) as session:
        await session.initialize()
        tools = await session.list_tools()
        result = await session.call_tool(
            "optimize_geometry",
            arguments={
                "input_pdb": "r.pdb",
                "charge": -1,
                "max_cycles": 50,
            },
        )
        print(result.content)
```

## Sandbox / safety notes

- The MCP server inherits the calling environment's PATH, conda env, and
  CUDA setup. Long-running tools (opt / tsopt / irc) launch the
  `pdb2reaction` CLI in a subprocess — the agent should set `timeout_seconds`
  on each tool call to bound runaway computations.
- Stage/pipeline output files land under the `out_dir` kwarg (defaults to a unique
  `tempfile.mkdtemp(prefix="p2r_mcp_<subcmd>_…")` so concurrent agent calls
  don't collide). Helper tools use their explicit output path. Expert
  `extra_args` can request additional non-managed behavior; it cannot override
  typed output paths, `--out-dir`, or the managed `--out-json/--no-out-json`
  toggle. Inspect the returned `argv` when applying a filesystem policy.
- The server does not modify `~/.bashrc` / login env, install software,
  or authenticate model registries. All MLIP weights / PDB inputs must already
  exist on disk.
