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
| `optimize_geometry` | `pdb2reaction opt` | Optimize one structure |
| `find_transition_state` | `pdb2reaction tsopt` | TS search (RS-P-RFO / Dimer / TRIM / RS-I-RFO) |
| `run_irc` | `pdb2reaction irc` | IRC integration from a TS geometry |
| `compute_frequencies` | `pdb2reaction freq` | Vibrational analysis + thermochemistry |
| `run_single_point` | `pdb2reaction sp` | Single-point MLIP energy + forces (+optional Hessian) |

### Scans / paths / pipeline

| MCP tool | wraps | purpose |
|---|---|---|
| `scan_1d` / `scan_2d` / `scan_3d` | `pdb2reaction scan` / `scan2d` / `scan3d` | Restraint-driven distance scans |
| `optimize_path` | `pdb2reaction path-opt` | Two-endpoint MEP optimization |
| `search_paths` | `pdb2reaction path-search` | Recursive reaction-pathway search |
| `run_full_pipeline` | `pdb2reaction all` | Configurable end-to-end pipeline; TS/IRC, thermo/freq, and DFT stages run only when their tool arguments enable them |
| `run_single_point_dft` | `pdb2reaction dft` | Single-point DFT via GPU4PySCF by default; CPU PySCF is available through `extra_args=["--engine", "cpu"]` |

### Structure / I/O helpers

| MCP tool | wraps | purpose |
|---|---|---|
| `extract_active_site` | `pdb2reaction extract` | Cut a sphere around a ligand |
| `add_element_info` | `pdb2reaction add-elem-info` | Repair PDB element columns |
| `fix_altloc` | `pdb2reaction fix-altloc` | Resolve PDB alternate locations |
| `plot_trajectory` | `pdb2reaction trj2fig` | Energy profile PNG / SVG / PDF / HTML / CSV |
| `plot_energy_diagram` | `pdb2reaction energy-diagram` | Categorical energy diagram |
| `detect_bond_changes` | `pdb2reaction bond-summary` | Bond-change diff between two structures, including PDB/mmCIF |

## `SubcmdResult` return schema

Every tool returns the same structured dict so the calling agent can dispatch on `status` without parsing stderr:

```python
{
    "schema_version": "1.1",         # pin to MCP_SUBCMD_RESULT_SCHEMA_VERSION
    "status": "ok" | "failed" | "summary_missing" | "summary_parse_error" | "summary_run_mismatch",
    "exit_code": int,                # subprocess exit code
    "out_dir": str | None,           # working directory the CLI wrote to (None for the I/O helpers)
    "summary": dict,                 # parsed summary.json; {} for the I/O helpers
    "stderr_tail": str,              # last ~60 lines of stderr
    "stdout_tail": str,              # last ~60 lines of stdout
    "hint": str | None,              # parsed `; recover: <hint>` suffix
    "argv": list[str],               # full argv that was executed (reproducible)
    "run_id": str,                   # invocation ownership token
}
```

`summary_run_mismatch` means the output directory contains a summary owned by
a different invocation; clients must not treat that stale summary as the
current result. `run_id` identifies the invocation used for this ownership
check.

`summary` and `out_dir` carry data for the 12 tools that run their subcommand against an output directory (the stage runners, the scans, `optimize_path`, `search_paths`, `run_full_pipeline`, `run_single_point_dft`). The 6 structure / I/O helper tools (`extract_active_site`, `add_element_info`, `fix_altloc`, `plot_trajectory`, `plot_energy_diagram`, `detect_bond_changes`) run with `out_dir=None` and write no summary.json: they always return `summary={}` and `out_dir=None`, with `status` taken straight from the exit code (`ok` / `failed`). Read their result from `stdout_tail`, and their failure reason from `stderr_tail` / `exit_code`.

For those 12 summary-writing tools, a failed subcommand additionally surfaces a structured exception envelope inside `summary`:

- `error_class_chain` — list of class names walking the MRO (e.g. `["CudaOutOfMemoryError", "RuntimeError", "Exception"]`)
- `error_module` — module path of the originating exception class
- `error_label` — high-level CLI stage label (e.g. `"IRC"`, `"optimization"`, `"frequency analysis"`, `"path optimization"`, fallback `"UnhandledError"`)

so an MCP client can pattern-match on the hierarchy instead of substring-matching `stderr_tail`.

## Charge and ordered-input contract

For tools that accept `charge` and `ligand_charge`, use the same precedence as
the CLI:

- PDB/mmCIF input: normally omit `charge` and set `ligand_charge="RES:Q,..."`; the
  total is derived from all selected amino acids, ions, and unknown ligands.
- XYZ without residue context: provide the explicit total `charge`. A valid
  GJF header already supplies charge/multiplicity unless deliberately
  overridden.
- Outside extraction, an explicit `charge` overrides inferred charge. For
  `run_full_pipeline` when extraction is active, it is instead an assertion
  against the extract-derived total; a mismatch aborts rather than silently
  overriding the extracted chemistry.

`search_paths` requires both `input_pdb` (reactant) and `product_pdb`; optional
`intermediate_pdbs=[...]` are inserted between them in the given reaction
order. Never try to invoke recursive path search with one structure.

For staged scans, `scan_1d` and `run_full_pipeline` take the first literal in
`scan_lists` and later literals in `additional_scan_stages`. They generate one
`--scan-lists` occurrence followed by all values, because repeating that CLI
flag is rejected.

## `find_transition_state` opt-mode kwarg

`find_transition_state` forwards `opt_mode` to `pdb2reaction tsopt --opt-mode`, which accepts six strings (four modes plus two aliases; default `"hess"`):

- `opt_mode="hess"` (default; alias `"rsprfo"`) — Banerjee restricted-step P-RFO TS opt
- `opt_mode="grad"` (alias `"dimer"`) — Hessian-guided Dimer TS opt
- `opt_mode="trim"` — Helgaker trust-region image-minimization TS opt
- `opt_mode="rsirfo"` — restricted-step I-RFO TS opt

## Sandbox / safety notes

IRC recovery knobs are exposed directly: `run_irc(step_size=0.05,
never_stop=True)` and `run_full_pipeline(irc_step_size=0.05,
irc_never_stop=True)`. Reduce the step first for an immediately stopping IRC;
`never_stop` ignores gradient and energy endpoint criteria and traces to the
cycle cap unless propagation fails. The full pipeline also exposes `flatten` and `refine_path`; recursive
refinement can split a poor path into extra segments and increase cost.

- The MCP server inherits the calling environment's PATH, conda env, and CUDA setup. Long-running tools (opt / tsopt / irc) launch the `pdb2reaction` CLI in a subprocess — set `timeout_seconds` on each call to bound runaway computations.
- The 12 stage / scan / pipeline tools take an `out_dir` kwarg and write their output (including `summary.json`) there; left unset, it defaults to a unique `tempfile.mkdtemp("p2r_mcp_<subcmd>_…")` so concurrent agent calls don't collide.
- The 6 structure / I-O helpers (`extract_active_site`, `add_element_info`, `fix_altloc`, `plot_trajectory`, `plot_energy_diagram`, `detect_bond_changes`) take no `out_dir`: each writes to the path given by its required `output_pdb` / `output_png` argument (`detect_bond_changes` writes no file). They return `status: "ok"` with `out_dir: null` and an empty `summary` — read their result from `stdout_tail` and the file they wrote.
- The server leaves `~/.bashrc` / login env untouched and installs no software.
  Normal named parameters write to `out_dir` or the explicit helper output
  path. `extra_args` is an expert escape hatch and can request additional CLI
  outputs or in-place behavior, so inspect the returned `argv` before treating
  a call as sandboxed. Structure inputs must exist locally. Backend weights may
  download on first load when network access/authentication permit it; pre-cache
  them for offline or sandboxed operation.

## See also
- Full MCP server doc: [`docs/mcp_server.md`](../../docs/mcp_server.md)
- Sample MCP-client config: [`examples/mcp_client_config.json`](../../examples/mcp_client_config.json)
- Server / tool source: `pdb2reaction/mcp/{server,_tools,_runner}.py`
