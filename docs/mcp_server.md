# pdb2reaction MCP server

This page explains how to call the 18 pdb2reaction tools from an AI agent over
MCP (Model Context Protocol): installation, the list of tools, and client
configuration. The server, `pdb2reaction-mcp` (alias `p2r-mcp`), speaks
JSON-RPC over stdio, so any [MCP](https://modelcontextprotocol.io/) client can
use it, including Claude Desktop, Claude Code, Cursor, Codeium, and agents built
on the official Python or TypeScript MCP SDKs.

## Install

```bash
pip install "pdb2reaction[mcp]"
```

This adds the `mcp[cli]` dependency and the `pdb2reaction-mcp` / `p2r-mcp`
commands.

## Tools

18 tools, one per CLI subcommand. Each tool returns a structured dict with:

- `schema_version`: version of the result format
- `execution_status`: `completed` | `failed`
- `scientific_status`: `success` | `partial` | `failed`
- `summary_status`: only `ok` comes with a non-empty `summary`
  - `ok`: the `summary` of this call was read
  - `not_required`: a structure / I/O helper, which writes no summary
  - `summary_missing`: no `summary.json` in `out_dir`
  - `summary_parse_error`: `summary.json` cannot be read as a JSON object
  - `summary_run_mismatch`: the file belongs to another run
- `exit_code`: exit code of the CLI process
- `out_dir`: output directory of the stage runners and the scan / path / pipeline tools, or null for the structure / I/O helpers
- `summary`: parsed `summary.json`; an empty object for the structure / I/O helpers
- `stderr_tail` / `stdout_tail`: last ~60 lines of process output
- `hint`: the recovery hint from a `; recover: <hint>` suffix of a CLI error message, if any
- `argv`: the full command line that was run
- `run_id`: UUID of this call

Each tool runs the CLI in a subprocess that keeps the caller's working
directory, so relative input paths keep their meaning. The tables below list
each tool's required arguments; the client receives every argument with its
type in the tool's input schema.

- Input paths such as `input_pdb`, `ts_pdb`, and `reactant_pdb` go to the command's `-i`, so they take any format that `-i` reads, such as XYZ.
- Optional arguments set the command's CLI options, for example `charge` (`-q`) and `max_cycles` (`--max-cycles`).

### Structured error envelope

When a stage runner or a scan / path / pipeline tool fails, its `summary` carries the error fields below, so agents can match the error class without parsing text. The structure / I/O helpers, and a run that stops before it writes `summary.json` (`summary_missing`), return no error fields; read `stderr_tail` and `hint` instead.

- `error`: the error message
- `error_type`: exception class name
- `error_class_chain`: the class and its parent classes, most specific first
- `error_module`: module that defines the exception class
- `error_label`: the high-level CLI stage label

### Stage runners

| MCP tool | Required arguments | CLI subcmd | Purpose |
|---|---|---|---|
| `optimize_geometry` | `input_pdb` | `pdb2reaction opt` | Optimize a single molecular geometry |
| `find_transition_state` | `ts_pdb` | `pdb2reaction tsopt` | TS search (RS-P-RFO / Dimer / TRIM / RS-I-RFO) |
| `run_irc` | `ts_pdb` | `pdb2reaction irc` | IRC integration from a TS geometry |
| `compute_frequencies` | `input_pdb` | `pdb2reaction freq` | Vibrational analysis + thermochemistry |
| `run_single_point` | `input_pdb` | `pdb2reaction sp` | Single-point energy + forces with the chosen `backend` (+optional Hessian) |

### Scans / paths / pipeline

| MCP tool | Required arguments | CLI subcmd | Purpose |
|---|---|---|---|
| `scan_1d` / `scan_2d` / `scan_3d` | `input_pdb`, `scan_lists` | `pdb2reaction scan` / `pdb2reaction scan2d` / `pdb2reaction scan3d` | Restraint-driven scans of distances, angles, or dihedrals |
| `optimize_path` | `reactant_pdb`, `product_pdb` | `pdb2reaction path-opt` | Two-endpoint MEP optimization |
| `search_paths` | `input_pdb`, `product_pdb` | `pdb2reaction path-search` | Recursive reaction-pathway search |
| `run_full_pipeline` | `reactant_pdb` | `pdb2reaction all` | End-to-end: extract → MEP → TS → IRC → freq → DFT |
| `run_single_point_dft` | `input_pdb` | `pdb2reaction dft` | Single-point DFT energy and atomic charges (GPU4PySCF or PySCF) |

### Structure / I/O helpers

| MCP tool | Required arguments | CLI subcmd | Purpose |
|---|---|---|---|
| `extract_active_site` | `complex_pdb`, `ligand_id`, `radius_angstrom`, `output_pdb` | `pdb2reaction extract` | Active-site model: residues near the ligand, capped with hydrogens |
| `add_element_info` | `input_pdb`, `output_pdb` | `pdb2reaction add-elem-info` | Repair PDB element columns |
| `fix_altloc` | `input_pdb`, `output_pdb` | `pdb2reaction fix-altloc` | Resolve PDB alternate locations |
| `plot_trajectory` | `input_trj_xyz`, `output_png` | `pdb2reaction trj2fig` | Energy profile figure (PNG default; also JPEG/SVG/PDF/HTML/CSV) |
| `plot_energy_diagram` | `energies`, `output_png` | `pdb2reaction energy-diagram` | State energy diagram from given energies |
| `detect_bond_changes` | `reactant_pdb`, `product_pdb` | `pdb2reaction bond-summary` | Bond changes between two structures (XYZ / PDB / mmCIF / GJF) |

### Charge and ordered inputs

For tools that take `charge` and `ligand_charge`, omit `charge` and give a
per-resname `ligand_charge` mapping for PDB inputs when you want the total
charge derived from residue names. XYZ input without PDB residue information
needs an explicit total `charge`. A valid GJF charge/multiplicity header supplies
both values unless you override them. If both `charge` and `ligand_charge` are
given, the explicit `charge` wins.

`search_paths` takes any ordered intermediates between `input_pdb` (reactant)
and `product_pdb` as `intermediate_pdbs`. For staged scans in `scan_1d` and
`run_full_pipeline`, put the first stage in `scan_lists` and later stages in
`additional_scan_stages`; the tool passes one `--scan-lists` flag to the CLI,
followed by all stage values.

## IRC and TS optimization settings

The IRC and TS arguments are the CLI options of the same name; the command
pages give their meaning and defaults.

- `run_irc` — `step_size`, `never_stop`, `irc_pos_def`: [`irc`](irc.md); `--irc-pos-def` is in the [generated reference](reference/commands/irc.md)
- `find_transition_state` — `opt_mode`: [`tsopt`](tsopt.md) `--opt-mode`, default `hess` (RS-P-RFO); see {ref}`--opt-mode by command <opt-mode-semantics>`
- `run_full_pipeline` — `irc_step_size`, `irc_never_stop`, `flatten`, `refine_path`: [`all`](all.md); `--irc-step-size` and `--irc-never-stop` are in the [generated reference](reference/commands/all.md)

## Client configuration

Client configuration schemas differ. The snippet below applies to clients that
accept a top-level `mcpServers` object.

- Claude Desktop — `~/Library/Application Support/Claude/claude_desktop_config.json` (macOS) / `%APPDATA%\Claude\claude_desktop_config.json` (Windows)
- Cursor — `~/.cursor/mcp.json`
- Claude Code — no file to edit: run `claude mcp add pdb2reaction -- pdb2reaction-mcp`; `claude mcp list` then shows `✔ Connected` for it
- Other clients — consult the client's own MCP-server docs

Once the client has started the server, its tool list shows the 18 tools.

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

See [`examples/mcp_client_config.json`](https://github.com/t-0hmura/pdb2reaction/blob/main/examples/mcp_client_config.json)
for a full example that sets environment variables (PATH / CUDA_VISIBLE_DEVICES).

VS Code instead uses a [top-level `servers` object](https://code.visualstudio.com/docs/agents/reference/mcp-configuration)
in `.vscode/mcp.json`:

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
import asyncio

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

async def main():
    server_params = StdioServerParameters(command="pdb2reaction-mcp")
    async with stdio_client(server_params) as (read, write):
        async with ClientSession(read, write) as session:
            await session.initialize()
            result = await session.call_tool(
                "optimize_geometry",
                arguments={
                    "input_pdb": "r.pdb",
                    "charge": -1,
                    "max_cycles": 50,
                },
            )
            print(result.content)

asyncio.run(main())
```

## Sandbox / safety notes

- The server inherits the caller's PATH, conda environment, and CUDA setup.
  Set `timeout_seconds` on long calls (opt / tsopt / irc) to stop runaway
  calculations (default: no timeout).
- The stage runners and the scan / path / pipeline tools write under `out_dir`. When it is not given, each
  call gets its own temporary directory (`p2r_mcp_<subcmd>_…`), so parallel
  calls do not collide.
- The structure / I/O helpers have no `out_dir` and write to their explicit output path.
  `extra_args` passes additional CLI flags but cannot override the typed output
  paths, `--out-dir`, or `--out-json/--no-out-json`; the returned `argv` shows
  every path the command was given.
- The server does not change `~/.bashrc` or the login environment, install
  software, or log in to model registries. MLIP weights and input PDB files must
  already be on disk.

## See also

* [JSON Output Reference](json-output.md) — the status fields and `summary.json` that the tools return
* [Troubleshooting](troubleshooting.md) — what to do when a run fails
* [Command Reference](reference/commands/index.md) — the CLI options behind each tool, for `extra_args`
