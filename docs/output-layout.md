# Output Directory Layout

Each `pdb2reaction` subcommand writes to its output directory under the filename conventions below, which agents and downstream scripts can rely on.

## Filename conventions

| Filename | Written by | Purpose |
|---|---|---|
| `summary.json` | `all`, `path-search`, and per-stage subcommands **only when `--out-json` is passed** (default `--no-out-json`) | Authoritative JSON envelope (see [JSON Output Reference](json-output.md)). Read this first. Pure utility subcommands (e.g. `fix-altloc`, `add-elem-info`, `bond-summary`) never emit it. |
| `result.json` | per-stage subcommands **only when `--out-json` is passed** — default `--no-out-json` (`opt`, `tsopt`, `freq`, `irc`, `sp`, `scan` / `scan2d` / `scan3d`, `path-opt`, `dft`, `extract`) | Alternate filename — identical payload to `summary.json`. Read `summary.json` for the single-filename convention; `result.json` carries the same content and can be deleted if you only consume `summary.json`. |
| `summary.log` | `path-search`, `all` | Human-readable run log (one row per segment / stage). |
| `final_geometry.xyz` | `opt`, `tsopt` | Optimized geometry (XYZ, full precision). |
| `mep.pdb` / `mep_trj.xyz` | `path-search` | Reaction path frames (PDB / XYZ). |
| `final_geometries_trj.xyz` / `hei.xyz` | `path-opt` | Standalone path-opt trajectory (full path) and highest-energy image (`.pdb` / `.gjf` companions when conversion is enabled). |
| `mep_plot.png` | `path-search` | Energy profile (PNG) of the MEP. (`all` promotes the styled `energy_diagram_MEP.png` to the root instead.) |
| `finished_irc_trj.xyz` / `forward_irc_trj.xyz` / `backward_irc_trj.xyz` | `irc` | IRC trajectories (full path plus per-branch; `.pdb` companions when a reference PDB is available). |
| `frequencies_cm-1.txt` | `freq` | Vibrational mode listing. |
| `*.gjf` | various (when `--convert-files`) | Gaussian-format companion structure. |

## Default `--out-dir`

| Subcommand | Default `--out-dir` |
|---|---|
| `all` | `./result_all/` |
| `opt` | `./result_opt/` |
| `tsopt` | `./result_tsopt/` |
| `freq` | `./result_freq/` |
| `irc` | `./result_irc/` |
| `dft` | `./result_dft/` |
| `scan` / `scan2d` / `scan3d` | `./result_scan*/` |
| `path-opt` / `path-search` | `./result_path_*/` |
| `sp` | `./result_sp/` |
| `extract` | `./` (writes `model.pdb`, or `model_<input>.pdb` for multiple inputs) |

Override with `--out-dir <path>` (or `-o`).

## Standalone vs `all`

A subcommand run on its own writes a **flat** result directory. The same writer, when orchestrated by `all`, nests into a structured tree. The two layouts differ by design:

- **Standalone subcommand** → flat `result_<subcmd>/` with the files above. There is no `segments/` and no `_work/`; those only appear when `all` coordinates several writers in one run.
- **Inside `all`, leaf writers nest unchanged.** A per-segment leaf output at `segments/seg_NN/<subcmd>/` is structurally identical to the standalone `result_<subcmd>/` — `all` simply hands the same writer a different output directory.
- **`path-search` / `path-opt` are the engine exception.** Run standalone, `path-search` is itself a deliverable (`result_path_search/` with its own `summary.log`, `mep.pdb`, `mep_trj.xyz`, `mep_plot.png`, `energy_diagram_MEP.png`). Inside `all`, its raw output is treated as scratch under `_work/path_opt/` (or `_work/path_search/` with `--refine-path`), and only the merged products (`mep.pdb`, `mep_trj.xyz`, `mep_w_ref.pdb`, `energy_diagram_MEP.png`) are promoted to the pipeline root.

The `all` tree therefore has three zones:

```text
result_all/
├─ summary.log · summary.json                 # authored at the root
├─ mep.pdb · mep_w_ref.pdb · mep_trj.xyz       # MEP deliverables promoted from the engine
├─ energy_diagram_MEP.png · energy_diagram_*.png
├─ segments/
│  └─ seg_NN/                                  # 2-digit per-reactive-segment deliverables
│     ├─ reactant.{pdb,xyz,gjf} · ts.* · product.*   # canonical R/TS/P
│     └─ ts/ · irc/ · freq/{R,TS,P}/ · dft/         # per-stage working files (--tsopt / --thermo / --dft)
└─ _work/                                      # pipeline scratch (safe to rm -rf)
   ├─ models/ · scan/ · add_elem_info/ · fix_altloc/
   └─ path_opt/                                # raw MEP-engine output (path_search/ with --refine-path True)
```

In TSOPT-only mode there is no MEP stage, so `_work/path_opt/` is absent and the deliverables live under `segments/seg_01/`. See [all](all.md) for the full per-mode breakdown.

## Agent recipe

```python
# Read whichever subcommand's output, single filename across the board.
import json
from pathlib import Path

summary = json.loads((Path(out_dir) / "summary.json").read_text())

if summary["status"] == "error":
    chain = summary.get("error_class_chain", [])
    if "OptimizationError" in chain:
        # retry with looser convergence threshold
        ...
    else:
        raise RuntimeError(summary["error"])
```

`summary.json` / `result.json` are written by `all` and `path-search`, and by per-stage subcommands **only when `--out-json` is passed** (default `--no-out-json`). When written, the envelope carries the schema version + status (and, on the failure path, the error class chain); do not assume a per-stage `summary.json` exists by default.
