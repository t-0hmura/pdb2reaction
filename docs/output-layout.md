# Output Directory Layout

Each `pdb2reaction` subcommand writes to its output directory under the filename conventions below, which agents and downstream scripts can rely on.

## Filename conventions

| Filename | Written by | Purpose |
|---|---|---|
| `summary.json` | `all` and `path-search` after their summary writer is reached | Authoritative aggregate JSON envelope (see [JSON Output Reference](json-output.md)). Early CLI/input validation may fail before it exists. |
| `summary.json` | successful per-stage/report runs when `--out-json` is passed (default `--no-out-json`); caught runtime errors may write a best-effort error envelope without the flag | Compatibility mirror of the leaf `result.json`. A successful writer return guarantees identical bytes; pure utilities such as `fix-altloc`, `add-elem-info`, and `bond-summary` do not use it. |
| `result.json` | same conditions as the per-stage `summary.json` (`opt`, `tsopt`, `freq`, `irc`, `sp`, scan variants, `path-opt`, `dft`, `extract`, `trj2fig`, `energy-diagram`) | Authoritative leaf/report envelope. It is published after the compatibility mirror; consume this file when distinguishing interrupted generations. |
| `summary.log` | `path-search`, `all` | Human-readable run log (one row per segment / stage). |
| `final_geometry.xyz` | `opt`, `tsopt` | Optimized geometry (XYZ, full precision). |
| `mep.pdb` / `mep.cif` / `mep_trj.xyz` | `path-search` | Reaction path frames; `.cif` is added for mmCIF/oversized-PDB topology when conversion is enabled. |
| `final_geometries_trj.xyz` / `hei.xyz` | `path-opt` | Standalone path-opt trajectory (full path) and highest-energy image (`.pdb` / `.cif` / `.gjf` companions when applicable and conversion is enabled). |
| `mep_plot.png` | `path-search` | Energy profile (PNG) of the MEP. (`all` promotes the styled `energy_diagram_MEP.png` to the root instead.) |
| `finished_irc_trj.xyz` / `forward_irc_trj.xyz` / `backward_irc_trj.xyz` | `irc` | IRC trajectories (full path plus per-branch; `.pdb` and, for bridge topology, `.cif` companions when a reference topology is available). |
| `frequencies_cm-1.txt` | `freq` | Vibrational mode listing. |
| `*.cif` / `*.gjf` | various (when `--convert-files`) | Identifier-preserving mmCIF or Gaussian-format companion structure, according to the input template. |

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
- **`path-search` / `path-opt` are the engine exception.** Run standalone, each is itself a deliverable: `path-search` → `result_path_search/` (`summary.log`, `mep.pdb`, optional `mep.cif`, `mep_trj.xyz`, `mep_plot.png`, `energy_diagram_MEP.png`), and `path-opt` → `result_path_opt/` (`final_geometries_trj.xyz`, `hei.xyz`). Inside `all`, the raw engine output is treated as scratch under `_work/path_opt/` (or `_work/path_search/` with `--refine-path`), and only the merged products (`mep.pdb`, optional `mep.cif`, `mep_trj.xyz`, `mep_w_ref.pdb` / `.cif`, `energy_diagram_MEP.png`) are promoted to the pipeline root.

The `all` tree therefore has three zones:

```text
result_all/
├─ summary.log · summary.json                 # authored at the root
├─ mep.{pdb,cif} · mep_w_ref.{pdb,cif} · mep_trj.xyz # CIF for bridge input
├─ energy_diagram_MEP.png · energy_diagram_*.png
├─ segments/
│  └─ seg_NN/                                  # 2-digit per-reactive-segment deliverables
│     ├─ reactant.{pdb,cif,xyz,gjf} · ts.* · product.* # canonical R/TS/P
│     └─ ts/ · irc/ · freq/{R,TS,P}/ · dft/         # per-stage working files (--tsopt / --thermo / --dft)
└─ _work/                                      # pipeline scratch (safe to rm -rf)
   ├─ models/ · scan/ · add_elem_info/ · fix_altloc/
   └─ path_opt/                                # raw MEP-engine output (path_search/ with --refine-path)
```

In TSOPT-only mode there is no MEP stage, so `_work/path_opt/` is absent and the deliverables live under `segments/seg_01/`. See [all](all.md) for the full per-mode breakdown.

## Agent recipe

```python
# Select the authoritative name for the command that produced out_dir.
import json
from pathlib import Path

subcommand = "opt"  # replace with the command you ran
primary = "summary.json" if subcommand in {"all", "path-search"} else "result.json"
summary = json.loads((Path(out_dir) / primary).read_text())

if summary["status"] == "error":
    error_type = summary.get("error_type", "RuntimeError")
    raise RuntimeError(f"{error_type}: {summary['error']}")
```

`summary.json` / `result.json` are written on successful per-stage runs only
when `--out-json` is passed (default `--no-out-json`). A caught runtime error
may emit a best-effort error envelope even without the flag; validation exits
or failures before output setup may emit none. When written, the envelope
carries schema version + status (and, on the error path, the class chain).
