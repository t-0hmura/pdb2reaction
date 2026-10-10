# Output Directory Layout

This page lists the main files that the commands write, the default output directories, and where each file goes inside `all`. Each command page lists all of its files under **Output files**.

## Filename conventions

| Filename | Written by | Purpose |
|---|---|---|
| `summary.json` | `all` and `path-search` | The aggregate JSON result (see [JSON Output Reference](json-output.md)). |
| `summary.json` | per-stage and report commands with `--out-json` (default `--no-out-json`) | A copy of the command's `result.json`; the two files are identical when the write completes. `fix-altloc`, `add-elem-info`, and `bond-summary` do not write it. |
| `result.json` | With `--out-json`: `opt`, `tsopt`, `freq`, `irc`, `sp`, `scan` / `scan2d` / `scan3d`, `path-opt`, `dft`, `extract`, `trj2fig`, `energy-diagram` | The JSON result of one command or report, also written when the run ends without converging. `extract`, `trj2fig`, and `energy-diagram` write it next to their first output file. |
| `run.log` | CLI and Colab runs, once the output directory exists | The command line (quoted for the shell) and the stdout/stderr printed while the command runs, including the lines that show how the run ended; the command's page lists them (for example [opt → Checking convergence](opt.md#checking-convergence)). Help, version, and dry-run calls do not create it, nor do the commands without an output directory (`extract`, `fix-altloc`, `add-elem-info`, `bond-summary`, `trj2fig`, `energy-diagram`). |
| `summary.log` | `path-search`, `all` | Text summary. The header gives `Scientific status`; the numbered sections give the barrier and bond changes of each segment and the output tree (see [all → Reading the run status](all.md#reading-the-run-status)). |
| `final_geometry.xyz` / `final_geometry.pdb` | `opt`, `tsopt` | Optimized geometry, the structure to pass to the next command; also written when the optimization does not converge. The `.xyz` is always written; the `.pdb` is for PDB/mmCIF input (`.gjf` for Gaussian input). |
| `mep_trj.pdb` / `mep_trj.cif` / `mep_trj.xyz` | `path-search` | Reaction path frames. The `.cif` is added for mmCIF or large-PDB input when file conversion is on. |
| `final_geometries_trj.xyz` / `hei.xyz` | `path-opt` | Reaction path frames (all images) and the highest-energy image, with `.pdb` / `.cif` / `.gjf` copies when the input format allows and file conversion is on. |
| `mep_plot.png` / `energy_diagram_MEP.png` | `path-search` | Energy profile of the minimum energy path (MEP): `mep_plot.png` along the path, `energy_diagram_MEP.png` as a state-energy diagram. Of the two, `all` places only `energy_diagram_MEP.png` at its root. |
| `finished_irc_trj.xyz` / `forward_irc_trj.xyz` / `backward_irc_trj.xyz` | `irc` | IRC (intrinsic reaction coordinate) trajectories (full path plus each branch), with a `.pdb` copy when a reference topology is available and a `.cif` copy for mmCIF or large-PDB input. |
| `finished_first.xyz` / `finished_last.xyz` | `irc` | Ends of the two branches: the endpoint candidates to optimize with [`opt`](opt.md). |
| `frequencies_cm-1.txt` | `freq` | Vibrational mode listing. |
| `*.pdb` / `*.cif` / `*.gjf` | commands with `--convert-files` (the default; `--no-convert-files` turns it off), and `extract` | Copies of the outputs in the input format, written next to them: PDB for PDB/mmCIF input, plus a `.cif` that keeps the original chain IDs, residue numbers, and insertion codes for mmCIF or large-PDB input; GJF for Gaussian input. `extract` has no toggle and always writes the `.cif` for mmCIF or large-PDB input. |

## Default `--out-dir`

| Subcommand | Default `--out-dir` |
|---|---|
| `all` | `./result_all/` |
| `opt` | `./result_opt/` |
| `tsopt` | `./result_tsopt/` |
| `freq` | `./result_freq/` |
| `irc` | `./result_irc/` |
| `dft` | `./result_dft/` |
| `scan` | `./result_scan/` |
| `scan2d` | `./result_scan2d/` |
| `scan3d` | `./result_scan3d/` |
| `path-opt` | `./result_path_opt/` |
| `path-search` | `./result_path_search/` |
| `sp` | `./result_sp/` |
| `extract` | `./` (writes `model.pdb`, or `model_<input>.pdb` for multiple inputs) |

Set another directory with `-o/--out-dir <path>`. `extract` takes one or more `-o/--output <file>` paths instead.

## Standalone vs `all`

A standalone subcommand writes a flat `result_<subcmd>/` directory, without `segments/` or `_work/`. Inside `all`, each post-processing stage uses the same file layout under `segments/seg_NN/` in `ts/`, `irc/`, `freq/`, and `dft/`.

- **`path-search` / `path-opt` are laid out differently.** Inside `all`, the MEP search runs `path-opt` by default and the recursive `path-search` with `--refine-path`; the raw output stays in `_work/path_opt/` or `_work/path_search/`, and only `mep_trj.*` and `energy_diagram_MEP.png` are placed at the root.

The tree below shows the main entries; [all → Output files](all.md#output-files) lists every file:

```text
result_all/
├─ summary.log · summary.json                 # run summary
├─ mep_trj.pdb · mep_trj.cif · mep_trj.xyz           # Core MEP coordinates
├─ mep_w_ref.{pdb,cif}                               # MEP merged into the full input (--write-ref-merge)
├─ energy_diagram_MEP.png · energy_diagram_*_all.png · irc_plot_all.png
├─ segments/
│  └─ seg_NN/                                  # One reaction step: seg_01, seg_02, ...
│     ├─ reactant.{pdb,cif,xyz,gjf} · ts.* · product.* # Optimized R, TS, and P (--tsopt)
│     └─ ts/ · irc/ · freq/{R,TS,P}/ · dft/         # per-stage working files (--tsopt / --thermo / --dft)
└─ _work/                                      # Intermediate files, including the TS candidates (HEI)
   ├─ models/ · scan/ · add_elem_info/ · fix_altloc/
   └─ path_opt/                                # MEP search and hei_seg_NN.* (path_search/ with --refine-path)
```

TS-only mode is a run on one TS candidate with `--tsopt` and no `-s/--scan-lists`. It has no MEP stage, so `_work/path_opt/` is absent and the deliverables live under `segments/seg_01/`.

## Notes

* **Runs that stop early**: a run that stops at the argument or input checks, before its output directory is set up, may write none of the files above.
* **`summary.json` / `result.json` without `--out-json`**: a successful per-stage run writes them only with `--out-json`. When a run stops on an exception after the output directory is set up, both are written even without the flag, with `"execution_status": "failed"` and an `"error_type"`.

## See Also

- [all](all.md#output-files) — the full `result_all/` tree and the energy diagrams
- [JSON Output Reference](json-output.md) — keys of `summary.json` and `result.json`, with Python and jq examples
- [Common options and selectors](cli-conventions.md) — `--out-dir`, `--convert-files`, and exit codes
- [Troubleshooting](troubleshooting.md) — common errors and fixes
