# `pdb2reaction all` — base orientation

## Purpose

`all` is the meta-command that orchestrates the stages selected by its input
mode and flags: optional extraction, scan/MEP or TS-only entry, and optional
TSOPT/IRC, frequency, and DFT post-processing. The MEP stage runs single-pass
`path-opt` by default; pass `--refine-path` to run recursive `path-search`
instead. `all`
resolves three input modes via flag context (see the three companion
mds: `all-endpoint-mep.md`, `all-scan-list.md`, `all-ts-only.md`).

Use `all` when you want a **single qsub-able invocation** that
produces R / TS / P / IM coordinates plus barrier candidates for one or
more path segments. Frequency and IRC checks determine which candidates are
validated elementary steps.

## Synopsis

```bash
pdb2reaction all -i <input(s)> [-c <substrate>] [-l 'RES:Q,...'] \
    [--scan-lists '...'] [--tsopt] [--thermo] [--dft] \
    [-b uma|orb|mace|aimnet2] [-o result_all/]
```

## Key flags (cross-mode)

> **Note:** Do not pass `--max-cycles` to `pdb2reaction all`; let each stage
> use its own default. Set it only when running a single-stage subcommand
> directly, such as `opt`, `tsopt`, or `path-opt`.

| Flag | Type | Default | Description |
|---|---|---|---|
| `-i, --input` | path(s) | required | One or more reaction-ordered structures, or a TS-candidate alone |
| `-c, --center` | str | (uses input as-is) | Residue names, IDs, PDB/mmCIF path, chain-qualified name (`'B:SAM'`), or exact named ID (`'B:SAM:321'`) |
| `-l, --ligand-charge` | str | none | Per-resname charges, e.g. `'SAM:1,GPP:-3'` (or a bare number = total). With `-c` it feeds extraction; without `-c` it derives the total from the full PDB/mmCIF model. |
| `-q, --charge` | int | derived from `-l` | Total system charge. With `-c` (extraction runs) `-q` acts as an **assertion**: it must **match** the extract-derived charge, otherwise the run aborts with `BadParameter` (this is also checked by `--dry-run`). With `-c` omitted (pre-carved input passed as-is) `-q` sets the total directly and **emits a warning**. Normally unnecessary — let `-l` derive the total. |
| `-m, --multiplicity` | int | 1 | Spin multiplicity (2S+1) |
| `-r, --radius` | float | 2.6 | Pocket radius (Å) when `-c` triggers extraction |
| `-s, --scan-lists` | one flag followed by one or more values | none | Staged distance scans (mode 2 — `all-scan-list.md`). Use one `-s` occurrence; each following Python literal is one sequential stage. Repeating the flag is rejected. |
| `--refine-path / --no-refine-path` | toggle | off | Recursive `path-search` when enabled; single-pass `path-opt` when disabled. Refinement can improve a poor TS seed but may split a bad path into unnecessary segments and greatly increase cost |
| `--thresh` | str | `gau` | Convergence preset for single-structure optimization and scan relaxation |
| `--thresh-gsm` | str | `gau_loose` | Convergence preset for the GSM string optimizer |
| `--thresh-dmf` | str/float | `tight` | DMF IPOPT dual-infeasibility tolerance: `tight`, `middle`, `loose`, or a positive float |
| `--tsopt / --no-tsopt` | toggle | off | Run TS optimization + IRC per reactive segment (also required to enter TS-only mode with a single `-i`) |
| `--flatten/--no-flatten` | flag | off | Enable surplus-imaginary-mode cleanup when TSOPT does not reach a first-order saddle |
| `--irc-step-size` | float | IRC default `0.10` | Forward a smaller EulerPC maximum step; try `0.05` when an IRC branch stops after only a few frames |
| `--irc-never-stop/--no-irc-never-stop` | flag | off | Ignore IRC gradient and energy stops and trace to the cycle cap; propagation failures still stop |
| `--reject-uphill / --no-reject-uphill` | toggle | off | Opt in to rejection above `1e-3` Hartree during Hessian/RFO post-IRC endpoint re-optimization only. At the emergency floor, the retained endpoint receives a final convergence check. It never affects TS optimization or path search. |
| `--thermo / --no-thermo` | toggle | off | Run freq + thermochemistry on R / TS / P |
| `--freq-symmetry-number` | int ≥ 1 | child YAML/default (normally 1) | Use one external rotational symmetry number for every R/TS/P frequency job. Point-group symmetry is not inferred. |
| `--dft / --no-dft` | toggle | off | Run DFT single point on R / TS / P |
| `--dft-func-basis` | str | `wb97m-v/def2-tzvpd` | DFT functional/basis (when `--dft`) |
| `-b, --backend` | str | `uma` | MLIP backend |
| `--workers`, `--workers-per-node` | int | `1`, `1` | UMA predictor workers. `workers > 1` plus an explicit `Analytical` Hessian raises `BackendError`; use one worker or `FiniteDifference`. Other built-in backends ignore these worker kwargs. |
| `--solvent` | str | none | xTB-ALPB solvent name (`water`, `methanol`, …) |
| `-o, --out-dir` | path | `./result_all/` | Top-level output directory |
| `--config` | path | none | YAML config applied before CLI flags |
| `--show-config` | flag | off | Print the resolved config and continue running |
| `--dry-run` | flag | off | Validate inputs and print the plan, then exit before scan/MEP/TSOPT/IRC/freq/DFT. With `-c/--center`, runs extraction in a temporary directory to validate the derived charge and electron parity, then deletes it. |
| `--help-advanced` | flag | — | Reveal hidden flags (freeze, advanced overrides) |

Run `pdb2reaction all --help-advanced` for the full list (it changes
between versions).

## Mode selection cheatsheet

| Inputs | Mode md | Behavior |
|---|---|---|
| Single `-i input.{xyz,pdb,cif,gjf}` + `--tsopt` (no `--scan-lists`) | `all-ts-only.md` | treat input as TS candidate; tsopt + IRC, plus R/TS/P freq only with `--thermo` |
| Single `-i input.{pdb,cif,xyz,gjf}` + `--scan-lists '...'` | `all-scan-list.md` | single reactant + staged distance scans; XYZ/GJF need numeric selectors and explicit/header charge because residue selection/extraction is unavailable |
| Multiple reaction-ordered PDB/mmCIF/XYZ/GJF inputs | `all-endpoint-mep.md` | multi-endpoint MEP |

A single `-i` without **either** `--scan-lists` or `--tsopt` raises
`BadParameter` (see `all.py`: "Provide at least two structures... or a
single structure with --scan-lists, or a single structure with --tsopt").

## Output tree (typical)

Three zones: deliverables at `<out_dir>/`, per-segment deliverables under `<out_dir>/segments/seg_NN/`, scratch under `<out_dir>/_work/`. `<work_path>` = `_work/path_opt/` (default) or `_work/path_search/` (`--refine-path`).

| Path | When | Content |
|---|---|---|
| `<out_dir>/summary.json` | pipeline reaches its summary writer | machine-readable per-stage results; early CLI/input validation can fail before this file exists |
| `<out_dir>/summary.log` | pipeline reaches its summary writer | human-readable text + dir tree; early CLI/input validation can fail before this file exists |
| `<out_dir>/mep_trj.xyz`; `mep.pdb`; bridge inputs also `mep.cif` | successful MEP/scan-list mode; companions additionally require `--convert-files` and topology | stitched MEP across segments |
| `<out_dir>/mep_w_ref.pdb`; bridge inputs also `.cif` | requested recursive full-template merge succeeds; CIF additionally needs bridge metadata | MEP merged into the full-system template |
| `<out_dir>/energy_diagram_MEP.png` | MEP/scan-list mode when diagram export succeeds | bare all-segment MEP energies |
| `<out_dir>/energy_diagram_{MLIP,G_MLIP,DFT,G_DFT_plus_MLIP}_all.png` | matching stages provide finite energies and PNG export succeeds | aggregated multi-segment diagrams |
| `<out_dir>/segments/seg_NN/{reactant,ts,product}.xyz`; topology companions `.pdb`/`.cif` | successful `--tsopt` + IRC/endpoint processing; companions require `--convert-files` and topology/bridge metadata | canonical post-IRC/endpoint-optimized R/TS/P (2-digit) |
| `<out_dir>/segments/seg_NN/ts/final_geometry.xyz`, `vib/imag_*_trj.xyz`; topology companions | TS optimizer reaches output; companions require `--convert-files` and topology | final TS attempt and imaginary-mode displacement |
| `<out_dir>/segments/seg_NN/irc/{forward,backward,finished}_irc_trj.xyz` | corresponding IRC branch/stitching reaches output | IRC trajectories |
| `<out_dir>/segments/seg_NN/freq/{R,TS,P}/{frequencies_cm-1.txt, thermoanalysis.yaml}` | `--thermo` and the corresponding frequency stage succeeds | per-state freq + thermo |
| `<out_dir>/segments/seg_NN/dft/{R,TS,P}/result.yaml` | `--dft` and the corresponding state calculation reaches output | per-state DFT; the in-pipeline subprocess does not request standalone `result.json` |
| `<out_dir>/segments/seg_NN/energy_diagram_{MLIP,G_MLIP,DFT,G_DFT_plus_MLIP}.png` | matching stages provide finite energies and PNG export succeeds | per-segment diagrams |
| `<out_dir>/_work/models/model_<stem>.pdb` and, for bridge inputs, `.cif` | `-c` given | extracted active-site clusters |
| `<work_path>/seg_NNN_<tag>/` | corresponding MEP string is attempted | per-string MEP scratch (3-digit, `_mep` / `_maxdepth` / `_bridge`) |
| `<work_path>/mep_seg_NN_trj.xyz`, `mep_seg_NN.{pdb,cif,gjf}` | corresponding segment output succeeds; companions require `--convert-files` and topology/template | per-segment MEP frames |
| `<work_path>/hei_seg_NN.{xyz,pdb,cif,gjf}` | a TS-candidate HEI is emitted; companions require `--convert-files` and topology/template | HEI candidate per segment (TS seed) |

## Output keys (summary.json — top level)

```python
import json
d = json.load(open("result_all/summary.json"))
print(d["status"])                    # "success" / "partial" / "failed"
print(d["pdb2reaction_version"])
print(d["charge"], d["spin"])
print(d["rate_limiting_step"])        # legacy key: highest local segment barrier
print(len(d["segments"]))             # number of path-segment records/candidates
for seg in d["segments"]:
    print(seg["index"], seg["barrier_kcal"], seg["delta_kcal"])
```

The lightweight `segments` records carry MEP barrier/delta/bond-change data.
Requested post-processing is recorded separately in `post_segments`; match
the two lists by `index` because bridge/skipped segments can make their list
positions differ. See the schema skill before assuming nested keys.

## Re-running individual stages

To rerun a specific stage (for example after a walltime hit), call the
standalone subcommands directly on the segment outputs `all` produced.
Each standalone subcommand resolves charge on its own, so an `.xyz` input
needs an explicit `-q <total_charge>`; take the value from the parent run's
`summary.json` (`d["charge"]`). Spin defaults to 1 (`-m` to override):

```bash
cd result_all
RUN_ID="$(date +%Y%m%d-%H%M%S)"
TOTAL_CHARGE=-1  # replace with summary.json["charge"]
pdb2reaction tsopt -i _work/path_opt/hei_seg_03.xyz -q "$TOTAL_CHARGE" -o "segments/seg_03/retry_${RUN_ID}/ts" -b uma
pdb2reaction irc   -i "segments/seg_03/retry_${RUN_ID}/ts/final_geometry.xyz" -q "$TOTAL_CHARGE" -o "segments/seg_03/retry_${RUN_ID}/irc" -b uma
pdb2reaction freq  -i "segments/seg_03/retry_${RUN_ID}/ts/final_geometry.xyz" -q "$TOTAL_CHARGE" -o "segments/seg_03/retry_${RUN_ID}/freq" -b uma
```

Inspect every retry `result.json` and the structures before deliberately
adopting it as a canonical segment result. A retry is kept in a unique
directory so it cannot silently overwrite the original partial run.

## Caveats

- `--scan-lists` is a Python literal-eval expression. Most
  shell-quoting trouble traces back to single vs double quotes. Put all
  stage literals after one flag occurrence; do not repeat `-s`.
- If `summary.json` shows `"status": "failed"` (or `"partial"`) for any segment, look
  at the corresponding `summary.log` block; per-stage errors are also
  duplicated into `segments/seg_NN/<stage>/result.json`.
- `segments/seg_NN/` is created when post-processing starts and may therefore
  be partial after a failed TSOPT/IRC/freq/DFT stage. Presence of the directory
  is not a success signal; check `summary.json` and the stage `result.json`.
  `all` starts IRC only after the TS stage reports `status=converged` and
  `n_imaginary_modes=1`; a failed or higher-order TS attempt stops the pipeline
  before IRC and endpoint optimization.
  MEP scratch remains under `<work_path>/seg_NNN_<tag>/` (3-digit;
  `_work/path_opt/` by default, `_work/path_search/` with
  `--refine-path`).
- All four built-in backends implement analytical Hessians. The special
  restriction is UMA's parallel predictor: an explicit analytical request
  with `--workers > 1` is an error, not an automatic finite-difference fallback.

## See also

- `all-endpoint-mep.md`, `all-scan-list.md`, `all-ts-only.md` — three
  invocation modes.
- `extract.md`, `path-search.md`, `tsopt.md`, `irc.md`, `freq.md`,
  `dft.md` — the underlying subcommands.
- `pdb2reaction-workflows-output/SKILL.md` — output schema and
  R/TS/P coordinate conventions.
- Defaults: `import pdb2reaction.core.defaults` (`OUT_DIR_ALL`, plus the
  per-stage `*_KW` dicts).
