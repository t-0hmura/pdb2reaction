# `pdb2reaction all` — base orientation

## Purpose

`all` is the meta-command that chains the entire workflow:
extract → path-opt → tsopt → irc → freq → (optional) dft. The MEP
stage runs single-pass `path-opt` by default; pass `--refine-path True`
to run recursive `path-search` instead. `all`
resolves three input modes via flag context (see the three companion
mds: `all-endpoint-mep.md`, `all-scan-list.md`, `all-ts-only.md`).

Use `all` when you want a **single qsub-able invocation** that
produces R / TS / P / IM coordinates plus barrier numbers for one or
more elementary steps.

## Synopsis

```bash
pdb2reaction all -i <input(s)> [-c <substrate>] [-l 'RES:Q,...'] \
    [--scan-lists '...'] [--tsopt True] [--thermo True] [--dft True] \
    [-b uma|orb|mace|aimnet2] [-o result_all/]
```

## Key flags (cross-mode)

| Flag | Type | Default | Description |
|---|---|---|---|
| `-i, --input` | path(s) | required | One or more reaction-ordered structures, or a TS-candidate alone |
| `-c, --center` | str | (uses input as-is) | Substrate selector: residue-name list `'RES1,RES2,...'`, residue-ID list `'A:44,B:321'`, or a PDB path. Chain-qualified residue *names* (`'B:SAM'`) are not supported — use the residue ID instead. |
| `-l, --ligand-charge` | str | none | Per-resname charges, e.g. `'SAM:1,GPP:-3'` (or a bare number = total). The per-resname mapping is honored **whether or not extraction runs**: with `-c` it feeds the extractor summary; with `-c` omitted (extraction skipped, pre-carved model passed as-is) the same mapping is applied to the **full input PDB** to derive the total. Prefer `-l` — `-q` is rarely needed. |
| `-q, --charge` | int | derived from `-l` | Total system charge. With `-c` (extraction runs) `-q` acts as an **assertion**: it must **match** the extract-derived charge, otherwise the run aborts with `BadParameter` (this is also checked by `--dry-run`). With `-c` omitted (pre-carved input passed as-is) `-q` sets the total directly and **emits a warning**. Normally unnecessary — let `-l` derive the total. |
| `-m, --multiplicity` | int | 1 | Spin multiplicity (2S+1) |
| `-r, --radius` | float | 2.6 | Pocket radius (Å) when `-c` triggers extraction |
| `-s, --scan-lists` | repeated | none | Staged distance scans (mode 2 — `all-scan-list.md`) |
| `--refine-path` | BOOL | `False` | MEP engine: `False` runs single-pass `path-opt`; `True` runs recursive `path-search`, often improving a poor TS seed but potentially splitting a bad path into unnecessary segments and greatly increasing cost |
| `--tsopt` | BOOL | `False` | Run TS optimization + IRC per reactive segment (also required to enter TS-only mode with a single `-i`) |
| `--flatten/--no-flatten` | flag | off | Enable surplus-imaginary-mode cleanup when TSOPT does not reach a first-order saddle |
| `--irc-step-size` | float | IRC default `0.10` | Forward a smaller EulerPC maximum step; try `0.05` when an IRC branch stops after only a few frames |
| `--irc-never-stop/--no-irc-never-stop` | flag | off | Ignore IRC energy-rise/plateau stops only; gradient/integrator convergence and max cycles remain active |
| `--thermo` | BOOL | `False` | Run freq + thermochemistry on R / TS / P |
| `--dft` | BOOL | `False` | Run DFT single point on R / TS / P |
| `--dft-func-basis` | str | `wb97m-v/def2-tzvpd` | DFT functional/basis (when `--dft`) |
| `-b, --backend` | str | `uma` | MLIP backend |
| `--solvent` | str | none | xTB-ALPB solvent name (`water`, `methanol`, …) |
| `-o, --out-dir` | path | `./result_all/` | Top-level output directory |
| `--config` | path | none | YAML config applied before CLI flags |
| `--show-config` | flag | off | Print the resolved config and continue running |
| `--dry-run` | flag | off | Validate inputs + print the execution plan, then exit without running |
| `--help-advanced` | flag | — | Reveal hidden flags (freeze, advanced overrides) |

Run `pdb2reaction all --help-advanced` for the full list (it changes
between versions).

## Mode selection cheatsheet

| Inputs | Mode md | Behavior |
|---|---|---|
| Single `-i input.{xyz,pdb,gjf}` + `--tsopt` (no `--scan-lists`) | `all-ts-only.md` | treat input as TS candidate; tsopt + irc + freq |
| Single `-i input.pdb` + `--scan-lists '...'` | `all-scan-list.md` | single reactant + staged distance scans |
| Multiple `-i 1.R.pdb [2.IM.pdb …] N.P.pdb` (reaction-ordered) | `all-endpoint-mep.md` | multi-endpoint MEP |

A single `-i` without **either** `--scan-lists` or `--tsopt` raises
`BadParameter` (see `all.py`: "Provide at least two structures... or a
single structure with --scan-lists, or a single structure with --tsopt
True").

## Output tree (typical)

Three zones: deliverables at `<out_dir>/`, per-segment deliverables under `<out_dir>/segments/seg_NN/`, scratch under `<out_dir>/_work/`. `<work_path>` = `_work/path_opt/` (default) or `_work/path_search/` (`--refine-path True`).

| Path | When | Content |
|---|---|---|
| `<out_dir>/summary.json` | always | machine-readable per-stage results |
| `<out_dir>/summary.log` | always | human-readable text + dir tree |
| `<out_dir>/mep.pdb`, `mep_trj.xyz` | `mep_trj.xyz` always; `mep.pdb` for PDB input | stitched MEP across segments (promoted to root) |
| `<out_dir>/mep_w_ref.pdb` | `-c` extraction with all-PDB inputs (auto) | MEP merged into the full-system template |
| `<out_dir>/energy_diagram_MEP.png` | always | bare all-segment MEP energies |
| `<out_dir>/energy_diagram_{UMA,G_UMA,DFT,G_DFT_plus_UMA}_all.png` | combos of `--thermo`/`--dft` | aggregated multi-segment diagrams |
| `<out_dir>/segments/seg_NN/{reactant,ts,product}.{pdb,xyz}` | always | canonical R/TS/P (2-digit) |
| `<out_dir>/segments/seg_NN/ts/final_geometry.{xyz,pdb}`, `vib/imag_*.pdb` | `--tsopt` | tsopt output |
| `<out_dir>/segments/seg_NN/irc/{forward,backward,finished}_irc_trj.xyz` | `--tsopt` | IRC trajectories |
| `<out_dir>/segments/seg_NN/freq/{R,TS,P}/{frequencies_cm-1.txt, thermoanalysis.yaml}` | `--thermo` | per-state freq + thermo |
| `<out_dir>/segments/seg_NN/dft/{R,TS,P}/result.{yaml,json}` | `--dft` | per-state DFT |
| `<out_dir>/segments/seg_NN/energy_diagram_{UMA,G_UMA,DFT,G_DFT_plus_UMA}.png` | combos of `--thermo`/`--dft` | per-segment diagrams |
| `<out_dir>/_work/models/model_<stem>.pdb` | `-c` given | extracted active-site clusters (one per `-i` input) |
| `<work_path>/seg_NNN_<tag>/` | always | per-string MEP scratch (3-digit, `_mep` / `_maxdepth` / `_bridge`) |
| `<work_path>/mep_seg_NN_trj.xyz`, `mep_seg_NN.{pdb,gjf}` | always | per-segment MEP frames |
| `<work_path>/hei_seg_NN.{xyz,pdb,gjf}` | always | HEI candidate per segment (TS seed) |

## Output keys (summary.json — top level)

```python
import json
d = json.load(open("result_all/summary.json"))
print(d["status"])                    # "success" / "partial" / "failed"
print(d["pdb2reaction_version"])
print(d["charge"], d["spin"])
print(d["rate_limiting_step"])        # which segment is rate-limiting
print(len(d["segments"]))             # number of elementary steps
for seg in d["segments"]:
    print(seg["index"], seg["barrier_kcal"], seg["delta_kcal"])
```

Per-segment fields include `barrier_kcal`, `delta_kcal`, `bond_changes`,
`structures` (paths to reactant/ts/product files), `tsopt` /
`irc` / `freq` / `dft` sub-objects with their own `status` and energies.

## Re-running individual stages

To rerun a specific stage (for example after a walltime hit), call the
standalone subcommands directly on the segment outputs `all` produced.
Each standalone subcommand resolves charge on its own, so an `.xyz` input
needs an explicit `-q <total_charge>`; take the value from the parent run's
`summary.json` (`d["charge"]`). Spin defaults to 1 (`-m` to override):

```bash
pdb2reaction tsopt -i _work/path_opt/hei_seg_03.xyz -q <total_charge> -o segments/seg_03/ts -b uma
pdb2reaction irc   -i segments/seg_03/ts/final_geometry.xyz -q <total_charge> -o segments/seg_03/irc -b uma
pdb2reaction freq  -i segments/seg_03/ts/final_geometry.xyz -q <total_charge> -o segments/seg_03/freq -b uma
```

The directory layout matches what `all` produces, so downstream
analysis scripts keep working.

## Caveats

- `--scan-lists` is a Python literal-eval expression. Most
  shell-quoting trouble traces back to single vs double quotes.
- If `summary.json` shows `"status": "failed"` (or `"partial"`) for any segment, look
  at the corresponding `summary.log` block; per-stage errors are also
  duplicated into `segments/seg_NN/<stage>/result.json`.
- The `segments/seg_NN/` deliverable directory is **only populated on success**
  for that segment. Failed segments leave scratch artifacts under
  `<work_path>/seg_NNN_<tag>/` (3-digit; `_work/path_opt/` by default,
  `_work/path_search/` with `--refine-path True`) but not the `segments/seg_NN/`
  2-digit deliverable copy.

## See also

- `all-endpoint-mep.md`, `all-scan-list.md`, `all-ts-only.md` — three
  invocation modes.
- `extract.md`, `path-search.md`, `tsopt.md`, `irc.md`, `freq.md`,
  `dft.md` — the underlying subcommands.
- `pdb2reaction-workflows-output/SKILL.md` — output schema and
  R/TS/P coordinate conventions.
- Defaults: `import pdb2reaction.core.defaults` (`OUT_DIR_ALL`, plus the
  per-stage `*_KW` dicts).
