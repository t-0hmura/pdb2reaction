# `pdb2reaction all` — base orientation

## Purpose

`all` is the meta-command that chains the entire workflow:
extract → path-search → tsopt → irc → freq → (optional) dft. It
resolves three input modes via flag context (see the three companion
mds: `all-endpoint-mep.md`, `all-scan-list.md`, `all-ts-only.md`).

Use `all` when you want a **single qsub-able invocation** that
produces R / TS / P / IM coordinates plus barrier numbers for one or
more elementary steps.

## Synopsis

```bash
pdb2reaction all -i <input(s)> [-c <substrate>] [-l 'RES:Q,...'] \
    [--scan-lists '...'] [--tsopt] [--thermo] [--dft] \
    [-b uma|orb|mace|aimnet2] [-o result_all/]
```

## Key flags (cross-mode)

| Flag | Type | Default | Description |
|---|---|---|---|
| `-i, --input` | path(s) | required | One or more reaction-ordered structures, or a TS-candidate alone |
| `-c, --center` | str | (uses input as-is) | Substrate selector: residue-name list `'RES1,RES2,...'`, residue-ID list `'A:44,B:321'`, or a PDB path. Chain-qualified residue *names* (`'B:SAM'`) are not supported — use the residue ID instead. |
| `-l, --ligand-charge` | str | none | Per-residue charges, e.g. `'SAM:1,GPP:-3'` |
| `-q, --charge` | int | derived from `-l` | Total cluster charge override |
| `-m, --multiplicity` | int | 1 | Spin multiplicity (2S+1) |
| `-r, --radius` | float | 2.6 | Pocket radius (Å) when `-c` triggers extraction |
| `-s, --scan-lists` | repeated | none | Staged distance scans (mode 2 — `all-scan-list.md`) |
| `--tsopt / --no-tsopt` | toggle | `--no-tsopt` | Run TS optimization + IRC per reactive segment (also required to enter TS-only mode with a single `-i`) |
| `--thermo / --no-thermo` | toggle | `--no-thermo` | Run freq + thermochemistry on R / TS / P |
| `--dft / --no-dft` | toggle | `--no-dft` | Run DFT single point on R / TS / P |
| `--dft-func-basis` | str | `wb97m-v/def2-tzvpd` | DFT functional/basis (when `--dft`) |
| `-b, --backend` | str | `uma` | MLIP backend |
| `--solvent` | str | none | xTB-ALPB solvent name (`water`, `methanol`, …) |
| `-o, --out-dir` | path | `./result_all/` | Top-level output directory |
| `--config` | path | none | YAML config applied before CLI flags |
| `--show-config` / `--dry-run` | flag | off | Print resolved config and exit |
| `--help-advanced` | flag | — | Reveal hidden flags (resume / checkpoint / freeze options) |

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

`<path_dir>` = `path_search/` (default) or `path_opt/` (`--refine-path False`).

| Path | When | Content |
|---|---|---|
| `<out_dir>/summary.json` | always | machine-readable per-stage results |
| `<out_dir>/summary.log` | always | human-readable text + dir tree |
| `<out_dir>/models/model_<stem>.pdb` | `-c` given | extracted active-site clusters (one per `-i` input) |
| `<out_dir>/seg_NN/{reactant,ts,product}.{pdb,xyz}` | always | canonical R/TS/P (2-digit, top-level) |
| `<out_dir>/mep_trj.xyz`, `mep.{pdb,gjf}` | always | stitched MEP across segments |
| `<out_dir>/energy_diagram_{UMA,G_UMA,DFT,G_DFT_plus_UMA}_all.png` | combos of `--thermo`/`--dft` | aggregated multi-segment diagrams |
| `<path_dir>/seg_NNN_<tag>/` | always | per-string MEP scratch (3-digit, `_mep` / `_maxdepth` / `_bridge`) |
| `<path_dir>/mep_seg_NN_trj.xyz`, `mep_seg_NN.{pdb,gjf}` | always | canonical per-segment MEP frames |
| `<path_dir>/hei_seg_NN.{xyz,pdb,gjf}` | always | HEI candidate per segment (TS seed) |
| `<path_dir>/energy_diagram_MEP.png` | always | bare MEP energies (path-search level) |
| `<path_dir>/post_seg_NN/ts/final_geometry.{xyz,pdb}`, `vib/imag_*.pdb` | always | tsopt output |
| `<path_dir>/post_seg_NN/irc/{forward,backward,finished}_irc_trj.xyz` | always | IRC trajectories |
| `<path_dir>/post_seg_NN/freq/{R,TS,P}/{frequencies_cm-1.txt, thermoanalysis.yaml}` | always | per-state freq + thermo |
| `<path_dir>/post_seg_NN/dft/{R,TS,P}/result.{yaml,json}` | `--dft` | per-state DFT |
| `<path_dir>/post_seg_NN/energy_diagram_{UMA,G_UMA,DFT,G_DFT_plus_UMA}.png` | combos of `--thermo`/`--dft` | per-segment diagrams |

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

## Resume / restart

`pdb2reaction all` supports `--resume` to skip stages whose outputs already exist. For finer-grained re-runs, call subcommands directly:

```bash
pdb2reaction tsopt -i path_search/hei_seg_03.xyz -o path_search/post_seg_03/ts -b uma
pdb2reaction irc   -i path_search/post_seg_03/ts/final_geometry.xyz -o path_search/post_seg_03/irc -b uma
pdb2reaction freq  -i path_search/post_seg_03/ts/final_geometry.xyz -o path_search/post_seg_03/freq -b uma
```

The directory layout matches what `all` produces, so downstream
analysis scripts keep working.

## Caveats

- `--scan-lists` is a Python literal-eval expression. Most
  shell-quoting trouble traces back to single vs double quotes.
- If `summary.json` shows `"status": "failed"` (or `"partial"`) for any segment, look
  at the corresponding `summary.log` block; per-stage errors are also
  duplicated into `post_seg_NN/<stage>/result.json`.
- The `seg_NN/` top-level directory is **only populated on success**
  for that segment. Failed segments leave artifacts under
  `path_search/seg_NNN_<tag>/` and `path_search/post_seg_NN/`
  (3-digit scratch) but not the top-level 2-digit copy.

## See also

- `all-endpoint-mep.md`, `all-scan-list.md`, `all-ts-only.md` — three
  invocation modes.
- `extract.md`, `path-search.md`, `tsopt.md`, `irc.md`, `freq.md`,
  `dft.md` — the underlying subcommands.
- `pdb2reaction-workflows-output/SKILL.md` — output schema and
  R/TS/P coordinate conventions.
- Defaults: `import pdb2reaction.defaults` (`OUT_DIR_ALL`, plus the
  per-stage `*_KW` dicts).