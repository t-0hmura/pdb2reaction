---
name: pdb2reaction-cli
description: Per-subcommand reference for pdb2reaction's 18 CLI subcommands (extract / path-search / tsopt / freq / irc / dft / scan / opt / sp / all / …). SKILL.md is a 1-line input→output cheatsheet; most subcommands also have their own md (`extract.md` / `tsopt.md` / …) for flags, validation, caveats. See `freeze-atoms.md` for cluster-boundary frozen-atom mechanics. TRIGGER on questions about a specific subcommand, flag, or shell invocation. SKIP for install / HPC / output-parsing / structure-format-editing questions.
---

# pdb2reaction CLI

## Cheatsheet (input → output)

Charge and multiplicity must be chemically resolved for every run even where
the CLI has neutral/singlet defaults. Use `-l 'RES:Q,...'` for PDB/mmCIF residue-based
charge derivation or a verified total `-q <int>`, and set `-m <int>` explicitly
for open-shell systems. The generic examples below use neutral singlet
`-q 0 -m 1`; replace those values rather than copying them blindly. `-b`
defaults to `uma`. Full flags: the matching `<sub>.md` next to this file, or
`pdb2reaction <sub> --help`.

| sub | role | minimal command | primary output |
|---|---|---|---|
| `all` | Orchestrate the selected optional extraction, path/scan, TS/IRC, thermo, and DFT stages | `pdb2reaction all -i 1.R.pdb 3.P.pdb -q 0 -m 1 --tsopt --thermo -o out` | successful segment deliverables + `out/summary.json` once the summary writer is reached |
| `all` (scan-list) | Single reactant + staged scans | `pdb2reaction all -i 1.R.pdb -q 0 -m 1 -s '[(a,b,1.6)]' --tsopt -o out` | as above |
| `all` (ts-only) | Pre-existing TS candidate | `pdb2reaction all -i ts.xyz -q 0 -m 1 --tsopt --thermo -o out` | `out/segments/seg_01/{ts,irc,freq}/...` + `out/segments/seg_01/structures/*.xyz` (`.pdb`/`.cif` only with topology) |
| `extract` | Active-site cluster cut | `pdb2reaction extract -i raw.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' -r 2.6 -o cluster.pdb` | `cluster.pdb` (`-o` is the output file path, not a directory) |
| `path-search` | Recursive MEP w/ bond-change segmentation | `pdb2reaction path-search -i 1.R.pdb 3.P.pdb -q 0 -m 1 -o out` | `out/hei_seg_NN.xyz` + `out/summary.json` |
| `path-opt` | Single-segment MEP refinement | `pdb2reaction path-opt -i 1.R.pdb 2.P.pdb -q 0 -m 1 -o out` | `out/final_geometries_trj.xyz` |
| `opt` | Geometry minimization (L-BFGS / RFO) | `pdb2reaction opt -i geom.pdb -q 0 -m 1 -o out` | `out/final_geometry.xyz` |
| `tsopt` | TS optimization (RS-P-RFO / Dimer) | `pdb2reaction tsopt -i ts.xyz -q 0 -m 1 -o out` | `out/final_geometry.xyz`; one imaginary mode is necessary but mode displacement and IRC connectivity are also required for TS validation |
| `freq` | Hessian + QRRHO thermo | `pdb2reaction freq -i geom.xyz -q 0 -m 1 -o out` | after successful evaluation, `out/frequencies_cm-1.txt`; `out/thermoanalysis.yaml` with `--dump` (`all --thermo` sets it) |
| `sp` | Single-point MLIP energy + forces (+optional Hessian) | `pdb2reaction sp -i geom.pdb -q 0 -m 1 -o out` | `out/forces.npy` (+ `out/hessian.npy` with `--hess`); energy printed to stdout; `out/result.json` + `out/summary.json` only with `--out-json` |
| `irc` | IRC from a TS | `pdb2reaction irc -i ts.xyz -q 0 -m 1 -o out` | `out/{forward,backward,finished}_irc_trj.xyz` |
| `dft` | Single-point DFT (PySCF / GPU4PySCF) | `pdb2reaction dft -i geom.pdb -q 0 -m 1 --func-basis 'wb97m-v/def2-tzvpd' -o out` | after a successful calculation, `out/result.yaml`; `out/result.json` with `--out-json` |
| `scan` | 1D distance scan w/ restraints | `pdb2reaction scan -i 1.R.pdb -q 0 -m 1 -s '[(a,b,1.6)]' -o out` | `out/scan_trj.xyz`, per-stage `stage_NN/result.xyz` |
| `scan2d` | 2D distance grid scan | `pdb2reaction scan2d -i 1.R.pdb -q 0 -m 1 -s '[(a,b,1.3,3.1),(c,d,1.2,3.2)]' -o out` | grid records + `out/surface.csv`; `out/scan2d_map.png` when interpolation/export succeeds |
| `scan3d` | 3D distance grid scan | `pdb2reaction scan3d -i 1.R.pdb -q 0 -m 1 -s '[(a,b,L,H),(c,d,L,H),(e,f,L,H)]' -o out` | grid records + `out/surface.csv`; `out/scan3d_density.html` when export succeeds |
| `trj2fig` | Energy profile from XYZ trj | `pdb2reaction trj2fig -i trj.xyz` | `energy.png` (default when no `-o`) |
| `energy-diagram` | Diagram from energy values | `pdb2reaction energy-diagram -i "[0.0, 21.5, -0.7]" --label-x "['R','TS','P']"` | `energy_diagram.png` |
| `add-elem-info` | Add PDB element column (cols 77-78) | `pdb2reaction add-elem-info -i raw.pdb -o fixed.pdb` | `fixed.pdb` |
| `fix-altloc` | Resolve PDB alternate locations | `pdb2reaction fix-altloc -i raw.pdb -o fixed.pdb` | `fixed.pdb` (single conformation per residue) |
| `bond-summary` | Diff bonds between consecutive structures | `pdb2reaction bond-summary -i reactant.pdb -i product.pdb` (or positional `R.pdb P.pdb`) | stdout text by default; JSON to stdout with `--json` |

For mmCIF or oversized-PDB input, geometry workflows keep a normalized `.pdb`
for communication between pipeline stages and, with conversion enabled, add
`.cif` companions carrying the original chain and residue identifiers. The
`.pdb` names in the table are
therefore workflow paths, not a claim that CIF metadata is discarded.

## Cross-cutting topic guides

| md | Topic |
|---|---|
| `freeze-atoms.md` | Cluster-boundary frozen atoms — cap hydrogens (`LKH/HL`), `--freeze-links`, `--freeze-atoms`, YAML `geom.freeze_atoms`. The three sources are unioned; use a chemically justified boundary and inspect it rather than assuming every cluster needs the same freeze set. |

## Common flag conventions

| Flag | Meaning |
|---|---|
| `-i, --input` | Input file(s); geometry workflows accept `.pdb`, `.cif` / `.mmcif`, `.xyz`, `.gjf` |
| `-q, --charge` | Total charge (integer) |
| `-l, --ligand-charge` | `'RES1:Q1,RES2:Q2'` per-residue mapping (PDB/mmCIF metadata) |
| `-m, --multiplicity` | Spin multiplicity (2S+1), default 1 |
| `-b, --backend` | MLIP backend: `uma` / `orb` / `mace` / `aimnet2` |
| `-o, --out-dir` | Output directory, subcommand-specific default |
| `--config` | YAML configuration file applied before CLI flags |
| `--show-config` | Print the resolved configuration, then **continue** with the full run |
| `--dry-run` | Validate options and print the execution plan, then exit before MLIP/DFT stages. Special case: `all -c/--center ... --dry-run` runs extraction in a temporary directory so it can validate the derived charge and electron parity; it does not run scan/MEP/TSOPT/IRC/freq/DFT. |
| `--help-advanced` | Reveal hidden / advanced flags |
| `--ref-pdb` | Reference PDB/mmCIF used to derive residue context for XYZ/GJF inputs while retaining their coordinates |
| `--solvent` | Experimental, computationally expensive xTB correction: `E_xTB(solvent) - E_xTB(vacuum)`. `none` disables it; see `pdb2reaction-install-backends/xtb.md`. |

For calculation commands, explicit `-q` takes precedence over
`-l 'RES:Q'` derivation, then config/defaults. This includes
`all -c/--center`: extraction derives the cluster charge, but explicit `-q`
sets the total and a mismatch is reported as a warning. A per-resname `-l` mapping can
derive charge from a residue-bearing PDB/mmCIF whether or not extraction runs; bare
XYZ needs an explicit total, while a valid GJF supplies its header value.

## Canonical recipes

### Multi-input MEP for a 1-step reaction

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    --out-dir result_mep
```

### Single-input scan-list (when only the reactant is available)

```bash
pdb2reaction all -i 1.R.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --scan-lists '[("CS1 SAM 320","C7 GPP 321",1.60)]' \
                 '[("H11 GPP 321","OE2 GLU 186",0.90)]' \
    --tsopt --thermo \
    --out-dir result_scan
```

### Validate a TS candidate without re-running path search

```bash
pdb2reaction tsopt -i ts_guess.xyz -q -1 -m 1 -b uma -o result_tsopt
pdb2reaction freq  -i result_tsopt/final_geometry.xyz -q -1 -m 1 -b uma -o result_freq
pdb2reaction irc   -i result_tsopt/final_geometry.xyz -q -1 -m 1 -b uma -o result_irc
```

### DFT//MLIP single point on the highest-local-barrier TS candidate

```bash
SUMMARY=result_mep/summary.json
RLS_SEG=$(python - "$SUMMARY" <<'PY'
import json, sys
with open(sys.argv[1], encoding="utf-8") as handle:
    print(int(json.load(handle)["rate_limiting_step"]["segment"]))
PY
)
printf -v RLS_DIR 'seg_%02d' "$RLS_SEG"
TS_FILE="result_mep/segments/${RLS_DIR}/ts.pdb"
test -f "$TS_FILE"
pdb2reaction dft -i "$TS_FILE" \
    -l 'SAM:1,GPP:-3' \
    --func-basis 'wb97m-v/def2-tzvpd' \
    --engine gpu
```

### Bond-change report between R and P

```bash
pdb2reaction bond-summary -i reactant.pdb -i product.pdb
```

## Cross-cutting caveats

| Pitfall | Fix |
|---|---|
| `--scan-lists` syntax error | The list is a Python literal-eval expression. Quote with single-quotes outside, double-quotes inside, and do not confuse the backtick (`` ` ``) with the backslash (`\`). |
| Wrong charge silently | Resolve protonation/oxidation state and inspect the residue charge breakdown. `--dry-run` reports the resolved value and checks electron parity but cannot prove that the chemistry is correct. For `all -c/--center`, it performs temporary extraction, validates the derived charge/parity, removes the temporary data, and exits before computational stages. |
| Forgetting to pin `-b` for production | The default is `-b uma`; specify `-b uma` / `-b orb` / `-b mace` / `-b aimnet2` explicitly so a future default change cannot silently re-route the run. |
| `--config` YAML ignored | YAML is read **after** built-in defaults but **before** explicit CLI flags. Anything also given on CLI overrides YAML. |
| `--help-advanced` flags differ between versions | They are subject to change; if a flag isn't in `--help`, check `--help-advanced` and version-pin if the workflow is shared. |
| OOM on the Hessian step | Try `--hessian-calc-mode FiniteDifference` to avoid the analytical autograd graph, use a justified frozen boundary/PHVA, or select a smaller backend model. The active-space Hessian itself remains dense, so benchmark memory rather than assuming any one switch is sufficient. |
| UMA `--workers > 1` with an explicit analytical Hessian | This raises `BackendError`; it never silently changes the requested method. Use `--workers 1` for `Analytical`, or explicitly select `FiniteDifference`. ORB/MACE/AIMNet2 do not use these worker flags. All four built-in backends implement analytical Hessians when used in a supported configuration. |

## Defaults

Most calculation defaults are exported from `pdb2reaction.core.defaults`; Click-only presentation defaults may live in the subcommand module. Use live `--help` plus the relevant config dict rather than assuming every flag is in one dict:

```bash
python -c "import pdb2reaction.core.defaults as d; print([n for n in dir(d) if n.endswith('_KW') or n.startswith('OUT_DIR')])"

# Examples:
python -c "import pdb2reaction.core.defaults as d; print(d.LBFGS_KW)"
python -c "import pdb2reaction.core.defaults as d; print(d.RSIRFO_KW)"
python -c "import pdb2reaction.core.defaults as d; print(d.IRC_KW)"
python -c "import pdb2reaction.core.defaults as d; print(d.UMA_CALC_KW)"
```

Each per-subcommand md points at the relevant `_KW` dict in the
"See also" section.

## See also

- `pdb2reaction-overview/SKILL.md` — what `pdb2reaction` is and when to
  use it.
- `pdb2reaction-structure-io/` — input file formats and charge / spin.
- `pdb2reaction-install-backends/` — `<tool>` / backend installation.
- `pdb2reaction-workflows-output/SKILL.md` — what comes out of each
  invocation, six canonical workflows, energy diagrams. The `summary.json`
  schema, R/TS/P canonical paths, and bond-change interpretation are in
  [`pdb2reaction-workflows-output/summary-json.md`](../pdb2reaction-workflows-output/summary-json.md).
- `pdb2reaction-hpc/SKILL.md` — running these recipes on PBS / SLURM.
