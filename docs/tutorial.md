# Tutorials

Use this page when you want example-first guidance.

## Recommended order

1. {ref}`bezA Enzyme Case Study <tutorial-case-study>`
2. {ref}`Smoke Test Matrix <tutorial-smoke-tests>`

```{important}
Treat TS outputs as initial **candidates**. For enzyme reactions, iterative refinement is common (endpoint quality, pocket definition, constraints, scan targets), and TS validation with both `freq` and `irc` is required before interpretation.
```

(tutorial-case-study)=
## bezA Enzyme Case Study

This tutorial explains how to read a full `pdb2reaction all` output using a real bezA enzyme example.

Packaging note: bundled `run.out` / `summary.yaml` paths are normalized to `./bezA_case_study` so the example is portable across environments.

### What You Will Learn

- how a single-structure + staged-scan workflow is configured,
- how to read `summary.log` / `summary.yaml` quickly,
- how to separate TS **candidates** from validated TS structures,
- what to tune next when multiple imaginary frequencies remain.

### Bundled Assets

- `_static/bezA_case_study/r_input.pdb`
- `_static/bezA_case_study/summary.log`
- `_static/bezA_case_study/summary.yaml`
- `_static/bezA_case_study/run.out`
- `_static/bezA_case_study/ts_seg_01.pdb`
- `_static/bezA_case_study/ts_seg_03.pdb`
- `_static/bezA_case_study/energy_diagram_MEP.png`
- `_static/bezA_case_study/energy_diagram_UMA_all.png`
- `_static/bezA_case_study/energy_diagram_G_UMA_all.png`
- `_static/bezA_case_study/irc_plot_all.png`

### Archived Command (Reference)

```bash
pdb2reaction -i r.pdb -c 'SAM,GPP,MG' --ligand-charge 'SAM:1,GPP:-3' \
  --scan-lists '[ ("CS1 SAM 320","GPP 321 C7",1.60) ]' \
               '[ ("GPP`321/H11","GLU`186/OE2",0.90) ]' \
  --tsopt True --thermo True --out-dir bezA_case_study
```

In the bundled tutorial assets, this initial structure is provided as `r_input.pdb`.

Pipeline mode: `path-search` with recursive refinement.

### Reproduce From This Repository

From the repository root:

```bash
pdb2reaction -i docs/_static/bezA_case_study/r_input.pdb \
  -c 'SAM,GPP,MG' --ligand-charge 'SAM:1,GPP:-3' \
  --scan-lists '[ ("CS1 SAM 320","GPP 321 C7",1.60) ]' \
               '[ ("GPP`321/H11","GLU`186/OE2",0.90) ]' \
  --tsopt True --thermo True --out-dir ./bezA_case_study
```

Note: this command creates `./bezA_case_study/` in your current working directory. It is separate from the bundled assets under `docs/_static/bezA_case_study/`.

This reproduces the same workflow type as the archived case study using the bundled initial structure.

### Key Results

From `summary.log`:

```text
Number of MEP images : 29
Number of segments   : 3
```

MEP-state summary (kcal/mol):
- `TS1`: `+14.90`
- `IM1_1`: `-36.65`
- `IM1_2`: `-38.68`
- `TS2`: `+0.21`
- `P`: `-47.13`

Per-segment TS checks:
- Segment 01: `n_imag = 2` (`-205.1 cm^-1` max imaginary)
- Segment 03: `n_imag = 3` (`-1067.2 cm^-1` max imaginary)

Interpretation: the workflow completed and produced coherent candidates, but TS refinement is still required before final mechanistic claims.

### Figures

#### Global MEP

![MEP energy diagram](_static/bezA_case_study/energy_diagram_MEP.png)

#### UMA energies (all segments)

![UMA energy diagram](_static/bezA_case_study/energy_diagram_UMA_all.png)

#### UMA Gibbs energies (all segments)

![UMA Gibbs energy diagram](_static/bezA_case_study/energy_diagram_G_UMA_all.png)

#### IRC overview

![IRC plot](_static/bezA_case_study/irc_plot_all.png)

### Reading Guide

#### What is already strong

- End-to-end outputs were generated: MEP, TSOPT, IRC, and thermochemistry.
- Bond-change detection captured chemically meaningful events in segment 01 and 03.
- Energy summaries are internally consistent across `summary.log` and `summary.yaml`.

#### What still needs iteration

- `n_imag > 1` means current TS structures are candidates, not validated TSs.
- Additional TS-focused optimization and validation are required.

### Required TS Validation Checklist

1. Confirm a single imaginary frequency per TS candidate (`freq`).
2. Confirm the imaginary mode follows the intended reaction coordinate.
3. Run IRC in both directions and verify endpoint minima.
4. Confirm endpoint connectivity matches intended R/P chemistry.

Related commands:
- [tsopt](tsopt.md)
- [freq](freq.md)
- [irc](irc.md)
- [path-search](path_search.md)

### Typical Next Adjustments

1. Retry TS refinement with `--opt-mode heavy`.
2. Enable `--flatten-imag-mode True`.
3. Increase MEP resolution with `--max-nodes` for difficult segments.
4. Revisit scan targets and initial structures for improved TS guesses.

### Next

- Continue to {ref}`Smoke Test Matrix <tutorial-smoke-tests>`.
- For setup and CLI basics, see [Getting Started](getting-started.md).

(tutorial-smoke-tests)=
## Smoke Test Matrix

This section maps `test/run.sh` cases to practical goals.

This guide explains the script by purpose so readers can choose checks without relying on internal run labels.

### Locations

- Script: `test/run.sh`
- Developer note: `test/test.md`
- Inputs: `test/*.pdb`, `test/*.xyz`, `test/*.gjf`

Run from repository root:

```bash
bash test/run.sh
```

### Purpose-Based Matrix

#### Format and baseline pipeline checks

- PDB full pipeline + TS/thermo check
  - `pdb2reaction -i r.pdb p.pdb -q -1 --tsopt True --thermo True`
- XYZ format check
  - `pdb2reaction -i r.xyz p.xyz -q -1`
- GJF metadata / template behavior
  - `pdb2reaction -i r.gjf p.gjf ...`

#### Single-input staged scan entry

- PDB scan entry
- GJF scan entry
- XYZ scan entry

These confirm selector parsing and staged-scan behavior.

#### TS-focused flows

- TS+thermo from PDB
- TS from PDB in `light` mode
- TS from GJF template

Use these when you want TS behavior checks without full multistructure MEP setup.

#### Standalone subcommand checks

- `dft` CPU quick check
- `freq` quick check
- `tsopt` quick check

Use these for targeted debugging when `all` is too broad.

### Complex-System Coverage

- The complex-system blocks near the end of `run.sh` cover complex pockets, scan2d/scan3d, DMF mode, and TSOPT mode variants.
- Keep these as heavier regression checks after baseline tests pass.

### Recommended Order

1. Run baseline format checks (PDB/XYZ/GJF two-structure cases).
2. Run staged scan entry cases and inspect `scan/stage_*` outputs.
3. Run TS-focused cases and inspect imaginary modes/IRC behavior.
4. Run complex-system blocks last.

### Outputs To Inspect First

- `<run_name>.out`
- `<run_name>/summary.log`
- `<run_name>/summary.yaml`
- `<run_name>/path_search/mep_plot.png`
- `<run_name>/path_search/post_seg_*/`

### Context

- `test/test.md`: concise developer smoke-test note.
- This section: user-facing selection and interpretation guide.

## Related pages

- [Getting Started](getting-started.md)
- [Troubleshooting](troubleshooting.md)
- [Concepts](concepts.md)
