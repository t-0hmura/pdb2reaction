---
name: pdb2reaction-workflows-output
description: Output parsing and multi-step workflow selection for pdb2reaction — `summary.json` schema, `seg_NN/` layout, R/TS/P/IM canonical paths, `bond_changes` interpretation, and the cluster + 1-step / multistep / scan-list / endpoint-MEP / TS-only / DFT//MLIP recipes plus energy-diagram extraction. TRIGGER on output parsing (`summary.json`, `result.json`, `seg_NN/`), extracting barriers / ΔE / Gibbs for a paper, or choosing between multi-input / scan-list / endpoint-MEP / TS-only modes. SKIP for single-subcommand syntax (CLI skill) or install / HPC questions.
---

# pdb2reaction Workflows and Output Parsing

## Six canonical workflows

### 1. Cluster + 1-step reaction (multi-input MEP)

You have R and P PDBs (from a published QM study). One step.

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    -o result_mep
```

Result: `result_mep/seg_NN/{reactant,ts,product}.pdb`,
`summary.json["segments"][0]["barrier_kcal"]`.

### 2. Multi-step recursive (multi-input MEP, recursive segmentation)

You have R and P. The mechanism is multi-step, but you don't have
intermediates handy. Path-search recursively segments by detecting
bond changes.

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb \
    -c '...' -l '...' \
    --tsopt --thermo \
    -o result_mep
```

The output `summary.json["n_segments"]` may be > 1 — that's the
recursion finding intermediates the inputs didn't contain.

### 3. Single-input scan-list (when only R is available)

You have just the reactant. Articulate the reaction as a sequence of
distance-restraint scans.

```bash
pdb2reaction all -i 1.R.pdb \
    -c '...' -l '...' \
    --scan-lists '[("CS1 SAM 320","C7 GPP 321",1.60)]' \
                 '[("H11 GPP 321","OE2 GLU 186",0.90)]' \
    --tsopt --thermo \
    -o result_scan
```

Each `--scan-lists` argument is one stage. See
[`pdb2reaction-cli/all-scan-list.md`](../pdb2reaction-cli/all-scan-list.md) for syntax details.

### 4. Endpoint-MEP with explicit intermediates

You have R, IM₁, IM₂, P from the literature. Pass them all in order:

```bash
pdb2reaction all -i 1.R.pdb 2.IM1.pdb 3.IM2.pdb 4.P.pdb \
    -c '...' -l '...' \
    --tsopt --thermo \
    -o result_mep_4pt
```

Recursive segmentation still runs *between* adjacent endpoints, so
you don't have to provide every elementary step.

### 5. TS-only validation (existing TS candidate)

You have a TS guess from another code or a prior run. Skip
extract / path-search:

```bash
pdb2reaction tsopt -i ts.xyz -q -1 -m 1 -b uma -o result_tsopt
pdb2reaction freq  -i result_tsopt/final_geometry.xyz -q -1 -m 1 -b uma -o result_freq
pdb2reaction irc   -i result_tsopt/final_geometry.xyz -q -1 -m 1 -b uma -o result_irc
```

Or use `pdb2reaction all` with a single `-i` (collapses to TS-only
mode automatically; see [`pdb2reaction-cli/all-ts-only.md`](../pdb2reaction-cli/all-ts-only.md)).

### 6. DFT//MLIP refinement

After any of the above, refine R / TS / P energies at DFT level:

```bash
pdb2reaction dft -i seg_01/reactant.pdb \
    -l 'SAM:1,GPP:-3' \
    --func-basis 'wb97m-v/def2-tzvpd' \
    --engine gpu \
    -o dft_R
pdb2reaction dft -i seg_01/ts.pdb       -l '...' --func-basis '...' -o dft_TS
pdb2reaction dft -i seg_01/product.pdb  -l '...' --func-basis '...' -o dft_P
```

Composite the energies with `energy-diagram` (see below).

## Output parsing

Schema, per-segment / post-segment keys, R/TS/P canonical paths,
programmatic extraction snippets, `bond_changes` interpretation, and
failed-run diagnostics live in [`summary-json.md`](summary-json.md).

## Energy diagrams

`pdb2reaction all` writes per-segment diagrams under each `seg_NN/`
(e.g. `seg_01/tsopt/energy_diagram_UMA.png`,
`energy_diagram_G_UMA.png`, `energy_diagram_DFT.png` when `--dft`,
`energy_diagram_G_DFT_plus_UMA.png` when `--dft --thermo`) and writes
the **aggregated** multi-segment diagrams at the top of the output
directory:

- `<out_dir>/energy_diagram_MEP.png` — bare MEP energies from the
  path-search string (MLIP, no thermochemistry; bundled with
  `path_search/` when `--refine-path` true, `path_opt/` otherwise).
- `<out_dir>/energy_diagram_UMA_all.png` — per-segment energies stitched
  across all `seg_NN/` for whichever backend was used (MLIP).
- `<out_dir>/energy_diagram_G_UMA_all.png` — same with QRRHO Gibbs
  thermochemistry (when `--thermo`).
- `<out_dir>/energy_diagram_G_DFT_plus_UMA_all.png` — DFT//MLIP combined
  diagram (when `--dft --thermo`).

To compose a custom diagram from energies of multiple runs, use
[`pdb2reaction-cli/energy-diagram.md`](../pdb2reaction-cli/energy-diagram.md):

```bash
pdb2reaction energy-diagram \
    -i 0.0 21.5 -0.7 2.2 -18.2 \
    --label-x R TS1 IM TS2 P \
    -o my_diagram.png
```

## Cross-references

- [`summary-json.md`](summary-json.md) — `summary.json` schema, R/TS/P
  paths, programmatic extraction, bond-change interpretation, failed-run
  diagnostics.
- [`pdb2reaction-cli/all.md`](../pdb2reaction-cli/all.md) and the three `all-*.md` mode files.
- `pdb2reaction-cli/{tsopt,freq,irc,dft}.md` — per-stage
  `result.json` schemas.
- [`pdb2reaction-cli/bond-summary.md`](../pdb2reaction-cli/bond-summary.md) — same bond-change algorithm,
  standalone.
- `pdb2reaction-structure-io/SKILL.md` — input file formats that feed
  these workflows.
