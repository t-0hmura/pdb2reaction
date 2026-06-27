---
name: pdb2reaction-workflows-output
description: Output parsing and multi-step workflow selection for pdb2reaction — `summary.json` schema, `seg_NN/` layout, R/TS/P/IM canonical paths, `bond_changes` interpretation, and the cluster + 1-step / multistep / scan-list / endpoint-MEP / TS-only / DFT//MLIP / stage-by-stage (subcommand-only, gate each stage) recipes plus energy-diagram extraction. TRIGGER on output parsing (`summary.json`, `result.json`, `seg_NN/`), extracting barriers / ΔE / Gibbs for a paper, choosing between multi-input / scan-list / endpoint-MEP / TS-only modes, or running the pipeline subcommand-by-subcommand with a success check at each stage (instead of one `all` run). SKIP for single-subcommand syntax (CLI skill) or install / HPC questions.
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

Result: `result_mep/segments/seg_NN/{reactant,ts,product}.pdb`,
`summary.json["segments"][0]["barrier_kcal"]`.

### 2. Multi-step recursive (multi-input MEP, recursive segmentation)

You have R and P. The mechanism is multi-step, but you don't have
intermediates handy. With `--refine-path True`, path-search recursively
segments by detecting bond changes (the default single-pass `path-opt`
does not).

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb \
    -c '...' -l '...' \
    --refine-path True \
    --tsopt --thermo \
    -o result_mep
```

With `--refine-path True`, the output `summary.json["n_segments"]` may be
> 1 — that's the recursion finding intermediates the inputs didn't contain.

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

With `--refine-path True`, recursive segmentation still runs *between*
adjacent endpoints, so you don't have to provide every elementary step;
the default single-pass `path-opt` assumes each adjacent pair is one
elementary step.

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

## Stage-by-stage execution (subcommand-only, gate each stage)

Run the pipeline as separate subcommands instead of one `pdb2reaction all` when you want
to **judge each stage's success before spending GPU time on the next** — e.g. confirm
path-search found the right segments / bond changes before optimizing a TS, or validate
the TS (one imaginary mode + correct IRC connectivity) before thermo / DFT. By default
`pdb2reaction all` runs this chain (the MEP stage is single-pass `path-opt`; pass
`--refine-path True` to swap in recursive `path-search`):

```
[extract] → path-opt → (per reactive seg) tsopt → freq → irc → [freq R/TS/P] → [dft] → energy-diagram
```

pdb2reaction runs the whole pipeline on a cluster model (the active-site cluster)
with a single MLIP backend, so pass the *same* `-l` / `-q` / `-m`
(and `--solvent` if used) on every stage. After each stage, read
its `result.json` / `summary.json` `status` and gate before continuing.

**Stage 0 — prep** (optional): if you start from full structures, cut the active-site
cluster with `-c 'RES,...' -r 2.6` (or pre-extract; see
[`pdb2reaction-cli/extract.md`](../pdb2reaction-cli/extract.md)). Most staged campaigns
start from already-prepared R/P cluster PDBs.

**Stage 1 — MEP (`path-search`)**

```bash
pdb2reaction path-search -i 1.R.pdb 3.P.pdb -l 'SAM:1,GPP:-3' -b uma -o ps/
```

**GATE** (`ps/summary.json`): `status == "success"`; inspect `n_segments` and EACH
segment's `bond_changes` — the intended bonds must be *formed AND broken* for the right
atoms ([`pdb2reaction-cli/bond-summary.md`](../pdb2reaction-cli/bond-summary.md),
`pdb2reaction-ts-strategy/SKILL.md`). Wrong segmentation / spurious changes → fix
chemistry or inputs **before** any TS work (don't optimize a TS for the wrong step).

**Stage 2 — per reactive segment: TS → validate → connectivity** (seed = `ps/hei_seg_NN.xyz`)

```bash
pdb2reaction tsopt -i ps/hei_seg_NN.xyz            -l 'SAM:1,GPP:-3' -b uma -o seg_NN/ts
pdb2reaction freq  -i seg_NN/ts/final_geometry.xyz -l 'SAM:1,GPP:-3' -b uma -o seg_NN/freq
pdb2reaction irc   -i seg_NN/ts/final_geometry.xyz -l 'SAM:1,GPP:-3' -b uma -o seg_NN/irc
```

**GATE** in order: tsopt `result.json` `status` is `converged` (or `completed`; not
`not_converged`) → freq `result.json` `n_imaginary == 1` (exactly one imaginary frequency)
whose mode moves the reacting atoms (0 or >1 → fix via fp64 / `--coord-type dlc` /
`--flatten`, see `pdb2reaction-ts-strategy/SKILL.md` §3, before trusting the barrier) → irc
`result.json` `status == "completed"` and forward/backward endpoints connect the **intended**
R and P (bond changes match this step). A TS that fails any gate is not this elementary step.

**Stage 3 — thermochemistry** (optional, = `all --thermo`): run `pdb2reaction freq` on
R / TS / P for the Gibbs/QRRHO profile.

**Stage 4 — DFT//MLIP** (optional, = `all --dft`):

```bash
pdb2reaction dft -i seg_NN/reactant.pdb -l 'SAM:1,GPP:-3' --func-basis 'wb97m-v/def2-tzvpd' -o seg_NN/dft/R   # repeat for ts, product
```

**GATE**: each `dft/<state>/result.{json,yaml}` shows `converged: true`.

**Stage 5 — energy diagram**

```bash
pdb2reaction energy-diagram -i "[0.0, 21.5, -0.7]" --label-x "['R','TS','P']" -o diagram.png
```

> Resuming after a walltime hit uses these same commands — see
> [`pdb2reaction-cli/all.md`](../pdb2reaction-cli/all.md) "Re-running individual stages".
> On any non-success status read `summary.log`, then `segments/seg_NN/<stage>/result.json`,
> before retrying. Large IRC/freq (n ≳ 4000 atoms) on a 16 GB GPU can OOM — run pysisyphus
> `bofill_update` on CPU.

## Output parsing

Schema, per-segment / post-segment keys, R/TS/P canonical paths,
programmatic extraction snippets, `bond_changes` interpretation, and
failed-run diagnostics live in [`summary-json.md`](summary-json.md).

## Energy diagrams

`pdb2reaction all` writes:

| Path | When | Content |
|---|---|---|
| `<out_dir>/segments/seg_NN/energy_diagram_UMA.png` | always | per-segment MLIP |
| `<out_dir>/segments/seg_NN/energy_diagram_G_UMA.png` | `--thermo` | + QRRHO Gibbs |
| `<out_dir>/segments/seg_NN/energy_diagram_DFT.png` | `--dft` | DFT//MLIP electronic |
| `<out_dir>/segments/seg_NN/energy_diagram_G_DFT_plus_UMA.png` | `--dft --thermo` | DFT//MLIP + Gibbs |
| `<out_dir>/energy_diagram_MEP.png` | always | bare MEP energies (promoted to root) |
| `<out_dir>/energy_diagram_UMA_all.png` | always | aggregated MLIP |
| `<out_dir>/energy_diagram_G_UMA_all.png` | `--thermo` | aggregated + Gibbs |
| `<out_dir>/energy_diagram_DFT_all.png` | `--dft` | aggregated DFT |
| `<out_dir>/energy_diagram_G_DFT_plus_UMA_all.png` | `--dft --thermo` | aggregated DFT + Gibbs |

In TS-only mode the per-segment diagrams land under `segments/seg_01/`.

To compose a custom diagram from energies of multiple runs, use
[`pdb2reaction-cli/energy-diagram.md`](../pdb2reaction-cli/energy-diagram.md):

```bash
pdb2reaction energy-diagram \
    -i "[0.0, 21.5, -0.7, 2.2, -18.2]" \
    --label-x "['R','TS1','IM','TS2','P']" \
    -o my_diagram.png
```

## See also
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
