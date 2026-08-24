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

Result: `result_mep/segments/seg_NN/{reactant,ts,product}.pdb`; bridged
mmCIF/oversized-PDB inputs also receive `.cif` companions with original IDs.
Barriers come from `post_segments` (TSOPT+IRC-refined): ΔE‡ =
`summary.json["post_segments"][0]["mlip"]["barrier_kcal"]`, and with `--thermo`
ΔG‡ = `summary.json["post_segments"][0]["gibbs_mlip"]["barrier_kcal"]`.
Read top-level `mlip_backend`, `mlip_model`, and `mlip_precision` for the exact
provenance.
`summary.json["rate_limiting_step"]` is a legacy key giving the highest independently referenced local barrier with an explicit
`method` label; `summary.json["segments"][*]["barrier_kcal"]` is the raw MEP band
(also mirrored as `mep_barrier_kcal`) and is reported only under an MEP label.

### 2. Multi-step recursive (multi-input MEP, recursive segmentation)

You have R and P. The mechanism is multi-step, but you don't have
intermediates handy. With `--refine-path`, path-search recursively
segments by detecting bond changes (the default single-pass `path-opt`
does not).

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --refine-path \
    --tsopt --thermo \
    -o result_mep
```

With `--refine-path`, the output `summary.json["n_segments"]` may be
greater than 1: recursion has proposed additional path segments/intermediate
regions. That is a candidate decomposition, not proof of extra elementary
steps; validate each reactive HEI with TSOPT, frequency, and IRC.

### 3. Single-input scan-list (when only R is available)

You have just the reactant. Express the reaction as a sequence of
distance-restraint scans.

```bash
pdb2reaction all -i 1.R.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --scan-lists '[("CS1 SAM 320","C7 GPP 321",1.60)]' \
                 '[("H11 GPP 321","OE2 GLU 186",0.90)]' \
    --tsopt --thermo \
    -o result_scan
```

The flag occurs once; each literal value that follows it is one stage. See
[`pdb2reaction-cli/all-scan-list.md`](../pdb2reaction-cli/all-scan-list.md) for syntax details.

### 4. Endpoint-MEP with explicit intermediates

You have R, IM₁, IM₂, P from the literature. Pass them all in order:

```bash
pdb2reaction all -i 1.R.pdb 2.IM1.pdb 3.IM2.pdb 4.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    -o result_mep_4pt
```

With `--refine-path`, recursive segmentation runs across the ordered
endpoint series and may propose narrower candidate segments. The default
single-pass `path-opt` treats each adjacent pair as one MEP unit, but does not
prove that it contains exactly one elementary step.

### 5. TS-only validation (existing TS candidate)

You have a TS guess from another code or a prior run. Skip
extract / path-search:

```bash
pdb2reaction tsopt -i ts.xyz -q -1 -m 1 -b uma -o result_tsopt
pdb2reaction freq  -i result_tsopt/final_geometry.xyz -q -1 -m 1 -b uma -o result_freq
pdb2reaction irc   -i result_tsopt/final_geometry.xyz -q -1 -m 1 -b uma -o result_irc
```

Or run `pdb2reaction all` with a single `-i` **plus `--tsopt`** — that
combination selects TS-only mode (`pipeline_mode: "tsopt-only"` in the
summary). A single `-i` needs either `--tsopt` or `--scan-lists`; on its
own it is rejected with a `BadParameter` error. See
[`pdb2reaction-cli/all-ts-only.md`](../pdb2reaction-cli/all-ts-only.md).

### 6. DFT//MLIP single-point energies

After any of the above, evaluate R / TS / P with DFT single points. The
structures to feed `dft` are the deliverables of the preceding run, which live
at `<out_dir>/segments/seg_NN/{reactant,ts,product}.pdb` (and `.cif` for a
bridged input; below, `<out_dir>` is
`result_mep` from workflow 1):

```bash
LIGAND_CHARGE='SAM:1,GPP:-3'
FUNC_BASIS='wb97m-v/def2-tzvpd'
pdb2reaction dft -i result_mep/segments/seg_01/reactant.pdb \
    -l "$LIGAND_CHARGE" \
    --func-basis "$FUNC_BASIS" \
    --engine gpu \
    -o dft_R
pdb2reaction dft -i result_mep/segments/seg_01/ts.pdb      -l "$LIGAND_CHARGE" --func-basis "$FUNC_BASIS" -o dft_TS
pdb2reaction dft -i result_mep/segments/seg_01/product.pdb -l "$LIGAND_CHARGE" --func-basis "$FUNC_BASIS" -o dft_P
```

Composite the energies with `energy-diagram` (see below).

## Stage-by-stage execution (subcommand-only, gate each stage)

Run the pipeline as separate subcommands instead of one `pdb2reaction all` when you want
to **judge each stage's success before spending GPU time on the next** — e.g. confirm
path-search found the right segments / bond changes before optimizing a TS, or validate
the TS (one imaginary mode + correct IRC connectivity) before thermo / DFT. By default
`pdb2reaction all` runs this chain (the MEP stage is single-pass `path-opt`; pass
`--refine-path` to swap in recursive `path-search`):

```
[extract] → [scan] → MEP → (per candidate seg) tsopt → irc → endpoint opt → [freq R/TS/P] → [dft] → energy-diagram
```

pdb2reaction runs the whole pipeline on a cluster model (the active-site cluster)
with a single MLIP backend, so pass the *same* `-l` / `-q` / `-m` on every
stage.
After each stage, read
its `result.json` / `summary.json` `status` and gate before continuing.

**Stage 0 — prep** (optional): if you start from full structures, cut the active-site
cluster with `-c 'RES,...' -r 2.6` (or pre-extract; see
[`pdb2reaction-cli/extract.md`](../pdb2reaction-cli/extract.md)). Most staged campaigns
start from already-prepared R/P cluster PDBs.

**Stage 1 — default MEP (`path-opt`)**

```bash
pdb2reaction path-opt -i 1.R.pdb 3.P.pdb \
    -l 'SAM:1,GPP:-3' -b uma --out-json -o mep/
```

**GATE** (`mep/result.json`): require `status == "converged"` when the
engine exposes convergence; `completed` only means it returned. Inspect
`final_geometries_trj.xyz`, its energy profile, and `hei.pdb` / `hei.cif`.
Confirm endpoint atom identity/order and use `bond-summary` on the intended
endpoint pair before spending time on TS work.

For a known sequence R → IM → P, run one `path-opt` directory per adjacent
pair. If recursive candidate discovery is intentionally requested, replace
this stage with:

```bash
pdb2reaction path-search -i 1.R.pdb -i 3.P.pdb \
    -l 'SAM:1,GPP:-3' -b uma -o ps/
```

Then gate `ps/summary.json`, inspect every candidate segment's
`bond_changes`, and use `ps/hei_seg_NN.pdb` as that segment's seed. Recursive
segmentation can add unnecessary candidates on a poor/noisy path and is not
the default behavior of `all`.

**Stage 2 — TS candidate → frequency → IRC** (default seed = `mep/hei.pdb`)

```bash
pdb2reaction tsopt -i mep/hei.pdb -l 'SAM:1,GPP:-3' -b uma --out-json -o seg_01/ts
pdb2reaction freq  -i seg_01/ts/final_geometry.pdb -l 'SAM:1,GPP:-3' -b uma --out-json -o seg_01/freq_TS
pdb2reaction irc   -i seg_01/ts/final_geometry.pdb -l 'SAM:1,GPP:-3' -b uma --out-json -o seg_01/irc
```

Feed a residue-bearing **PDB or mmCIF** companion here, not a bare `.xyz`:
`-l/--ligand-charge` derives the total from residue names and is rejected on
bare `.xyz`/`.gjf`. Stages write an internal PDB companion whenever topology is
available; a bridged mmCIF/oversized-PDB run also writes CIF with original IDs.
If only `.xyz` remains, either pass verified total `-q`, or keep `-l` and add
`--ref-pdb <topology.pdb-or-cif>`.

Pass `--out-json` on each standalone subcommand: it writes the `result.json` that the gate
below reads (these subcommands default to `--no-out-json`).

**GATE** for first-order certification: tsopt numerical `optimization_status`
is `converged`; terminal `hessian_status` is `completed`; `saddle_validation` is
`first_order`; and the one imaginary displacement follows the reacting atoms.
The composite `all` workflow can also continue warning-labelled **diagnostic**
IRC from a numerically converged `higher_order` result when a validated negative
root exists, but that continuation is not first-order certification. Numerical
non-convergence, zero modes, failed/skipped PHVA, or no valid negative root stops
`all` after retaining TS artifacts. Then require standalone freq `result.json`
`n_imaginary == 1` before trusting the barrier → irc
`result.json` `status == "completed"` **plus** both enabled directions'
`*_converged` / `*_energy_increased` fields and frame counts are acceptable;
then confirm the first/last endpoints connect the intended R and P (bond
changes match this step). `completed` alone only means the runner returned.
A TS that fails any gate is not validated as this elementary step.

**Stage 3 — optimize and identify the raw IRC ends**

Standalone IRC writes `finished_first.xyz` and `finished_last.xyz`; it does
**not** create the canonical `reactant.pdb` / `product.pdb` files produced by
`all`. Optimize the two actual files and keep neutral names until identity is
verified:

```bash
pdb2reaction opt -i seg_01/irc/finished_first.xyz --ref-pdb 1.R.pdb \
    -l 'SAM:1,GPP:-3' -b uma --out-json -o seg_01/end_first
pdb2reaction opt -i seg_01/irc/finished_last.xyz --ref-pdb 1.R.pdb \
    -l 'SAM:1,GPP:-3' -b uma --out-json -o seg_01/end_last
```

**GATE**: both opt results must be `converged`. Compare bond patterns and
coordinates with the MEP left/right endpoints to assign which optimized end is
R and which is P. Do not infer chemical identity from `first` / `last` or from
energy alone. If you need the canonical
`segments/seg_NN/{reactant,ts,product}` layout and automatic orientation, use
`all` instead of manually inventing those paths.

**Stage 4 — thermochemistry** (optional, = `all --thermo`): run freq on the
actual TS and the two optimized endpoints:

```bash
pdb2reaction freq -i seg_01/end_first/final_geometry.pdb -l 'SAM:1,GPP:-3' -b uma --dump --out-json -o seg_01/freq_first
pdb2reaction freq -i seg_01/end_last/final_geometry.pdb  -l 'SAM:1,GPP:-3' -b uma --dump --out-json -o seg_01/freq_last
```

Use the identity assignment from Stage 3 when labeling R/P.

**GATE**: the requested thermochemistry fields for both endpoints must be
finite. R/P imaginary counts are optional minimum-certification diagnostics,
not a thermochemistry gate.

**Stage 5 — DFT//MLIP** (optional, = `all --dft`):

```bash
pdb2reaction dft -i seg_01/end_first/final_geometry.pdb -l 'SAM:1,GPP:-3' --func-basis 'wb97m-v/def2-tzvpd' --out-json -o seg_01/dft/first
pdb2reaction dft -i seg_01/ts/final_geometry.pdb        -l 'SAM:1,GPP:-3' --func-basis 'wb97m-v/def2-tzvpd' --out-json -o seg_01/dft/TS
pdb2reaction dft -i seg_01/end_last/final_geometry.pdb  -l 'SAM:1,GPP:-3' --func-basis 'wb97m-v/def2-tzvpd' --out-json -o seg_01/dft/last
```

**GATE**: each `dft/<state>/result.json` shows `converged: true`; label the two
endpoint energies using the Stage 3 identity assignment.

**Stage 6 — energy diagram**

```bash
pdb2reaction energy-diagram -i "[0.0, 21.5, -0.7]" --label-x "['R','TS','P']" -o diagram.png
```

> Resuming after a walltime hit uses these same commands — see
> [`pdb2reaction-cli/all.md`](../pdb2reaction-cli/all.md) "Re-running individual stages".
> On any non-success status read `summary.log`, then `segments/seg_NN/<stage>/result.json`,
> before retrying. Large dense Hessians can exceed GPU memory; PHVA for a
> justified frozen boundary or finite differences may reduce model/autograd
> peak memory, but the active-space Hessian itself remains dense.

## Output parsing

Schema, per-segment / post-segment keys, R/TS/P canonical paths,
programmatic extraction snippets, `bond_changes` interpretation, and
failed-run diagnostics live in [`summary-json.md`](summary-json.md).

## Energy diagrams

When the requested stages supply finite energies and PNG export succeeds,
`pdb2reaction all` writes:

| Path | When | Content |
|---|---|---|
| `<out_dir>/segments/seg_NN/energy_diagram_MLIP.png` | `--tsopt`, finite energies, successful PNG export | per-segment MLIP |
| `<out_dir>/segments/seg_NN/energy_diagram_G_MLIP.png` | `--thermo`, finite energies, successful PNG export | + QRRHO Gibbs |
| `<out_dir>/segments/seg_NN/energy_diagram_DFT.png` | `--dft`, finite energies, successful PNG export | DFT//MLIP electronic |
| `<out_dir>/segments/seg_NN/energy_diagram_G_DFT_plus_MLIP.png` | `--dft --thermo`, finite energies, successful PNG export | DFT//MLIP + Gibbs |
| `<out_dir>/energy_diagram_MEP.png` | MEP modes, finite energies, successful PNG export | bare MEP energies (promoted to root) |
| `<out_dir>/energy_diagram_MLIP_all.png` | complete MLIP segment energies and successful PNG export | aggregated MLIP |
| `<out_dir>/energy_diagram_G_MLIP_all.png` | complete Gibbs segment energies and successful PNG export | aggregated + Gibbs |
| `<out_dir>/energy_diagram_DFT_all.png` | complete DFT segment energies and successful PNG export | aggregated DFT |
| `<out_dir>/energy_diagram_G_DFT_plus_MLIP_all.png` | complete DFT//MLIP Gibbs segment energies and successful PNG export | aggregated DFT + Gibbs |

In TS-only mode the successfully exported per-segment diagrams land under
`segments/seg_01/`. An `energy_diagrams` metadata entry describes an attempted
or available diagram; verify the referenced path before treating the PNG as an
artifact.

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
