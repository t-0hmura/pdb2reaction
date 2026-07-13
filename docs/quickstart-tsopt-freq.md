# Quickstart: `pdb2reaction all --tsopt` (TS-only mode)

## Goal

Validate an existing TS candidate end-to-end without running extract or the MEP (path-opt) stage. `pdb2reaction all --tsopt` chains `tsopt → irc → freq → (dft)` from a single input structure and emits canonical reactant / TS / product geometries and thermochemistry.

## Prerequisites

- pdb2reaction installed (see [Installation](installation.md))
- One TS candidate geometry: `.pdb` (preferred — carries residue / charge info) or `.xyz`
- Charge: exactly one of `-q/--charge INT`, `--ligand-charge/-l 'RES:Q,...'`, or a `.gjf` header. For `.xyz` inputs `-q` and `-m` are mandatory — pass `--ref-pdb cluster.pdb` if you want `-l 'RES:Q'` resolution.
- TS-only mode activates when **all three** hold: exactly one `-i` input, no `--scan-lists`, and `--tsopt` (or `--tsopt True`). Otherwise the CLI raises `BadParameter` at the input gate (`Provide at least two structures with -i/--input in reaction order, or use a single structure with --scan-lists, or a single structure with --tsopt True.`)

## Minimal command

```bash
pdb2reaction all -i ts_candidate.pdb -l 'SAM:1,GPP:-3' \
    --tsopt --thermo -o ./result_ts_only
```

For an XYZ TS candidate, supply `-q` and `-m` explicitly:

```bash
pdb2reaction all -i ts_candidate.xyz -q -1 -m 1 -b uma \
    --tsopt --thermo -o ./result_ts_only
```

`--tsopt` activates the validation chain; `--thermo` adds ZPE / Gibbs corrections from the freq stage. Both stages run on the same backend (UMA by default).

### (Optional) Add DFT single-points

```bash
pdb2reaction all -i ts_candidate.pdb -l 'SAM:1,GPP:-3' \
    --tsopt --thermo --dft --dft-func-basis 'wb97m-v/def2-tzvpd' \
    -o ./result_ts_only
```

> **VRAM warning:** `--dft` launches GPU4PySCF single-point jobs and can OOM on GPUs with < 24 GB VRAM for clusters above ~200 atoms. If you hit `CUDA out of memory`, drop `--dft` and run `pdb2reaction dft` separately with a smaller basis or trimmed cluster, or move the DFT step to a larger-VRAM node. The `[dft]` extra must also be installed (see [Installation](installation.md)).

## Expected output

A successful run produces:

```text
result_ts_only/
├── summary.log                                # Human-readable summary
├── summary.json                               # status: success | partial | failed
└── segments/
    └── seg_01/                                # TS-only deliverables
        ├── reactant.pdb                        # Canonical R/TS/P (TS-only mode)
        ├── ts.pdb
        ├── product.pdb
        ├── ts/
        │   ├── final_geometry.{xyz,pdb}
        │   └── vib/imag_*_trj.xyz             # Imaginary-mode animation
        ├── irc/
        │   └── {forward,backward,finished}_irc_trj.xyz
        ├── freq/{R,TS,P}/
        │   ├── frequencies_cm-1.txt
        │   └── thermoanalysis.yaml
        └── dft/{R,TS,P}/                      # --dft only
            ├── result.yaml                    # always (when --dft)
            └── result.json                    # --out-json opt-in
```

## Inspecting the result

Walk these in order; each step has a fast pass/fail check before you move on.

**1. Top-level verdict** — open `result_ts_only/summary.json`:

- `status` should be `"success"` (`"partial"` = stages ran but one validator flagged; `"failed"` = a stage errored).
- `rate_limiting_step.barrier_kcal` and `segments[0].delta_kcal` are the headline ΔE‡ and ΔE in kcal/mol.
- `post_segments[0].gibbs_uma.barrier_kcal` / `.delta_kcal` are the same numbers with ZPE + thermal corrections applied (ΔG‡, ΔG at 298.15 K, 1 atm).

**2. Imaginary mode at the saddle** — `post_segments[0].ts_imag`:

- `n_imag` must be exactly `1`. `nu_imag_max_cm` (negative cm⁻¹) is the imaginary wavenumber.
- A chemistry-relevant TS sits in **|ν| ≈ 100–2000 cm⁻¹**. If `|nu_imag_max_cm|` is below ~50 cm⁻¹, the "imag" mode is a soft delocalised motion, not a reactive coordinate — IRC will wander. Re-pick the TS guess, tighten `--thresh-post`, or raise `irc.imag_below` in the YAML (default `0.0` cm⁻¹) to skip very-soft modes when initializing IRC.
- Visualise the mode: `pymol result_ts_only/segments/seg_01/ts/vib/imag_*_trj.xyz` — the animation should swing precisely the bond(s) you expect to break/form, not show whole-molecule/residue tumbling.

**3. IRC connectivity** — open the IRC trajectory in PyMOL:

```bash
pymol result_ts_only/segments/seg_01/irc/finished_irc_trj.xyz
```

The merged trajectory (forward + backward) should land on the intended reactant and product wells. Cross-check `segments[0].bond_changes` in `summary.json`: a non-empty list of `Bond formed` / `Bond broken` entries with sensible distance changes (e.g. `C12-O14 : 3.7 Å --> 1.4 Å`) is the chemistry-level confirmation.

**4. Endpoint minima and thermochemistry** — for each of R, TS, P:

- `result_ts_only/segments/seg_01/freq/{R,TS,P}/frequencies_cm-1.txt` — R and P should list zero negative frequencies; TS exactly one (matching step 2).
- `result_ts_only/segments/seg_01/freq/{R,TS,P}/thermoanalysis.yaml` — fields are `electronic_energy_ha`, `zpe_correction_ha`, `sum_EE_and_ZPE_ha`, and `sum_EE_and_thermal_free_energy_ha` (the absolute Gibbs energy in hartree, at `temperature_K: 298.15`, `pressure_atm: 1.0`). Subtract R from TS for ΔG‡; p2r already reports the difference in `gibbs_uma.barrier_kcal`.

**5. Visual structure check** — load the canonical R/TS/P PDBs:

```bash
pymol result_ts_only/segments/seg_01/reactant.pdb result_ts_only/segments/seg_01/ts.pdb result_ts_only/segments/seg_01/product.pdb
```

In PyMOL: `align` the three states, label the reactive atoms (`label name C12+O14+C2+C17, name`), and confirm bond-length deltas match `bond_changes`.

**Troubleshoot:**

| Symptom | Likely cause | Fix |
|---|---|---|
| `post_segments[0].ts_imag.n_imag == 0` | TS guess collapsed to a minimum | Re-do the TS guess with `path-search`; `all` passes its MEP tangent internally for bounded recovery. An ordinary TS-only run without path information cannot identify the intended neighboring saddle |
| `n_imag >= 2` | Near-degenerate negative modes | Add `--flatten` to flatten extras; see [tsopt](tsopt.md) for `hessian_dimer.flatten_max_iter` |
| `segments[0].bond_changes` is empty (`""` or `"(no covalent changes detected)"`) or IRC reaches the wrong endpoint | Imaginary mode not along the intended coordinate, or TS connects two essentially identical wells | Visualize `segments/seg_01/ts/vib/imag_*_trj.xyz` in PyMOL; if the mode is wrong, re-pick the TS guess |
| `freq/{R,P}/frequencies_cm-1.txt` shows residual imaginary modes | IRC endpoint is not a true minimum | Tighten convergence (`--thresh-post baker`) or extend IRC max cycles in YAML; see [freq](freq.md) |

## Tips

- For finer control over `tsopt` parameters (`--opt-mode`, `--max-cycles`, Hessian options), run the standalone subcommand — see [tsopt](tsopt.md).
- `--hessian-calc-mode Analytical` is recommended when VRAM permits (default: `FiniteDifference`).
- Inspect the full option surface with `pdb2reaction all --help-advanced`.

## Next step

- Multi-structure MEP route: [Quickstart: `pdb2reaction all`](quickstart-all.md)
- Single-structure scan route: [Quickstart: `pdb2reaction all --scan-lists`](quickstart-scan.md)
- Full option references: [all](all.md), [tsopt](tsopt.md), [irc](irc.md), [freq](freq.md), [dft](dft.md)
