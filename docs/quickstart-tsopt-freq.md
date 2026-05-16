# Quickstart: `pdb2reaction all --tsopt` (TS-only mode)

## Goal

Validate an existing TS candidate end-to-end without running extract or path-search. `pdb2reaction all --tsopt` chains `tsopt → irc → freq → (dft)` from a single input structure and emits canonical reactant / TS / product geometries plus thermochemistry.

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
├── seg_01/                                    # Canonical R/TS/P (TS-only mode)
│   ├── reactant.pdb
│   ├── ts.pdb
│   └── product.pdb
└── tsopt_single/                              # tsopt → irc → freq orchestration
    ├── ts/
    │   ├── final_geometry.{xyz,pdb}
    │   └── vib/imag_*_trj.xyz                 # Imaginary-mode animation
    ├── irc/
    │   └── {forward,backward,finished}_irc_trj.xyz
    ├── freq/{R,TS,P}/
    │   ├── frequencies_cm-1.txt
    │   └── thermoanalysis.yaml
    └── dft/{R,TS,P}/                          # --dft only
        ├── result.yaml                        # always (when --dft)
        └── result.json                        # --out-json opt-in
```

**What to check (in execution order):**

1. `summary.json` — `status` is `"success"`; inspect `segments[0].barrier_kcal` and `segments[0].delta_kcal`.
2. `post_segments[0].ts_imag.n_imag` — should be exactly `1` for a first-order saddle.
3. `tsopt_single/irc/{forward,backward}_irc_trj.xyz` — open in PyMOL; trajectories should reach the expected reactant and product wells.
4. `segments[0].bond_changes` — non-empty dict listing bonds broken / formed along the IRC.
5. `tsopt_single/freq/{R,TS,P}/frequencies_cm-1.txt` — R and P should have no imaginary frequencies; TS exactly one.

**Troubleshoot:**

| Symptom | Likely cause | Fix |
|---|---|---|
| `post_segments[0].ts_imag.n_imag == 0` | TS guess collapsed to a minimum | Re-do the TS guess (e.g. via `path-search`); TS-only mode cannot recover a missing saddle |
| `n_imag >= 2` | Near-degenerate negative modes | Add `--flatten` to flatten extras; see [tsopt](tsopt.md) for `hessian_dimer.flatten_max_iter` |
| `segments[0].bond_changes == {}` or IRC reaches the wrong endpoint | Imaginary mode not along the intended coordinate, or TS connects two essentially identical wells | Visualize `tsopt_single/ts/vib/imag_*_trj.xyz` in PyMOL; if the mode is wrong, re-pick the TS guess |
| `freq/{R,P}/frequencies_cm-1.txt` shows residual imaginary modes | IRC endpoint is not a true minimum | Tighten convergence (`--thresh-post baker`) or extend IRC max cycles in YAML; see [freq](freq.md) |

## Tips

- For finer control over `tsopt` parameters (`--opt-mode`, `--max-cycles`, Hessian options), run the standalone subcommand — see [tsopt](tsopt.md).
- `--hessian-calc-mode Analytical` is recommended when VRAM permits (default: `FiniteDifference`).
- Inspect the full option surface with `pdb2reaction all --help-advanced`.

## Next step

- Multi-structure MEP route: [Quickstart: `pdb2reaction all`](quickstart-all.md)
- Single-structure scan route: [Quickstart: `pdb2reaction all --scan-lists`](quickstart-scan.md)
- Full option references: [all](all.md), [tsopt](tsopt.md), [irc](irc.md), [freq](freq.md), [dft](dft.md)
