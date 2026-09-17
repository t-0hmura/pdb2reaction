# Quickstart: `pdb2reaction all --tsopt` (TS-only mode)

## Goal

Validate an existing TS candidate without the MEP stage. `pdb2reaction all --tsopt` runs `tsopt → irc`; `--thermo` adds `freq`, and `--dft` adds DFT single-points. PDB/mmCIF inputs are extracted only when `-c` is supplied.

## Prerequisites

- pdb2reaction installed (see [Installation](installation.md))
- One TS candidate geometry: PDB/mmCIF, XYZ, or GJF. PDB/mmCIF carries residue metadata.
- Charge: use `-q`, `-l`, a GJF header, or a configuration file; see [charge precedence](cli-conventions.md#charge-specification). Multiplicity comes from `-m`, then YAML `calc.spin`, then the GJF header, then `1`.
- TS-only mode requires one input, no `--scan-lists`, and `--tsopt`. Two or more inputs use the MEP route; one input with `--scan-lists` uses the scan route.

## Minimal command

```bash
pdb2reaction all -i ts_candidate.pdb -l 'SAM:1,GPP:-3' \
    --tsopt --thermo -o ./result_ts_only
```

For an XYZ singlet TS candidate, `-m` may be omitted:

```bash
pdb2reaction all -i ts_candidate.xyz -q -1 -b uma \
    --tsopt --thermo -o ./result_ts_only
```

`--tsopt` activates the validation chain; `--thermo` adds ZPE / Gibbs corrections from the freq stage. Both stages run on the same backend (UMA by default).

### (Optional) Add DFT single-points

```bash
pdb2reaction all -i ts_candidate.pdb -l 'SAM:1,GPP:-3' \
    --tsopt --thermo --dft --dft-func-basis 'wb97m-v/def2-tzvpd' \
    -o ./result_ts_only
```

> **VRAM warning:** `--dft` launches GPU4PySCF single-point jobs. Memory use
> depends on the structure, basis, functional, precision, and software stack;
> pilot a representative state and monitor peak memory on the target node. If
> it runs out of memory, drop `--dft` and run `pdb2reaction dft` separately
> with a smaller basis or trimmed cluster, or move the DFT step to a larger
> node. The `[dft]` extra must also be installed (see
> [Installation](installation.md)).

## Expected output

A successful run produces:

```text
result_ts_only/
├── summary.log                                # Run summary
├── summary.json                               # status: success | partial | failed
└── segments/
    └── seg_01/                                # TS-only deliverables
        ├── reactant.pdb                        # Canonical R/TS/P (TS-only mode)
        ├── ts.pdb
        ├── product.pdb
        ├── ts/
        │   ├── final_geometry.{xyz,pdb}
        │   └── vib/imag_*_trj.xyz             # Imaginary-mode trajectory
        ├── irc/
        │   └── {forward,backward,finished}_irc_trj.xyz
        ├── freq/{R,TS,P}/
        │   ├── frequencies_cm-1.txt
        │   └── thermoanalysis.yaml
        └── dft/{R,TS,P}/                      # --dft only
            └── result.yaml                    # always (when --dft)
```

## Inspecting the result

1. **Completion:** read `scientific_status` and `scientific_status_reasons` in `summary.json`. They report completion of required calculations and numerical optimizations. Frequency counts and chemical connectivity are separate checks.
2. **TS mode:** `post_segments[0].ts_imag.n_imag` should be `1` under the recorded criterion; `nu_imag_max_cm` gives its wavenumber. Visualize `segments/seg_01/ts/vib/imag_*_trj.xyz` and confirm the intended bond motion. Magnitude alone does not establish chemical relevance. The opt-in `irc.imag_below` filter (default `0.0` cm⁻¹; accepts `ν <= imag_below`) should be lowered only after a system-specific noise analysis.
3. **Connectivity and structures:** inspect `segments/seg_01/irc/finished_irc_trj.xyz`, the canonical `reactant.pdb`, `ts.pdb`, `product.pdb`, and `segments[0].bond_changes`. Check atom correspondence, geometry, and the expected bond changes. TS-only mode labels the higher-energy IRC endpoint as R; inspect `endpoint_assignment` before assigning chemical direction. See [all](all.md).
4. **Endpoint modes:** inspect `segments/seg_01/freq/{R,TS,P}/frequencies_cm-1.txt`. The complete signed spectrum is retained; `n_negative_modes` counts every negative sign. Residual R/P imaginary modes do not block thermochemistry but need inspection before identifying minima.
5. **Energies:** `rate_limiting_step.barrier_kcal` and `segments[0].delta_kcal` give ΔE‡ and ΔE. `post_segments[0].gibbs_mlip.barrier_kcal` / `.delta_kcal` give ΔG‡ and ΔG. The state-specific `thermoanalysis.yaml` records `electronic_energy_ha`, `zpe_correction_ha`, `sum_EE_and_ZPE_ha`, and `sum_EE_and_thermal_free_energy_ha` at the reported temperature and pressure (defaults: 298.15 K, 1 atm). Subtract the chosen R state from TS or P; do not infer chemical R/P identity from labels alone.

Energy differences are in kcal/mol; the state energies in `thermoanalysis.yaml` are in hartree.

| Observation | Next check |
|---|---|
| `n_imag == 0` | Improve the TS guess or MEP. A TS-only run without path information cannot identify the intended neighboring saddle; the default saddle-recovery budget is 0. |
| `n_imag >= 2` | Inspect every imaginary mode. Re-optimize with `all --thresh-post gau_tight` or `tsopt --thresh gau_tight`; `--flatten` is an explicit option for surplus modes. First-order classification requires one mode under the selected criterion. |
| Empty `bond_changes` or an unexpected endpoint | Inspect the TS mode and IRC; the path may connect the wrong wells. |
| Residual R/P imaginary modes | Inspect the endpoint geometry and modes; optionally tighten endpoint optimization or extend IRC. See [freq](freq.md). |

## Tips

- For finer control over `tsopt` parameters (`--opt-mode`, `--max-cycles`, Hessian options), run the standalone subcommand — see [tsopt](tsopt.md).
- Keep the default `FiniteDifference` unless analytical autograd has been validated for the chosen backend/model and system; its speed and memory cost are setup-dependent.
- Inspect the full option surface with `pdb2reaction all --help-advanced`.

## Next step

- Multi-structure MEP route: [Quickstart: `pdb2reaction all`](quickstart-all.md)
- Single-structure scan route: [Quickstart: `pdb2reaction all --scan-lists`](quickstart-scan.md)
- Full option references: [all](all.md), [tsopt](tsopt.md), [irc](irc.md), [freq](freq.md), [dft](dft.md)
