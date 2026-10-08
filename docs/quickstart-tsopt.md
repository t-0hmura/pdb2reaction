# Quickstart: `pdb2reaction all --tsopt` (TS-only mode)

TS-only mode checks one transition-state (TS) candidate without a minimum energy path (MEP) search. `pdb2reaction all --tsopt` optimizes the TS, follows the intrinsic reaction coordinate (IRC) in both directions, and optimizes the two endpoints, the reactant (R) and the product (P). `--thermo` adds vibrational analysis and thermochemistry, and `--dft` adds DFT single points on R, TS, and P.

---

## What it is for

* **Refining a candidate from a scan or an MEP**: optimize the top of a scan or the highest-energy image (HEI) of an MEP into a TS.
* **Checking a candidate made another way**: confirm that a structure from another program, or one built by hand, is a TS (n_imag = 1) that connects the intended R and P.
* **Checking an MLIP TS before DFT**: confirm a TS from the machine-learning interatomic potential (MLIP) before you refine it with [DFT](dft-backend.md).

## Minimal command

Pass one TS candidate with `--tsopt`. The bundled examples have no TS candidate, so the command below uses the HEI from the [`all` quickstart](quickstart-all.md) run; for your own reaction, pass your own candidate.

```bash
pdb2reaction all -i result_all/_work/path_opt/hei_seg_01.pdb -l 'SAM:1,GPP:-3' \
    --tsopt --thermo -o ./result_ts_only
```

The run succeeded when the `====== Pipeline summary ======` block near the end of the console shows `Scientific status: success`; `summary.json` holds the same value in `scientific_status`.

For an XYZ file, such as a candidate from another program, give the total charge of the model with `-q` (0 for the model of the `all` quickstart).

```bash
pdb2reaction all -i ts_candidate.xyz -q 0 \
    --tsopt --thermo -o ./result_ts_only
```

### (Optional) Add DFT single-points

`--dft` adds DFT single points on R, TS, and P, and `--func-basis` sets the functional and basis (default `wb97m-v/def2-svp`).

```bash
pdb2reaction all -i result_all/_work/path_opt/hei_seg_01.pdb -l 'SAM:1,GPP:-3' \
    --tsopt --thermo --dft --func-basis 'wb97m-v/def2-tzvpd' \
    -o ./result_ts_only
```

## Before you run

* **Input**: one TS candidate as PDB/mmCIF, XYZ, or GJF. A cluster model is cut out only when `-c` is given; otherwise the structure is used as is.
* **Charge and multiplicity**: give the [charge](cli-conventions.md#charge-specification) with `-q`, `-l` (PDB/mmCIF only), YAML `calc.charge`, or a GJF header. The multiplicity comes from `-m`, then YAML `calc.spin`, then the GJF header, then 1.
* **When TS-only mode runs**: one input, `--tsopt`, and no `--scan-lists`. Two or more inputs run the [MEP search](quickstart-all.md), and one input with `--scan-lists` runs a [scan](quickstart-scan.md).

## Expected output

A successful run writes:

```text
result_ts_only/
├── summary.log                     # Run summary
├── summary.json                    # Results, with scientific_status
└── segments/
    └── seg_01/
        ├── reactant.pdb            # R/TS/P structures (.xyz for XYZ input, .gjf for GJF input)
        ├── ts.pdb
        ├── product.pdb
        ├── energy_diagram_MLIP.png # R–TS–P energy diagram (energy_diagram_G_MLIP.png with --thermo)
        ├── ts/
        │   ├── final_geometry.{xyz,pdb}
        │   └── vib/imag_*_trj.xyz  # Animation of each imaginary mode
        ├── irc/
        │   └── {forward,backward,finished}_irc_trj.xyz
        ├── freq/{R,TS,P}/          # --thermo
        │   ├── frequencies_cm-1.txt
        │   └── thermoanalysis.yaml
        └── dft/{R,TS,P}/           # --dft
            └── result.yaml
```

## Checking the result

1. **Completion**: `scientific_status` is `success` when every requested stage converged; otherwise it is `partial` or `failed`, with the [reasons](json-output.md#execution-and-requested-stage-completion) in `scientific_status_reasons`. Two checks are left for you: that the imaginary mode moves the bonds that form or break, and that the endpoints are the intended R and P.
2. **TS mode**: a successful TS optimization gives one imaginary mode along the reaction coordinate. The console then prints `[tsopt] Converged (n_imag=1).`, and `summary.json` records the count in `post_segments[0].tsopt.n_imaginary_modes`. Open `segments/seg_01/ts/vib/imag_*_trj.xyz` in a viewer and check that the mode moves the bonds that form or break.
3. **Endpoints**: open `segments/seg_01/irc/finished_irc_trj.xyz` and the R/TS/P structures (`reactant.pdb`, `ts.pdb`, `product.pdb`), and read `segments[0].bond_changes`. The endpoints should be the intended R and P. Even if the IRC does not converge, the result is usable when the endpoint optimizations reach the intended R and P.
4. **Endpoint frequencies**: with `--thermo`, `segments/seg_01/freq/{R,TS,P}/frequencies_cm-1.txt` lists every frequency with its sign. R and P should have no imaginary mode (no value below −5.00 cm⁻¹).

| Result | What to try |
|---|---|
| n_imag = 0 | Start from a better candidate, such as the HEI of an MEP or the top of a scan; TS-only mode has no path to guide it. |
| n_imag ≥ 2 | Watch every imaginary mode. Re-optimize with `--flatten`, or tighten convergence with `all --thresh-post gau_tight` (stricter than the default [`baker`](tsopt.md#how-it-works)) or `tsopt --thresh gau_tight`. |
| `bond_changes` is empty, or an endpoint is not the intended one | Check the TS mode and the IRC; the path may connect other minima. |
| R or P keeps an imaginary mode | Check the endpoint geometry and the mode. Tighten the endpoint optimization with `--thresh-post gau_tight`, or extend the IRC with `--irc-max-cycles` (default 125). |

If the TS is still not found, see {ref}`Check the TS <mechanism-check-ts>` and {ref}`When the TS search fails <ts-search-fails>`.

## Notes

* **Cap hydrogens with XYZ or GJF**: their parent atoms are frozen only when a PDB of the same atoms is given with `--ref-pdb`; see {ref}`Freeze atoms and restrain distances <freeze-atoms-and-restraints>`.
* **Hessian mode**: keep the default `--hessian-calc-mode FiniteDifference`. Set `--hessian-calc-mode Analytical` only after checking its speed, memory use, and results for your backend, model, and system.
* **Energies**: `post_segments[0].mlip.barrier_kcal` is ΔE‡ (TS − R) and `.delta_kcal` is ΔE (P − R) in kcal/mol, from the optimized TS and endpoints. Since there is no MEP, `segments[0].barrier_kcal` and `.delta_kcal` hold the same values. With `--thermo`, `post_segments[0].gibbs_mlip.barrier_kcal` and `.delta_kcal` give ΔG‡ and ΔG; with `--dft`, `post_segments[0].dft.barrier_kcal` and `.delta_kcal` give the DFT values.
* **R and P labels**: without an MEP the direction of the reaction is unknown, so TS-only mode labels the higher-energy IRC endpoint R and the lower one P, and records this rule in `endpoint_assignment` in `summary.json`. The labels are not the chemical direction; the barrier from P is `barrier_kcal − delta_kcal`.
* **Imaginary modes at R or P**: thermochemistry is computed even when R or P keeps an imaginary mode, and that mode is left out of ZPE and G. Check the mode before treating the structure as a minimum.
* **When IRC runs**: `all` goes on to IRC only when the TS optimization converged, its final Hessian finished, and n_imag ≥ 1. With n_imag ≥ 2 the IRC follows one imaginary mode as a diagnostic; it does not make the structure a first-order saddle point.
* **More `tsopt` settings**: for `--opt-mode`, `--max-cycles`, and Hessian options, run [`tsopt`](tsopt.md) on its own. `pdb2reaction all --help-advanced` lists every option of `all`.
* **DFT**: to optimize the TS itself with DFT, and for the DFT extra and GPU memory, see [Refine an MLIP TS with DFT](dft-backend.md).

## Next steps

- [Refine an MLIP TS with DFT](dft-backend.md): refine and check the TS with DFT
- [Tips for studying reaction mechanisms](mechanism-tips.md): check the TS, and what to try when the TS search fails
- [`tsopt`](tsopt.md), [`irc`](irc.md), [`freq`](freq.md): run each stage on its own
- [Quickstart: `pdb2reaction all`](quickstart-all.md): build an MEP from R and P
- [Quickstart: scan](quickstart-scan.md): build a path from one structure
- [`all`](all.md), [`dft`](dft.md): full option references
- [Troubleshooting](troubleshooting.md): find an error message or symptom and its fix
