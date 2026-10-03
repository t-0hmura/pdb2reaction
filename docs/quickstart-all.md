# Quickstart: `pdb2reaction all`

## Overview

`pdb2reaction all` builds a reaction path from the reactant (R) and product (P) in one run. It cuts out a cluster model around the substrates and searches the minimum energy path (MEP) between R and P. With `--tsopt --thermo --dft`, the same run continues to transition-state (TS) optimization, an intrinsic reaction coordinate (IRC) calculation, frequencies, and DFT single points.

The commands below use the bundled example of the GPP (geranyl pyrophosphate) C6-methyltransferase BezA in [`examples/`](https://github.com/t-0hmura/pdb2reaction/tree/main/examples): `1.R.pdb` is the reactant and `3.P.pdb` the product. Get it with `git clone https://github.com/t-0hmura/pdb2reaction && cd pdb2reaction/examples`. For your own reaction, replace them with your full-system structures.

### What it is for

* **A first run of the whole workflow**: run every stage once on the bundled example.
* **The MEP between R and P**: get the path and its highest-energy image (HEI), the TS candidate.
* **TS, IRC, frequencies, and DFT in the same run**: add `--tsopt --thermo --dft` to check the TS candidate.

## Minimal command

Give R and P in reaction order, the residues to cut out around (`-c`), and the ligand charges (`-l`).

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 --out-dir ./result_all
```

The run succeeded when the `====== Pipeline summary ======` block near the end of the console shows `Scientific status: success`; `summary.json` holds the same value in `scientific_status`.

### (Optional) Add post-processing in the same run

`--tsopt` adds TS optimization and IRC for each [reactive segment](glossary.md) (here `seg_01`), `--thermo` adds frequencies and thermochemistry, and `--dft` adds DFT single points on R, TS, and P. `--thermo` and `--dft` require `--tsopt`.

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 --tsopt --thermo --dft --out-dir ./result_all
```

## Before you run

The structures need every hydrogen atom, and R and P must list the same atoms in the same order; see [Before you run: the input structures](getting-started.md#before-you-run-the-input-structures).

## Output files

The minimal command writes:

```text
result_all/
├── summary.log                  # Run summary
├── summary.json                 # Results, with scientific_status
├── mep_trj.pdb                  # MEP over all segments
├── energy_diagram_MEP.png       # MEP energy profile over all segments
└── _work/                       # Intermediate files, including the HEI (TS candidate); kept after the run
    └── path_opt/                # MEP search (path_search/ with --refine-path, the recursive MEP search)
        ├── hei_seg_01.{xyz,pdb} # Highest-energy image of segment 1
        └── summary.json         # MEP search results
```

The minimal command stops after the MEP search and does not create `segments/`. With `--tsopt`, a reactive segment adds `segments/seg_NN/` with the R/TS/P structures (`reactant.pdb`, `ts.pdb`, `product.pdb`), `ts/`, and `irc/`; `--thermo` also adds `freq/`.

## Checking the result

1. **Completion**: `scientific_status` is `success` when every requested stage converged; otherwise it is `partial` or `failed`, with the [reasons](json-output.md#execution-and-requested-stage-completion) in `scientific_status_reasons`. With `--tsopt`, two checks are left for you: that the imaginary mode moves the bonds that form or break, and that the endpoints are the intended R and P.
2. **TS candidate**: open `_work/path_opt/hei_seg_01.pdb`, the HEI of the first segment. With `--tsopt`, also open the optimized TS, `segments/seg_01/ts.pdb`.
3. **Energy profile**: `energy_diagram_MEP.png` should show a clear barrier between R and P.
4. **TS (with `--tsopt`)**: a successful TS optimization gives one imaginary mode along the reaction coordinate. The console then prints `[tsopt] Converged (n_imag=1).`, and `summary.json` records the count in `post_segments[].tsopt.n_imaginary_modes`. Open `segments/seg_01/ts/vib/imag_*_trj.xyz` in a viewer and check that the mode moves the bonds that form or break.
5. **Endpoints (with `--tsopt`)**: open `segments/seg_01/irc/finished_irc_trj.xyz` and the optimized endpoints `segments/seg_01/reactant.pdb` and `product.pdb`, and check that they are the intended R and P. Even if the IRC does not converge, the result is usable when the endpoint optimizations reach the intended R and P.

For how `all` judges each stage, see [Reading the run status](all.md#reading-the-run-status).

## Notes

* **DFT and GPU memory**: `--dft` needs the DFT extra; for its installation and GPU memory, see the Notes of [Refine an MLIP TS with DFT](dft-backend.md#notes).
* **Barriers in `summary.json`**: `segments[].barrier_kcal` is the barrier on the MEP, before TS optimization. With `--tsopt`, `post_segments[].mlip.barrier_kcal` is the barrier from the optimized TS and endpoints; `--thermo` adds `post_segments[].gibbs_mlip.barrier_kcal` and `--dft` adds `post_segments[].dft.barrier_kcal`.
* **`rate_limiting_step`**: `rate_limiting_step.barrier_kcal` is the highest barrier among the segments, compared at the highest level that every segment has (`DFT//MLIP_Gibbs` > `DFT` > `MLIP_Gibbs` > `MLIP` > `MEP`); `rate_limiting_step.method` names that level.

## Next steps

- [Quickstart: scan](quickstart-scan.md): start from one structure when there is no product structure
- [Quickstart: TS-only mode](quickstart-tsopt.md): optimize and check a TS candidate you already have
- [Building the cluster model](model-setup.md): trim the model, or extend it when residues are missing
- [Tips for studying reaction mechanisms](mechanism-tips.md): plan the calculations, and what to try when the TS search fails
- [Refine an MLIP TS with DFT](dft-backend.md): refine and check the TS with DFT
- [`all`](all.md): full option reference (also `pdb2reaction all --help-advanced`)
- [JSON Output Reference](json-output.md): the fields of `summary.json`
- [Troubleshooting](troubleshooting.md): find an error message or symptom and its fix
