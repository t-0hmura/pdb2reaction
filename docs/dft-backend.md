# Refine an MLIP TS with DFT

Once MLIP has found a reasonable pathway, pdb2reaction can take its TS straight into a DFT TS optimization. It runs the TS optimization → IRC → endpoint optimization → frequency workflow with GPU-accelerated DFT through GPU4PySCF.

The MLIP pathway search is the main tool; DFT is an add-on that checks the TS candidate you found with MLIP.

---

## What it is for

- **Refine the TS at the DFT level**: run TS optimization → IRC → endpoint optimization → frequencies with DFT (`-b dft`).
- **Add DFT energies to an MLIP run**: run DFT single points on the MLIP R, TS, and P (`--dft`).

## Workflow

1. **Explore with MLIP**: generate pathways, try variants, and pick the most promising TS candidate.
2. **Refine with DFT**: run TS-only mode on that TS with `-b dft`. The command is example 2 in [Examples](#examples).
3. **Check**: as in an MLIP run, a successful TS optimization gives one imaginary mode along the reaction coordinate. `Scientific status: success` under `====== Pipeline summary ======` shows that every requested stage converged. Then check the mode and the IRC endpoints as in [Checking the result](quickstart-tsopt.md#checking-the-result).

## Examples

### 1. Search on a small model with MLIP

Run the MLIP search on a model small enough for DFT. The TS it writes, `result_mlip/segments/seg_01/ts.pdb`, is the input of example 2. `seg_01` is the first reaction segment; if there are several, pick the one for the step you want to refine.

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -r 0 --selected-resn '44,63,186' --tsopt -o ./result_mlip
```

[Make the model smaller](model-setup.md#make-the-model-smaller) explains how `-r 0 --selected-resn` builds this model.

### 2. Refine the TS with DFT

Pass the TS from example 1 as the only input, which selects TS-only mode; this needs the DFT extra.

```bash
pdb2reaction all -i result_mlip/segments/seg_01/ts.pdb -l 'SAM:1,GPP:-3' \
    --tsopt --thermo -b dft -o ./result_dft
```

## Keep the model under about 300 atoms

DFT optimization is practical up to roughly 300 atoms. Compare 300 with the atom count including the cap hydrogens, N + M from the console lines in [Check the model](model-setup.md#check-the-model). If you build a model of that size before the MLIP search, the TS it gives can go to DFT as is.

For ways to trim the model, see [Make the model smaller](model-setup.md#make-the-model-smaller); to fix the boundary atoms of a model you trimmed by hand, see {ref}`Freeze atoms and restrain distances <freeze-atoms-and-restraints>`.

## `-b dft` and `--dft`

| Option | What DFT computes | When to use it |
|---|---|---|
| `-b dft` | Every calculation in the run (MEP search, TS optimization, IRC, endpoint optimization, frequencies) | Refine and check a TS candidate at the DFT level |
| `--dft` | Single points on the R, TS, and P from the MLIP run (`all` only) | Get DFT energies on MLIP geometries |

`-b dft` works in 11 commands: `all`, `opt`, `tsopt`, `irc`, `freq`, `scan`, `scan2d`, `scan3d`, `path-opt`, `path-search`, and `sp`.

## Output files

With `-b dft`, the output has the same layout as an MLIP run in the same mode; for TS-only mode, see [Expected output](quickstart-tsopt.md#expected-output). `--dft` adds these files:

| File | Content |
|---|---|
| `segments/seg_NN/dft/{R,TS,P}/result.yaml` | DFT single-point result for each state |
| `segments/seg_NN/energy_diagram_DFT.png` | DFT energy diagram on the MLIP geometries |
| `segments/seg_NN/energy_diagram_G_DFT_plus_MLIP.png` | DFT energy plus the MLIP thermal correction (with `--thermo`) |
| `energy_diagram_DFT_all.png`, `energy_diagram_G_DFT_plus_MLIP_all.png` | The same diagrams over all segments, at the top of the output directory |

## Main options

| Option | Description | Default |
|---|---|---|
| `-b, --backend dft` | Use DFT as the calculator (GPU4PySCF; CPU PySCF with `--dft-engine cpu`). | `uma` |
| `--func-basis TEXT` | Functional and basis as `FUNCTIONAL/BASIS`. Applies to both `-b dft` and `--dft`. | `wb97m-v/def2-svp` |
| `--solvent TEXT`, `--solvent-model [pcm\|smd]` | Implicit solvent for `-b dft`. A solvent name selects SMD; add `--solvent-model pcm` for PCM. | `none` (gas phase) |
| `--dft/--no-dft` | Add DFT single points on R, TS, and P (`all` only). | `--no-dft` |
| `--dft-solvent TEXT`, `--dft-solvent-model [pcm\|smd]` | Implicit solvent for the `--dft` single points only; without `--dft`, the run stops with an error. | `none`, `smd` |

The other DFT options are listed in the [`all` reference](reference/commands/all.md).

> **Note:** In YAML, the `-b dft` settings go under `calc.dft`, and the `--dft` single points read the top-level `dft` section, as the `dft` command does. In either block, `pyscf` passes attributes to PySCF objects by name, for example `mf: {level_shift: 0.2}` for a hard-to-converge SCF (`calc.dft.pyscf` or `dft.pyscf`). See [YAML Reference](yaml-reference.md#calc) and the {ref}`dft section <dft-section>`.

## Notes

- **DFT extra**: install it with `pip install "pdb2reaction[dft]"` for the `cu130` or `cu132` PyTorch wheel, or with `pip install "pdb2reaction[dft-cuda12]"` for `cu126`. Without a GPU, add `--dft-engine cpu`.
- **Charge**: removing residues changes the total charge. Check `Total active site model charge` in the console output of example 1 before the DFT run.
- **Cap hydrogens**: `ts.pdb` keeps the cap hydrogens (`LKH`/`HL`) of the model, so their parent atoms stay frozen in the DFT run as in the MLIP run.
- **Combinations**: `-b dft` and `--dft` cannot be used together, and the run stops with an error at startup. To add DFT single points after a `-b dft` run, run `pdb2reaction sp -b dft` or `pdb2reaction dft` as a separate job. `--dft` and `--thermo` require `--tsopt`.
- **File and key names**: with `-b dft`, the diagrams keep the file names `energy_diagram_MLIP.png` and `energy_diagram_G_MLIP.png` (with `--thermo`), and the `summary.json` blocks keep the names `mlip` and `gibbs_mlip`. Both hold the DFT values, and the plot title shows DFT.
- **Memory and threads**: `-b dft` and `--dft` use the same low-memory mode, `--dft-nprocs`, and `--dft-memory` as [`pdb2reaction dft`](dft.md#main-options). If GPU memory runs out, trim the model; if `--dft` runs out of memory, drop `--dft` and run `pdb2reaction dft` separately.
- **SCF not converged**: an SCF that does not converge even from a fresh guess stops the run with `PySCF SCF did not converge with either the reused density or a fresh guess.`. Density fitting (`--no-dft-low-memory`), when memory allows, or a [YAML level shift](#main-options) can help it converge.
- **Stepwise grid**: `--scf-stepwise-grid` (on by default) first converges the first SCF of a run on grid level 1 with a tolerance of 1e-6, then on the requested grid and tolerance from that density; later SCFs start from the previous density as usual. The gain depends on the system; small systems can get slightly slower, and `--no-scf-stepwise-grid` turns it off. It is skipped when the requested grid level is 1 or lower, or when `pyscf.grids.atom_grid` sets the grid directly. If the coarse stage does not converge, the normal SCF runs.
- **SCF checkpoints**: they are off by default because the files can be very large. `--save-scf-checkpoint` saves them, and `--scf-checkpoint PATH` chooses the file. Without a path, commands other than `all` write `<out-dir>/_work/dft_scf/state.chk`, and `all` saves one file per state. A checkpoint is used only when its method, atom order, and coordinates match the current structure.

## See also

- [`dft`](dft.md): DFT single point with population analysis
- [`sp`](sp.md): single-point energy and forces with any backend
- [Quickstart: TS-only mode](quickstart-tsopt.md): check a TS candidate with `all --tsopt`
- [Building the cluster model](model-setup.md): build, trim, and extend the active-site model
- [Installation](installation.md#step-by-step-installation): step 7 installs the DFT extra
- [MLIP Backends](backends.md): choosing a backend
- [Troubleshooting](troubleshooting.md): what to do when a run fails
