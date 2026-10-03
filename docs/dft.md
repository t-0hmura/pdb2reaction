# `dft` (DFT single point)

## Overview

`dft` runs a **DFT (density functional theory) single point** on one structure with GPU4PySCF (GPU) or PySCF (CPU). It reports the **energy** and **atomic charges** from Mulliken, meta-Löwdin, and IAO (intrinsic atomic orbital) population analyses. The default method is ωB97M-V/def2-svp. It needs the DFT extra; see step 7 of [Step-by-step installation](installation.md#step-by-step-installation).

`all --dft` runs DFT single points on the reactant (R), transition state (TS), and product (P) of an MLIP (machine-learning interatomic potential) run. `-b dft` uses DFT for every calculation of a command; see [Refine an MLIP TS with DFT](dft-backend.md).

### What it is for

* **DFT energies on MLIP geometries**: single points on the R, TS, and P optimized with an MLIP.
* **Charge distribution**: per-atom charges, and spin densities for open shells.
* **Energies in implicit solvent**: PySCF's PCM (polarizable continuum model) or SMD (solvation model based on density) with `--solvent`.

---

## Examples

### 1. GPU single point

Compute the energy and charges of a neutral singlet on the GPU.

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 --out-dir ./result_dft
```

The console prints `E_total (Hartree): …` and `E_total (kcal/mol): …`, and `result_dft/result.yaml` has `energy.converged: true`.

### 2. Tighter SCF and a larger basis

Tighten the SCF (self-consistent field) and use a larger basis.

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 \
  --func-basis 'wb97m-v/def2-tzvpd' --scf-tol 1e-10 --scf-max-cycles 200 \
  --out-dir ./result_dft_tight
```

### 3. CPU only

Run with CPU PySCF on a machine without a GPU.

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 --dft-engine cpu --out-dir ./result_dft_cpu
```

### 4. Total charge from ligand charges

Without `-q`, `-l` gives the formal charges of the ligands, and `dft` adds the charges of the amino-acid residues and ions in the PDB to get the total charge; the console prints the breakdown.

```bash
pdb2reaction dft -i input.pdb -l 'SAM:1,GPP:-3' -m 1 --out-dir ./result_dft_ligand
```

---

## How it works

1. **Reading the structure**:
PDB, mmCIF, XYZ, and GJF inputs are read, and the coordinates passed to PySCF are saved as `input_geometry.xyz`. For an XYZ or GJF input, `--ref-pdb` gives the PDB/mmCIF topology that `-l` needs; `dft` writes no PDB, mmCIF, or GJF file.
2. **SCF**:
`--func-basis` sets the functional and basis, and `--dft-engine` selects GPU4PySCF (`gpu`, the default) or PySCF (`cpu`). A closed shell runs RKS and an open shell UKS. Low-memory mode, on by default, builds J and K directly without density fitting; on the GPU, a closed shell then uses GPU4PySCF's low-memory RKS. `--no-dft-low-memory` uses density fitting instead.
3. **Charges and the result file**:
After the SCF, `dft` computes Mulliken, meta-Löwdin, and IAO charges and spin densities and writes them with the energy (Hartree and kcal/mol) to `result.yaml`. An analysis that fails gives `null` in its column and a warning.

---

## Output files

`dft` writes these files to `--out-dir`:

```text
result_dft/
├─ input_geometry.xyz   # Geometry passed to PySCF
├─ result.yaml          # Energy, convergence, engine, per-atom charges and spin densities
├─ result.json          # Machine-readable summary (with --out-json)
└─ summary.json         # Copy of result.json; read result.json (with --out-json)
```

* **`energy`** in `result.yaml`: `hartree`, `kcal_per_mol`, `converged`, `used_gpu`, `used_lowmem`, and `engine`, which is `gpu4pyscf(rks_lowmem)`, `gpu4pyscf`, or `pyscf(cpu)`.
* **`charges [index, element, mulliken, lowdin, iao]`**: one row per atom; `index` starts at 0. The console prints the same table.
* **`spin_densities [index, element, mulliken, lowdin, iao]`**: the same layout. It is always written (all zeros for a closed shell), and the console prints it only for an open shell.
* **`result.json`** holds the charges and spin densities as `mulliken`, `lowdin`, and `iao` arrays, together with the charge, multiplicity, functional, basis, and SCF settings; see [JSON Output Reference](json-output.md#dft).

---

## Main options

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path | (required) | Input structure (`.pdb`, `.cif`, `.xyz`, `.gjf`, ...) |
| `-q, --charge` | integer | `None` | Total charge. Required unless `-l`, YAML `calc.charge`, or a `.gjf` input gives it |
| `-m, --multiplicity` | integer | `1` | Spin multiplicity (2S+1); a `.gjf` input supplies its own |
| `-l, --ligand-charge` | text | `None` | Per-residue formal charges (e.g. `'SAM:1,GPP:-3'`) or one total ligand charge. Needs PDB/mmCIF input or `--ref-pdb` |
| `--func-basis` | text | `wb97m-v/def2-svp` | Functional and basis as `FUNCTIONAL/BASIS` |
| `--scf-tol` | float | `1e-9` | SCF convergence threshold (Hartree) |
| `--scf-max-cycles` | integer | `100` | Maximum number of SCF iterations |
| `--dft-grid-level` | integer | `3` | Integration grid level (PySCF `grids.level`) |
| `--dft-engine` | `gpu` / `cpu` | `gpu` | GPU4PySCF or CPU PySCF |
| `--dft-low-memory/--no-dft-low-memory` | flag | `True` | Build J and K directly; `--no-dft-low-memory` uses density fitting |
| `--solvent` | text | `none` | Solvent name for PySCF PCM/SMD (e.g. `water`); `none` is the gas phase |
| `--solvent-model` | `pcm` / `smd` | `smd` | Implicit-solvent model |
| `--dft-nprocs` | integer | auto | PySCF CPU threads (detected from the scheduler and the host) |
| `--dft-memory` | text | auto | PySCF host RAM limit (e.g. `64GB`); this is not GPU memory |
| `-o, --out-dir` | path | `./result_dft/` | Output directory |

See the [generated CLI reference](reference/commands/dft.md) for every option.

> **Note:** In YAML (`--config`), the [`dft`](yaml-reference.md#dft-section) section holds the same settings. `dft.pyscf` passes attributes to PySCF objects by name, for example `pyscf: {mf: {level_shift: 0.2}}` for a hard-to-converge SCF. The charge and multiplicity go in `calc.charge` and `calc.spin`, as in the other commands; `calc.spin` is the multiplicity 2S+1, not PySCF's 2S. `-q`, `-l`, and `-m` given on the command line come first, then YAML, then a `.gjf` header.

---

## Notes

* **Basis cost**: `def2-tzvpd` costs much more than `def2-svp`. There is no fixed limit on atoms or GPU memory; the cost depends on the number of basis functions, the elements, the functional, the grid, and the GPU. Run one representative structure first and watch the peak memory; if it runs out, use a smaller basis or a GPU with more memory.
* **GPU**: if GPU4PySCF cannot run, `dft` stops with an error that suggests `--dft-engine cpu`; it does not switch to the CPU by itself. On a new GPU generation, an out-of-memory or unsupported-kernel error can come from the GPU4PySCF and CuPy versions rather than from memory, so check the versions and the traceback first.
* **CPU**: with `--dft-engine cpu`, how large a system is practical depends on the method and the machine, so time one representative single point.
* **Scratch space**: PySCF writes temporary files to `$PYSCF_TMPDIR`. On compute nodes with a small `/tmp`, point it to a disk with enough free space before the run.
* **Machines other than x86**: prebuilt GPU4PySCF wheels may not support them; build GPU4PySCF from [source](https://github.com/pyscf/gpu4pyscf).
* **Auxiliary basis**: with `--no-dft-low-memory`, PySCF uses its default auxiliary basis for the chosen basis; you need not give one.
* **IAO analysis** can fail on difficult systems.
* **Solvent**: `--solvent` here is PySCF's PCM or SMD, not the xTB solvent correction of the MLIP backends. `--solvent-model` accepts only lowercase `pcm` or `smd`. With PCM, a solvent name that `dft` does not know needs its dielectric constant in YAML, `dft.pyscf.with_solvent.eps`.
* **SCF not converged**: `dft` prints `WARNING: SCF did not converge to the requested tolerance.`, still writes `result.yaml` with `converged: false`, and exits with code 1. In low-memory mode it suggests retrying with `--no-lowmem`, the alias of `--no-dft-low-memory`, when memory allows.
* **Multiplicity** below 1 is rejected.
* **Earlier results**: a new run first removes `result.yaml`, `result.json`, and `summary.json` left in the output directory.
* **Exit codes**: see {ref}`Exit codes <exit-codes>`.

---

## See also

* [Refine an MLIP TS with DFT](dft-backend.md) — `-b dft` and `--dft` in a workflow, DFT settings, and GPU memory
* [sp](sp.md) — single-point energy and forces with any backend, including `-b dft`
* [all](all.md) — the full workflow; `--dft` adds DFT single points on R, TS, and P
* [MLIP Backends](backends.md) — choosing a backend
* [Troubleshooting](troubleshooting.md) — what to do when a run fails
