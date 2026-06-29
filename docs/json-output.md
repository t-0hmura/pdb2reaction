# JSON Output Reference

pdb2reaction emits machine-readable JSON for programmatic use by scripts and agents.

## `--out-json` flag

Every MLIP-based subcommand supports `--out-json / --no-out-json` (default: off).
When enabled, a `result.json` file is written to the output directory alongside the normal outputs.

```bash
pdb2reaction opt -i r.pdb -q -1 --out-json --out-dir result_opt
cat result_opt/result.json | python -m json.tool
```

The `all` and `path-search` commands always write `summary.json` (no `--out-json` flag needed).

### `summary.json` mirror

`write_result_json` mirrors every per-stage `result.json` payload to `summary.json` alongside it. Agent scripts can read a single filename (`summary.json`) across every subcommand; the `result.json` written next to it contains the same payload.

## Common envelope

Every `result.json` (and the mirrored `summary.json`) automatically includes:

| Field | Type | Description |
|-------|------|-------------|
| `schema_version` | string | Envelope schema version; current value lives at `pdb2reaction.core.utils.RESULT_JSON_SCHEMA_VERSION` — pin against that constant rather than the literal in this doc. A version bump signals a structural change. |
| `command` | string | Subcommand name (e.g. `"opt"`) |
| `pdb2reaction_version` | string | Package version |
| `status` | string | Value depends on the subcommand (see each section below): e.g. `converged` / `not_converged` (opt, tsopt), `completed` (irc, freq), `ok` / `partial` (bond-summary), `success` / `partial` (all, path-search), and `error` on failure |
| `elapsed_seconds` | float | Wall-clock time (seconds) |
| `environment` | object | Hardware info (see below) |

### Error envelope (when `status == "error"`)

| Field | Type | Description |
|-------|------|-------------|
| `error` | string | `str(exc)` of the original exception |
| `error_type` | string | Exception class name |
| `error_class_chain` | list[string] | Full MRO class names so agents can match the hierarchy without parsing text |
| `error_module` | string | Module the exception class was defined in |
| `error_label` | string | High-level CLI stage label |

**`environment`**:

| Field | Type | Example |
|-------|------|---------|
| `device` | string | `"cuda"` or `"cpu"` |
| `gpu_name` | string | `"<gpu model>"` |
| `gpu_vram_gb` | float | `<vram in GB>` |
| `cuda_version` | string | `"<cuda version>"` |
| `cpu` | string | `"<cpu model>"` |
| `n_cpus` | int | `<int>` |
| `ram_gb` | float | `<ram in GB>` |

## Error handling

If a job fails (e.g., crash, OOM, convergence failure leading to `sys.exit`), `result.json` is still written with `"status": "error"` plus an `"error_type"` field describing the failure class. Check the `.out` log file for a full traceback. The authoritative signal is `status == "error"`, not the absence of `result.json`.

For jobs that complete but did not converge, `result.json` is written with `"status": "not_converged"` and the final force/step values, allowing an AI agent to decide whether to retry with more cycles.

## Subcommand schemas

### `opt`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"converged"` or `"not_converged"` |
| `energy_hartree` | float | Final energy (Hartree) |
| `n_opt_cycles` | int | Optimization cycles completed |
| `opt_mode` | string | `"grad"` or `"hess"` |
| `backend` | string | MLIP backend (`"uma"`, `"orb"`, `"mace"`, `"aimnet2"`) |
| `charge` | int | System charge |
| `spin` | int | Spin multiplicity |
| `model` | string | MLIP model identifier |
| `n_atoms` | int | Total atoms |
| `n_freeze_atoms` | int | Frozen atoms |
| `solvent` | string | Implicit solvent or `"none"` |
| `thresh` | string | Convergence threshold preset |
| `max_cycles` | int | Maximum allowed cycles |
| `input_file` | string | Input filename |
| `final_max_force` | float | Last max gradient (Hartree/Bohr) |
| `final_rms_force` | float | Last RMS gradient |
| `final_max_step` | float | Last max displacement (Bohr) |
| `final_rms_step` | float | Last RMS displacement |
| `convergence_thresholds` | object | `{max_force_thresh, rms_force_thresh, max_step_thresh, rms_step_thresh}` (Hartree/Bohr) |
| `files` | object | Output file map |

### `tsopt`

All fields from `opt`, plus:

| Field | Type | Description |
|-------|------|-------------|
| `n_imaginary_modes` | int | Number of imaginary frequencies |
| `imaginary_frequencies_cm` | float[] | Imaginary frequencies (cm⁻¹, negative) |
| `opt_mode` | string | `"rsirfo"` or `"dimer"` |

The `files` object may include `imaginary_mode_files` (list of vib file paths).
Convergence details are available for rsirfo mode; dimer mode also reports `status: "converged"` or `"not_converged"` (from `runner.is_converged`), but provides `n_opt_cycles` only and omits the per-cycle force/step convergence keys that rsirfo reports.

### `freq`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"completed"` |
| `n_modes` | int | Total normal modes |
| `n_imaginary` | int | Imaginary frequency count |
| `frequencies_cm` | float[] | All frequencies (cm⁻¹) |
| `imaginary_frequencies_cm` | float[] | Negative frequencies only |
| `thermochemistry` | object\|null | Thermodynamic data (see below) |
| `backend` | string | MLIP backend |
| `charge` | int | System charge |
| `spin` | int | Spin multiplicity |
| `model` | string | Model identifier |
| `n_atoms` | int | Total atoms |
| `n_freeze_atoms` | int | Frozen atoms |
| `solvent` | string | Implicit solvent or `"none"` |
| `temperature_K` | float | Temperature (K) |
| `pressure_atm` | float | Pressure (atm) |
| `input_file` | string | Input filename |
| `files` | object | `{"frequencies_txt": "frequencies_cm-1.txt"}` |

**`thermochemistry`** (null if thermoanalysis unavailable):

| Field | Type | Unit |
|-------|------|------|
| `electronic_energy_ha` | float | Hartree |
| `zpe_correction_ha` | float | Hartree |
| `thermal_correction_energy_ha` | float | Hartree |
| `thermal_correction_enthalpy_ha` | float | Hartree |
| `thermal_correction_free_energy_ha` | float | Hartree |
| `sum_EE_and_ZPE_ha` | float | Hartree |
| `sum_EE_and_thermal_energy_ha` | float | Hartree |
| `sum_EE_and_thermal_enthalpy_ha` | float | Hartree |
| `sum_EE_and_thermal_free_energy_ha` | float | Hartree |
| `E_thermal_cal_per_mol` | float | cal/mol |
| `Cv_cal_per_mol_K` | float | cal/(mol K) |
| `S_cal_per_mol_K` | float | cal/(mol K) |

### `irc`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"completed"` |
| `n_frames_forward` | int | Forward IRC frames |
| `n_frames_backward` | int | Backward IRC frames |
| `n_frames_total` | int | Total frames |
| `energy_reactant_hartree` | float | Reactant energy |
| `energy_ts_hartree` | float | TS energy |
| `energy_product_hartree` | float | Product energy |
| `forward_converged` | bool \| null | Forward IRC converged? `null` when the integrator did not expose the flag |
| `backward_converged` | bool \| null | Backward IRC converged? `null` when the integrator did not expose the flag |
| `backend` | string | MLIP backend |
| `charge` | int | System charge |
| `spin` | int | Spin multiplicity |
| `model` | string | Model identifier |
| `n_freeze_atoms` | int | Frozen atoms |
| `solvent` | string | Implicit solvent or `"none"` |
| `bond_changes` | object | `{formed: [...], broken: [...]}` of element-prefixed 1-based atom-pair strings (e.g. `"C7-O12"`); key is omitted when the comparison fails or `finished_first.xyz`/`finished_last.xyz` are absent. |
| `step_length` | float | IRC step length (Bohr) |
| `max_cycles` | int | Maximum IRC steps |
| `input_file` | string | Input filename |
| `files` | object | Trajectory files (xyz + pdb) |

### `scan`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"completed"` |
| `charge` | int | System charge |
| `spin` | int | Spin multiplicity |
| `backend` | string | MLIP backend |
| `model` | string | Model identifier |
| `solvent` | string | Implicit solvent or `"none"` |
| `preopt` | bool | Pre-optimization performed? |
| `max_step_size_angstrom` | float | Max bond-length step per increment (Å) |
| `n_stages` | int | Number of scan stages |
| `stages` | object[] | Per-stage data (see below) |
| `files` | object | Output files |

**`stages[]`**:

| Field | Type | Description |
|-------|------|-------------|
| `index` | int | 1-based stage index |
| `n_steps` | int | Steps in this stage |
| `converged` | bool | Constrained optimization converged? |
| `pairs_1based` | list | Atom pairs (1-based) |
| `initial_distances_angstrom` | list | Starting distances |
| `target_distances_angstrom` | list | Target distances |
| `final_energy_hartree` | float | Energy at last step |
| `energies_hartree` | float[] | Per-step energies |
| `bond_changes` | object | `{"changed": bool \| null, "summary": str}` from `has_bond_change` (free-text summary; `null`/`""` when the comparison did not run). |

### `scan2d` / `scan3d`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"completed"` |
| `charge` | int | System charge |
| `spin` | int | Spin multiplicity |
| `backend` | string | MLIP backend |
| `model` | string | Model identifier |
| `solvent` | string | Implicit solvent or `"none"` |
| `max_step_size_angstrom` | float | Max bond-length step per increment (Å, `scan2d` only) |
| `n_grid_points` | int | Total grid points |
| `grid_shape` | int[] | Grid dimensions (only when running fresh; absent under `scan3d --csv`) |
| `pair1`, `pair2` (,`pair3`) | object | `{i, j, low, high}` with optional `label_i`, `label_j`. `scan3d`: present only when running fresh; absent under `--csv` re-plot |
| `min_energy_hartree` | float | Surface minimum energy |
| `files` | object | CSV + plot files |

### `path-opt`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"converged"` / `"not_converged"` / `"completed"` |
| `converged` | bool | Convergence flag |
| `mep_mode` | string | `"dmf"` or `"gsm"` |
| `backend` | string | MLIP backend |
| `charge` | int | System charge |
| `spin` | int | Spin multiplicity |
| `model` | string | Model identifier |
| `solvent` | string | Implicit solvent or `"none"` |
| `image_energies_hartree` | float[] | All image energies |
| `n_images` | int | Image count |
| `hei_index` | int | Highest-energy image index |
| `hei_energy_hartree` | float | HEI energy (Hartree) |
| `barrier_kcal` | float | Forward barrier (kcal/mol) |
| `delta_kcal` | float | Reaction energy (kcal/mol) |
| `files` | object | Trajectory + HEI files |

### `path-search`

`path-search` has no `--out-json` flag; it **always** writes a `summary.json` to the output directory with the shared envelope (`command`, `pdb2reaction_version`, `environment`) plus:

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"success"` / `"partial"` |
| `n_segments` | int | Recursive MEP segment count |
| `segments` | object[] | Per-segment `index`, `tag`, `kind`, `barrier_kcal`, `delta_kcal`, `bond_changes` (list of `{title: [entries]}` dicts; bridge segments emit `""`). |
| `energy_diagrams` | object[] | Per-segment labeled energy profiles (kcal/mol) |
| `mlip_backend` | string | Backend/model identifier |
| `charge` | int | System charge |
| `spin` | int | Spin multiplicity |

See also the extended [`summary.json` section](#summary-json-path-search-all) for the additional fields that `all` layers on top.

### `dft`

> **Note:** `dft` writes `result.json` only on a successful SCF (exit 0). A
> non-converged SCF returns exit code 3 and skips `result.json`; SCF status
> is encoded by the `converged: bool` field plus the exit code, not by a
> `status` field. The generic `not_converged` status does not apply to `dft`;
> an unhandled exception, however, still writes the standard `error` envelope
> (`result.json` + `summary.json`).

| Field | Type | Description |
|-------|------|-------------|
| `converged` | bool | SCF converged? |
| `charge` | int | System charge |
| `spin` | int | Spin multiplicity |
| `n_atoms` | int | Atom count |
| `grid_level` | int | DFT grid level |
| `conv_tol` | float | SCF convergence tolerance |
| `max_cycle` | int | Maximum SCF cycles |
| `input_file` | string | Input filename |
| `energy_hartree` | float | DFT energy |
| `energy_kcal_per_mol` | float | DFT energy (kcal/mol) |
| `xc_functional` | string | XC functional |
| `basis_set` | string | Basis set |
| `engine` | string | Effective engine label (`"gpu4pyscf(rks_lowmem)"`, `"gpu4pyscf"`, or `"pyscf(cpu)"`) |
| `used_gpu` | bool | GPU acceleration used? |
| `used_lowmem` | bool | `gpu4pyscf.dft.rks_lowmem.RKS` actually used? (False on open-shell, CPU, `--no-lowmem`, or pre-`rks_lowmem` GPU4PySCF) |
| `charges` | object | `{mulliken, lowdin, iao}` per-atom arrays |
| `spin_densities` | object | `{mulliken, lowdin, iao}` per-atom arrays |
| `files` | object | `{"result_yaml": "result.yaml", "input_geometry_xyz": "input_geometry.xyz"}` |

### `extract`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"ok"` |
| `n_atoms_raw` | int | Atoms in input PDB |
| `n_atoms_extracted` | int | Atoms after extraction |
| `total_charge` | float | Computed total charge |
| `protein_charge` | float | Protein charge |
| `ligand_total_charge` | float | Ligand charge sum |
| `ion_total_charge` | float | Ion charge sum |
| `ion_charges` | list | `[[name, charge], ...]` |
| `unknown_residue_charges` | object | `{resname: charge}` |
| `n_link_hydrogens` | int | Cap hydrogens added at severed C/N bonds |
| `exclude_backbone` | bool | Whether backbone atoms were excluded |
| `include_h2o` | bool | Whether crystallographic waters were included |
| `ligand_charge_input` | string | User-supplied `--ligand-charge` mapping |
| `center` | string | Center residue |
| `radius` | float | Extraction radius (angstrom) |
| `input_files` | string[] | Input PDB paths |
| `files` | object | Output PDB / cluster filenames |

### `trj2fig`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"ok"` |
| `n_frames` | int | Number of trajectory frames |
| `min_energy_hartree` | float | Minimum energy across frames |
| `max_energy_hartree` | float | Maximum energy across frames |
| `backend` | string | MLIP backend |
| `files` | object | Output plot files |

### `energy-diagram`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"ok"` |
| `n_points` | int | Number of energy data points |
| `files` | object | Output diagram files |

### `bond-summary`

When `--json` is enabled, `bond-summary` prints JSON to **stdout** (no `result.json` file is written; redirect stdout if you need to persist it). This is **unlike** the MLIP-based subcommands above, which all write a `result.json` file into `out_dir`:

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"ok"` (every pair compared without error) or `"partial"` (some pairs failed) |
| `comparisons` | object[] | Per-pair comparison with `structure_a` (string), `structure_b` (string), `bonds_formed` (int count), `bonds_broken` (int count). |

(summary-json-path-search-all)=
## `summary.json` (`path-search` / `all`)

The `all` and `path-search` commands write `summary.json` with a richer structure:

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"success"` / `"partial"` / `"failed"` (`all`); `"success"` / `"partial"` (`path-search`) |
| `n_segments` | int | Segment count |
| `segments` | object[] | Per-segment `index`, `tag`, `kind`, `barrier_kcal`, `delta_kcal`, `bond_changes` (list of `{title: [entries]}` dicts produced by `_bond_changes_block`; bridge segments emit `""`). |
| `energy_diagrams` | object[] | Energy profiles with labels and kcal/mol values |
| `mlip_backend` | string | Model identifier |
| `charge` | int | System charge |
| `spin` | int | Spin multiplicity |
| `environment` | object | Hardware info |

The `all` command additionally includes:

| Field | Type | Description |
|-------|------|-------------|
| `rate_limiting_step` | object | RLS segment index and barrier |
| `overall_reaction_energy_kcal` | float | Overall reaction energy |
| `post_segments` | list | Per-segment TS/IRC/freq/DFT results |
| `key_output_files` | object | Curated output file listing |

## Usage examples

### Python

```python
import json

with open("result_opt/result.json") as f:
    result = json.load(f)

if result["status"] == "converged":
    print(f"Energy: {result['energy_hartree']:.6f} Hartree")
else:
    print(f"Not converged after {result['n_opt_cycles']} cycles")
    print(f"Max force: {result['final_max_force']:.6f}")
```

### jq

```bash
# Check convergence
jq '.status' result.json

# Get barrier from path-opt
jq '.barrier_kcal' result.json

# List imaginary frequencies from tsopt
jq '.imaginary_frequencies_cm' result.json

# Get thermochemistry from freq
jq '.thermochemistry.sum_EE_and_thermal_free_energy_ha' result.json
```

## See Also

- [CLI Conventions](cli-conventions.md) — `--out-json` / `--no-out-json` flag conventions and exit codes
- [YAML Reference](yaml-reference.md) — configuration inputs whose values surface in these schemas
- [all](all.md), [path-search](path-search.md) — subcommands that always emit `summary.json`
- [opt](opt.md), [tsopt](tsopt.md), [freq](freq.md), [irc](irc.md), [scan](scan.md), [path-opt](path-opt.md), [dft](dft.md), [extract](extract.md) — subcommands that emit `result.json` only under `--out-json`
