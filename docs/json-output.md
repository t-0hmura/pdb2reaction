# JSON Output Reference

This page lists the fields (keys) of the `result.json` and `summary.json` files written with `--out-json`, split into the fields every command shares and the fields of each command.

## `--out-json` flag

`opt`, `sp`, `tsopt`, `freq`, `irc`, `scan`, `scan2d`, `scan3d`, `path-opt`,
`dft`, `extract`, `trj2fig`, and `energy-diagram` support
`--out-json / --no-out-json` (default: off). When enabled, `result.json` and `summary.json` are written
beside the normal outputs; the two files have the same content, so read
`result.json`.

```bash
pdb2reaction opt -i r.pdb -q -1 --out-json --out-dir result_opt
cat result_opt/result.json | python -m json.tool
```

Open `result_opt/result.json` and read `execution_status` and `scientific_status` first; `opt`, `tsopt`, and `path-opt` also record `optimization_status`.
The `summary.json` that `all` and `path-search` write without an
`--out-json` flag is [a different file with its own structure](#summary-json-path-search-all).

## Common envelope

Every result file carries the fields below; fields marked optional appear only
when the command has the corresponding data:

| Field | Type | Description |
|-------|------|-------------|
| `schema_version` | string | Schema version of the file; a new version signals a structural change. |
| `command` | string | Single commands record the subcommand name (e.g. `"opt"`); the `all` / `path-search` summaries record the full command line |
| `pdb2reaction_version` | string | Package version |
| `execution_status` | string | Execution completion: `completed` / `failed`. |
| `scientific_status` | string | Result usability: `success` / `partial` / `failed`. |
| `run_id` | string | Optional UUID of the current invocation; written when the MCP server starts the command and in every `all` run, including its stages. |
| `elapsed_seconds` | float | Optional wall-clock time; omitted by commands that do not record timing |
| `environment` | object | Hardware info (see below) |
| `mlip_backend` | string \| null | Optional backend identifier; null when a plot-only command did not evaluate a calculator |
| `mlip_model` | string \| null | Optional exact model/checkpoint identifier, kept separate from the backend |
| `mlip_model_label` | string \| null | Optional publication-facing model label derived from the exact identifier |
| `mlip_task` | string \| null | Optional backend task used by a multi-domain model; `mlip_model` remains the exact identifier |
| `mlip_precision` | string \| null | Effective public precision token (`fp32` or `fp64`); null for a custom calculator whose dtype is controlled by user code |

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

### Execution and requested-stage completion

Every result reports `execution_status` and `scientific_status`; multi-stage and scan results also list the outcome of each stage in the fields below.

| Field | Type | Description |
|-------|------|-------------|
| `execution_status` | string | `completed` or `failed`. In `all`, an endpoint optimization after IRC that does not converge leaves it `completed`; one that stops on an error makes it `failed`. |
| `scientific_status` | string | `success` when every requested stage converged; otherwise `partial` or `failed`. In `all`, a TS with n_imag ≥ 2 gives `partial`, and a TS with n_imag = 0 stops the run before IRC, so `success` means n_imag = 1. Standalone `tsopt` does not look at n_imag; read `n_imaginary_modes`. |
| `scientific_status_reasons` | string[] | Reasons for unusable or missing stages; omitted on clean success. |
| `expected_item_ids` / `observed_item_ids` | string[] | Expected and observed stage identifiers, used to detect missing work. |
| `stage_outcomes` | object[] | One entry per stage with `stage`, `item_id`, `required`, `executed`, `converged`, `usable`, `reason`, and `artifacts`. |
| `point_outcomes` | object[] | Scan points with `point_id`, `executed`, `converged`, `energy_valid`, `artifact_written`, `seed_eligible`, and `reason`. |

### Error envelope (when `execution_status == "failed"`)

| Field | Type | Description |
|-------|------|-------------|
| `error` | string | `str(exc)` of the original exception |
| `error_type` | string | Exception class name |
| `error_class_chain` | list[string] | Names of the exception class and all its parent classes, so agents can match the hierarchy without parsing text |
| `error_module` | string | Module the exception class was defined in |
| `error_label` | string | High-level CLI stage label |

## Error handling

When a run stops on an exception after the output directory is set up,
`result.json` and `summary.json` are written even without `--out-json`, with
`"execution_status": "failed"` and an `"error_type"`; for a failure before that
point, see [Notes](#notes).

When an optimization ends without converging, `result.json` records
`"optimization_status": "not_converged"`. Optimizer results include the final
force/step or cycle fields that apply; DFT and Dimer omit the fields they do not
have. For what to change before a retry, see
{ref}`Troubleshooting › Calculation / convergence <calculation-convergence-problems>`.

An optimizer may also report `"optimization_status": "stalled"`: the energy stopped decreasing over the configured window (an energy plateau) while the force/step convergence criteria remained unmet. A stall is a kind of non-convergence, never `converged`; `stop_reason` records why the run stopped.

## Subcommand schemas

### `sp`

| Field | Type | Description |
|-------|------|-------------|
| `stage` | string | `"sp"` |
| `input` | string | Input path used for the calculation |
| `backend` / `model` | string / string \| null | MLIP backend and model; the same values appear in the common `mlip_*` fields |
| `custom_calculator` | string \| null | `filename:factory` for `--calc-file`, otherwise null |
| `charge` / `spin` | int / int | Total charge and multiplicity |
| `n_atoms` | int | Atom count |
| `energy_au` | float | Single-point energy (Hartree) |
| `forces_path` | string | Path to `forces.npy` |
| `hessian_path` | string \| null | Path to `hessian.npy`, or null without `--hess` |
| `elapsed` | string | Human-readable elapsed-time text |

### `opt`

| Field | Type | Description |
|-------|------|-------------|
| `optimization_status` | string | `"converged"`, `"not_converged"`, or `"stalled"` (energy plateau; see above) |
| `stop_reason` | string | Present only when the optimizer stopped early without converging (`stalled` or `not_converged`); records why, e.g. the energy-plateau range/window and the failed criteria |
| `energy_hartree` | float | Final energy (Hartree) |
| `n_opt_cycles` | int | Optimization cycles completed |
| `opt_mode` | string | `"grad"`, `"hess"`, `"lbfgs"`, or `"rfo"` |
| `backend` | string | Calculator backend (`"uma"`, `"orb"`, `"mace"`, `"aimnet2"`, `"dft"`, or `"custom"` with `--calc-file`); for `dft`, `model` is `FUNCTIONAL/BASIS`, the engine is recorded separately, and MLIP precision is null |
| `charge` | int | System charge |
| `spin` | int | Spin multiplicity |
| `model` | string | MLIP model identifier, or `FUNCTIONAL/BASIS` for `dft` |
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
| `convergence_thresholds` | object | `{max_force_thresh, rms_force_thresh, max_step_thresh, rms_step_thresh}`. Force/step units: Hartree/Bohr and Bohr (Cartesian), Hartree/rad and rad (angular). |
| `files` | object | Output file map |
| `rigid_projection` | object | Optional; present when `--flatten` runs. See [projection provenance](#rigid-projection-provenance). |

### `tsopt`

All fields from `opt`, plus:

| Field | Type | Description |
|-------|------|-------------|
| `optimization_status` | string | Numerical optimizer outcome: `"converged"`, `"not_converged"`, or `"stalled"`; independent of saddle order |
| `saddle_validation` | string | `"first_order"`, `"higher_order"`, `"no_imaginary"`, or `"unavailable"` from terminal exact PHVA (partial Hessian vibrational analysis) |
| `hessian_status` | string | `"completed"`, `"failed"`, `"skipped"`, or `"unavailable"`; `hessian_error` gives the failure reason |
| `reaction_mode_index` | int\|null | 0-based index, in the PHVA frequency list, of the imaginary mode that `all` follows in IRC: the mode tracked by the optimizer when it is imaginary, otherwise the lowest imaginary mode. `reaction_mode_source` records which (`"mep-reference-overlap"` or `"lowest-imaginary"`); `null` when there is no imaginary mode |
| `n_imaginary_modes` | int\|null | Number of imaginary frequencies; `null` if PHVA was not run |
| `n_negative_modes` | int\|null | Number of negative frequencies of any size, including those within the cutoff; a diagnostic next to `n_imaginary_modes`, `null` if PHVA was not run |
| `imaginary_frequencies_cm` | float[]\|null | Imaginary frequencies (cm⁻¹, negative); `null` when PHVA was not run |
| `frequency_zero_cutoff_cm` | float | Cutoff in cm⁻¹ (default `5.0`); with the default, only ν < −5.00 cm⁻¹ counts as imaginary |
| `imaginary_mode_criterion` | string | Name of the counting rule; the value is `"frequency_cutoff_cm"` |
| `imaginary_frequency_threshold_cm` | float | The same cutoff with a negative sign (default `-5.0`) |
| `opt_mode` | string | `"rsprfo"` (default), `"rsirfo"`, `"trim"`, or `"dimer"` |
| `opt_mode_requested` | string | Requested CLI preset (`grad`, `hess`, or an explicit algorithm) |
| `optimizer` | string | Effective optimizer algorithm used by the run |
| `reference_mode_file` | string\|null | Advanced path-mode file supplied through `--ref-mode`; normally generated and passed by `all` |
| `safeguards` | object | Hessian-TS diagnostics, including exact saddle checks, final target-mode identity/overlap, stop reason, and any explicitly enabled mode-loss/recovery activity. These recovery paths are inactive by default. |
| `rigid_projection` | object | Rigid-mode and exact-Hessian provenance; see [projection provenance](#rigid-projection-provenance) |

The `files` object may include `imaginary_mode_files` (list of vib file paths) and `hessian_npy` (absolute path of the `--dump-hess` file).
A successful TS optimization gives one imaginary mode along the reaction
coordinate: `saddle_validation: "first_order"` with `n_imaginary_modes: 1`.
The final PHVA runs after the optimizer converges or stops on an energy plateau
(`stalled`); otherwise it is recorded as skipped. `optimization_status` and
`saddle_validation` are independent, so a converged run can end with
`saddle_validation: "higher_order"`; it is not a first-order TS. For when `all`
goes on to IRC, see [tsopt › Reading the TS result](tsopt.md#reading-the-ts-result).
Dimer results have the same TS
fields but no per-cycle force/step details and no `safeguards` object.

### `freq`

| Field | Type | Description |
|-------|------|-------------|
| `n_modes` | int | Total normal modes |
| `n_imaginary` | int | Imaginary frequency count (n_imag): frequencies below −`freq.zero_cutoff_cm` |
| `n_negative_modes` | int | Number of negative frequencies of any size, including those within the cutoff |
| `frequencies_cm` | float[] | All frequencies (cm⁻¹) |
| `imaginary_frequencies_cm` | float[] | The frequencies counted in `n_imaginary` |
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
| `files` | object | `{"frequencies_txt": "frequencies_cm-1.txt"}`; includes `hessian_npy` (absolute path) when `--dump-hess` wrote a file |
| `rigid_projection` | object | Rigid-mode and Hessian provenance; also written to `thermoanalysis.yaml` with `--dump` |

**`thermochemistry`**:

| Field | Type | Unit |
|-------|------|------|
| `point_group` | string | Automatically detected molecular point group |
| `point_group_source` | string | `"auto"` or conservative `"auto-fallback"` |
| `symmetry_number` | int | External rotational symmetry number |
| `symmetry_number_source` | string | `"auto"`, `"auto-fallback"`, `"config"`, or `"override"` |
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

IRC reports `execution_status` and `scientific_status` and keeps the stop reason and trajectory of each direction; `all` reports the endpoint optimizations that follow under `endpoint_opt`. The stitched path runs from `finished_first` (the forward end) through the TS to `finished_last` (the backward end); which end is R and which is P is not decided by the order.

| Field | Type | Description |
|-------|------|-------------|
| `n_frames_forward` | int | Forward IRC frames |
| `forward_short_branch` / `backward_short_branch` | bool | Branch produced at most three frames without reaching the cycle cap; diagnostic only |
| `n_frames_backward` | int | Backward IRC frames |
| `n_frames_total` | int | Total frames |
| `energy_first_hartree` | float | Energy of `finished_first` |
| `energy_ts_hartree` | float | TS energy |
| `energy_last_hartree` | float | Energy of `finished_last` |
| `endpoint_energy_orientation` | string | `"finished_first_to_finished_last"` |
| `forward_requested` / `backward_requested` | bool | Whether each direction was requested |
| `forward_integration_converged` / `backward_integration_converged` | bool \| null | Whether the direction stopped because the RMS-gradient stationarity criterion fired; diagnostic only, and always `false` under `--never-stop`, which bypasses that criterion. Combine it with `*_downhill_departure_valid` to check that the branch both left the TS downhill and met that criterion. |
| `forward_downhill_departure_valid` / `backward_downhill_departure_valid` | bool \| null | Whether the branch established a downhill departure from the TS |
| `forward_integration_stop_reason` / `backward_integration_stop_reason` | string \| null | Non-empty only for a numerical propagation failure |
| `forward_energy_increased` | bool \| null | Final forward step exceeded `irc.energy_increase_thresh` (default `0` Hartree: any rise) |
| `backward_energy_increased` | bool \| null | Final backward step exceeded `irc.energy_increase_thresh` (default `0` Hartree: any rise) |
| `backend` | string | MLIP backend |
| `charge` | int | System charge |
| `spin` | int | Spin multiplicity |
| `model` | string | Model identifier |
| `never_stop` | bool | Whether opt-in physical endpoint-stop bypass mode was enabled |
| `never_stop_energy_bypasses` | int | Number of energy-rise or one-step energy-change stop events actually bypassed |
| `n_freeze_atoms` | int | Frozen atoms |
| `solvent` | string | Implicit solvent or `"none"` |
| `bond_changes` | object | Directed first→last `{formed: [...], broken: [...]}` of element-prefixed 1-based atom-pair strings (e.g. `"C7-O12"`); key is omitted when the comparison fails or `finished_first.xyz`/`finished_last.xyz` are absent. |
| `bond_changes_direction` | string | `"finished_first_to_finished_last"` when `bond_changes` is present |
| `step_length` | float | IRC step length (Bohr) |
| `max_cycles` | int | Maximum IRC steps |
| `input_file` | string | Input filename |
| `files` | object | Trajectory files (XYZ, plus PDB/CIF versions when available) |
| `rigid_projection` | object | Rigid-mode and initial-Hessian provenance; see [projection provenance](#rigid-projection-provenance) |

### `scan`

| Field | Type | Description |
|-------|------|-------------|
| `scan_opt_mode` | string | Optimizer preset used for the constrained relaxations |
| `scan_optimizer` | string | Effective optimizer identity (`lbfgs` or `rfo`) |
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
| `bond_changes` | object | `{"changed": bool \| null, "summary": str}` (free-text summary; `null`/`""` when the comparison did not run). |

### `scan2d` / `scan3d`

| Field | Type | Description |
|-------|------|-------------|
| `charge` | int \| null | System charge; null for plot-only `scan3d --csv` |
| `spin` | int \| null | Spin multiplicity; null for plot-only `scan3d --csv` |
| `backend` | string \| null | MLIP backend; null for plot-only `scan3d --csv` |
| `model` | string \| null | Model identifier; null for plot-only `scan3d --csv` |
| `solvent` | string \| null | Implicit solvent or `"none"`; null for plot-only `scan3d --csv` because imported energies have no calculator provenance |
| `max_step_size_angstrom` | float | Max bond-length step per increment (Å, `scan2d` only) |
| `n_grid_points` | int | Grid rows excluding `is_preopt=true` |
| `execution_status` | string | Execution-level completion state |
| `n_points_attempted` | int | Fresh-run grid points attempted, excluding preoptimization |
| `n_points_usable` | int | Fresh-run points eligible for scientific reuse |
| `point_outcomes` | object[] | Per-point convergence, energy, artifact, and eligibility data |
| `grid_points` | object[] | Explicit grid-index, distances, energy, convergence, and `geometry_file` mapping |
| `current_output_paths` | string[] | CSV/HTML/PNG files and grid geometries written by this run; files left from an earlier run are not listed |
| `grid_shape` | int[] | Grid dimensions (only when running fresh; absent under `scan3d --csv`) |
| `pair1`, `pair2` (,`pair3`) | object | `{i, j, low, high}` with optional `label_i`, `label_j`. `scan3d`: present only when running fresh; absent under `--csv` re-plot |
| `min_energy_hartree` | float | Surface minimum energy |
| `files` | object | CSV + plot files |

The outcome-count fields are emitted by fresh scans. Plot-only `scan3d --csv`
does not report an attempted count and may omit usable count when the imported
CSV lacks complete provenance.

### `path-opt`

| Field | Type | Description |
|-------|------|-------------|
| `optimization_status` | string | `"converged"` / `"not_converged"` / `"completed"` |
| `converged` | bool \| null | Convergence flag: `true` / `false` from the engine's own convergence signal, `null` when it exposed none (`optimization_status` is then `"completed"`, never a success claim) |
| `mep_mode` | string | `"dmf"` or `"gsm"` |
| `backend` | string | MLIP backend |
| `charge` | int | System charge |
| `spin` | int | Spin multiplicity |
| `model` | string | Model identifier |
| `solvent` | string | Implicit solvent or `"none"` |
| `preopt` | bool | Whether endpoint pre-optimization was enabled |
| `reactant_energy_hartree` | float | First-image energy (Hartree) |
| `product_energy_hartree` | float | Last-image energy (Hartree) |
| `image_energies_hartree` | float[] | All image energies |
| `n_images` | int | Image count |
| `hei_index` | int | Highest-energy image index |
| `hei_energy_hartree` | float | HEI energy (Hartree) |
| `barrier_kcal` | float | Forward barrier (kcal/mol) |
| `delta_kcal` | float | Reaction energy (kcal/mol) |
| `files` | object | Trajectory + HEI files |

### `path-search`

`path-search` has no `--out-json` flag. It writes `summary.json`; its fields
are listed in [`summary.json` (`path-search` / `all`)](#summary-json-path-search-all).

### `dft`

> **Note:** With `--out-json`, `dft` writes `result.json` and `summary.json` for
> both converged and non-converged SCF attempts, recording
> `scientific_status: "failed"` and `converged: false` on non-convergence.
> A non-converged SCF exits with code 1.

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
| `used_lowmem` | bool | Low-memory GPU4PySCF solver actually used? (False on open-shell, CPU, or `--no-dft-low-memory`) |
| `lowmem_requested` | bool | Whether low-memory mode was requested |
| `dft_settings` / `dft_resources` | object | Canonical scientific settings and effective host resources |
| `effective_ecp` | string/object \| null | Effective ECP passed to PySCF |
| `solvent` / `solvent_model` | string | Effective native implicit-solvent settings |
| `charges` | object | `{mulliken, lowdin, iao}` per-atom arrays |
| `spin_densities` | object | `{mulliken, lowdin, iao}` per-atom arrays |
| `files` | object | `{"result_yaml": "result.yaml", "input_geometry_xyz": "input_geometry.xyz"}` |

### `extract`

| Field | Type | Description |
|-------|------|-------------|
| `n_atoms_raw` | int | Atoms in selected residues before backbone/truncation filtering (not the whole input structure) |
| `n_atoms_extracted` | int | Selected atoms kept after truncation, before cap-H addition |
| `total_charge` | float | Computed total charge |
| `protein_charge` | float | Protein charge |
| `ligand_total_charge` | float | Ligand charge sum |
| `ion_total_charge` | float | Ion charge sum |
| `ion_charges` | list | `[[name, charge], ...]` |
| `unknown_residue_charges` | object | `{resname: charge}` |
| `n_link_hydrogens` | int | Cap hydrogens added at carbon-parent truncation bonds; the model has `n_atoms_extracted` + `n_link_hydrogens` atoms |
| `exclude_backbone` | bool | Whether backbone atoms were excluded |
| `include_h2o` | bool | Whether crystallographic waters were included |
| `ligand_charge_input` | string \| null | User-supplied mapping, or null when omitted |
| `center` | string | Center residue |
| `radius` | float | Extraction radius (angstrom) |
| `input_files` | string[] | Original input PDB/mmCIF paths |
| `files` | object | Output PDB / cluster filenames |

### `trj2fig`

| Field | Type | Description |
|-------|------|-------------|
| `n_frames` | int | Number of trajectory frames |
| `min_energy_hartree` | float | Minimum energy across frames |
| `max_energy_hartree` | float | Maximum energy across frames |
| `energy_source` | string | `"trajectory_comment"` or `"mlip_recomputed"` |
| `mlip_backend` / `mlip_model` / `mlip_model_label` / `mlip_task` / `mlip_precision` | string \| null | Effective recomputation provenance; all are null in trajectory-comment mode |
| `energy_provenance` | string[] | Per-frame energy source provenance |
| `energy_unit` | string | Stored energy unit (`hartree`) |
| `backend` | string or null | MLIP backend only when frame energies were recomputed; null in comment-energy mode |
| `charge` / `multiplicity` | int or null | Charge and multiplicity used when energies were recomputed, otherwise null |
| `solvent` / `solvent_model` | string or null | Recomputed-calculator solvent settings, otherwise null |
| `output_files` | string[] | Canonical ordered paths for every output; preserves files with the same basename in different directories |
| `files` | object | Basename-to-path map; when two outputs share a basename only one is kept, so prefer `output_files` |

### `energy-diagram`

| Field | Type | Description |
|-------|------|-------------|
| `n_points` | int | Number of energy data points |
| `files` | object | Output diagram files |

### `bond-summary`

With `--json`, `bond-summary` prints JSON to **stdout** and, unlike the subcommands above, writes no `result.json` file:

| Field | Type | Description |
|-------|------|-------------|
| `execution_status` / `scientific_status` | string / string | `completed` / `success` when every pair was compared. If any pair could not be compared, `execution_status` is `failed` and `scientific_status` is `partial` or `failed`. |
| `comparisons` | object[] | Per-pair comparison with `structure_a` (string), `structure_b` (string), `bonds_formed` (int count), `bonds_broken` (int count). |

### Rigid projection provenance

`freq`, `irc`, and `tsopt` results include a `rigid_projection` object; `opt` includes it when `--flatten` runs. `freq --dump` also writes the same object to `thermoanalysis.yaml`.

| Field | Type | Description |
|-------|------|-------------|
| `treatment` | string | Fixed rigid-mode treatment: `"constrained"` |
| `algorithm` | string | Name of the projection method |
| `effective_rank` | int | Number of rigid directions removed from the Hessian of the movable atoms |
| `full_rigid_rank` | int | Rank of the rigid motions of the whole system before the frozen atoms are taken into account |
| `frozen_constraint_rank` | int | Rank removed because the frozen atoms must stay in place |
| `svd_rtol` | float | Relative SVD tolerance used for the rank decision |
| `active_atom_count` / `frozen_atom_count` | int | Movable and frozen atom counts |
| `active_atoms` / `frozen_atoms` | int[] | 0-based indices of the movable and frozen atoms |
| `hessian_space` | string | `"full"` or `"active"` input Hessian space |
| `hessian_source` / `source` | string | Hessian provenance. `freq`/`irc` use `hessian_source`: `"file"` (`--read-hess`), `"cache"` (earlier stage in the same run), or `"fresh"`; `opt`/`tsopt` use `source`. |
| `hessian_shape` / `raw_hessian_shape` | int[2] | Input Hessian shape. `freq`/`irc` use `hessian_shape`; `opt`/`tsopt` use `raw_hessian_shape`. |
| `near_zero_mode_count` / `near_zero_frequencies_cm` | int / float[] | Number and values of the modes within ±`frequency_zero_cutoff_cm` (5.00 cm⁻¹ by default); these modes are also in the full frequency list |

`constrained` removes only the rigid motions of the whole system that keep the frozen atoms in place; see [freq](freq.md#rigid-modes-with-frozen-boundaries).

(summary-json-path-search-all)=
## `summary.json` (`path-search` / `all`)

The `all` and `path-search` commands write `summary.json` with a richer structure:

| Field | Type | Description |
|-------|------|-------------|
| `execution_status` / `scientific_status` | string / string | Execution completeness and completion of requested numerical/calculation stages. |
| `scientific_status_reasons` | string[] | Reasons for missing or unusable requested results; omitted on success. |
| `pipeline_stop` | object \| absent | Present only on an early stop. `stage` is `post` (`reason` `no_segments` / `no_reactive_segment`), `before_irc` (a TSOPT reason, plus `segment` and `tsopt_result`), or `endpoint_opt` (`endpoint_execution_failed` and endpoint-specific `failures`). Rendered in `summary.log` as `Pipeline stop`. |
| `expected_item_ids` / `observed_item_ids` | string[] | Expected and observed stage identifiers. |
| `config` | object | Effective settings. `mep_mode` identifies GSM/DMF; `ts_opt_mode` and `endpoint_opt_mode` identify the configured post-processing presets. Generic `opt_mode*` keys record the effective CLI values. `path_opt_mode` is the single-structure optimizer used for endpoint preoptimization (see `preopt`), not the MEP path algorithm. |
| `scan` | object \| absent | Preliminary scan status, stage outcomes and diagnostics in scan-seeded `all` runs. |
| `n_segments` | int | Segment count |
| `search_max_depth` | int | Effective recursion cap; `0` means subdivision was disabled |
| `path_optimizers` | string[] | Single-structure optimizers actually used during path preparation/refinement (`lbfgs`, `rfo`); includes scan and alignment work in `all`. Also present in `path-opt` `result.json` |
| `preopt_requested` / `preopt_converged` | bool / bool \| null | Whether endpoint preoptimization ran, and whether every endpoint converged; `null` when any endpoint reported no readable signal. In `all`, `preopt_converged` counts toward `scientific_status`; with `--tsopt`, it no longer counts once the TS and both endpoint optimizations have converged for every reactive segment |
| `segments` | object[] | Per-segment `index` (1-based; `all` writes segment 1 to `segments/seg_01/`), `tag` (a label; its number can differ from `index`), `kind` (`"seg"` for a reactive segment, `"kink"` for a [kink](path-search.md#how-it-works) with no covalent bond change, `"bridge"` for a [bridge segment](path-search.md#how-it-works), `"tsopt"` in TS-only mode), `converged` (whether every optimization that built the segment converged; `null` when no convergence signal was readable), `barrier_kcal`, `delta_kcal`, `bond_changes` (list of `{title: [entries]}` dicts; bridge segments emit `""`). `barrier_kcal` is the barrier on the MEP before TS optimization (`null` for a kink); in TS-only mode it is TS − R, where R is the higher-energy IRC endpoint (a name for the direction, not the chemical direction of the reaction; see [all › Notes](all.md#notes)). |
| `energy_diagrams` | object[] | Energy profiles with labels and kcal/mol values |
| `mlip_backend` | string | Backend identifier |
| `mlip_model` | string \| null | Model identifier, recorded separately from the backend |
| `mlip_model_label` | string \| null | Publication-facing model label |
| `mlip_task` | string \| null | Backend task for a multi-domain model |
| `mlip_precision` | string \| null | Effective `fp32` / `fp64` token; null for custom calculators |
| `charge` | int | System charge |
| `spin` | int | Spin multiplicity |
| `environment` | object | Hardware info |
| `references` | object[] | Methods actually used by the run, as `{method, citation, doi}` records. The same reference set is grouped at the end of `summary.log` and final stdout immediately before elapsed time. |

The `all` command additionally includes:

| Field | Type | Description |
|-------|------|-------------|
| `rate_limiting_step` | object | Highest local barrier among the reactive segments, taken at the highest-level method available for every segment (`DFT//MLIP_Gibbs` > `DFT` > `MLIP_Gibbs` > `MLIP` > `MEP`), with explicit `method` and raw `mep_barrier_kcal`. It is not a microkinetic rate-limiting-step assignment. |
| `overall_reaction_energy_kcal` | float | Overall reaction energy |
| `overall_reaction_energy_method` | string | Method of the overall reaction energy (`MEP`, `MLIP`, `MLIP_Gibbs`, `DFT`, or `DFT//MLIP_Gibbs`) |
| `post_segments` | list | Per-segment TS/IRC/freq/DFT results |
| `post_segments[].tsopt.n_imaginary_modes` / `.imaginary_frequencies_cm` | int / float[] | n_imag of the optimized TS and its imaginary frequencies (cm⁻¹, negative) |
| `post_segments[].mlip` / `.gibbs_mlip` / `.dft` / `.gibbs_dft_mlip` | object | R, TS, and P energies at one level, with `barrier_kcal`, `delta_kcal`, and `energies_kcal`: MLIP electronic energies (`--tsopt`), MLIP Gibbs energies (`--thermo`), DFT energies (`--dft`), and DFT//MLIP Gibbs energies (`--thermo` with `--dft`) |
| `post_segments[].tsopt.energy_valid` / `.structure_valid` | bool | `energy_valid`: the final TS energy is a finite number. `structure_valid`: the final TS structure file exists and its coordinates are finite |
| `post_segments[].tsopt.n_opt_cycles` / `.max_cycles` | int / int\|null | TS optimization cycles executed and configured limit. These are reported for both converged and normally non-converged runs. |
| `post_segments[].irc` / `.endpoint_assignment` / `.endpoint_opt` | object | IRC stop diagnostics, endpoint orientation, and endpoint-OPT convergence, respectively. `endpoint_opt.reactant` and `.product` report `optimization_status`, `n_opt_cycles`, `max_cycles`, and any `stop_reason`. How the IRC stopped and whether the bond changes match do not enter `scientific_status`; whether the optimized endpoints are the intended R and P is for you to check. |
| `post_segments[].thermo_symmetry` | object | Child-reported point-group and rotational-symmetry provenance by state. Only R/TS/P states with valid symmetry-number provenance are included; missing states are omitted, and the field is absent only when no state has valid provenance. |
| `current_output_paths` | string[] | Sorted paths relative to `--out-dir`, limited to artifacts claimed by the current invocation. |
| `key_output_files` | object | Current-run output index: root filename → description; each `seg_NN` entry is `{description, files}` with paths relative to that segment directory. |

## Usage examples

### Python

The script reads `result_opt/result.json`, which `opt --out-json` writes.

```python
import json

with open("result_opt/result.json") as f:
    result = json.load(f)

status = result.get("optimization_status")
if result["execution_status"] == "failed":
    raise RuntimeError(f"{result['error_type']}: {result['error']}")
elif status == "converged":
    print(f"Energy: {result['energy_hartree']:.6f} Hartree")
elif status in {"not_converged", "stalled"}:
    print(f"Not converged after {result['n_opt_cycles']} cycles")
    print(f"Max force: {result['final_max_force']:.6f}")
else:
    print(f"Status: {status}")
```

### jq

```bash
# Check convergence
jq '{execution_status, scientific_status}' result.json

# Get barrier from path-opt
jq '.barrier_kcal' result.json

# List imaginary frequencies from tsopt
jq '.imaginary_frequencies_cm' result.json

# Get thermochemistry from freq
jq '.thermochemistry.sum_EE_and_thermal_free_energy_ha' result.json

# Get each segment's barrier from all (after --tsopt)
jq '.post_segments[] | {index, barrier_kcal: .mlip.barrier_kcal}' result_all/summary.json
```

## Notes

- A run that fails while the CLI options or the input are being checked (before the output directory is set up) can stop without writing any JSON. Treat a nonzero exit code as a failure, and read stderr or the job log for the message.
- `all` and `path-search` write `summary.json` only once the run reaches its summary step, so an early input error leaves no file.
- The per-command sections also have their own `backend` / `model` fields. Scripts that read several commands should use `mlip_backend` / `mlip_model` / `mlip_precision`, which mean the same thing in every command.
- [How the IRC stopped](irc.md#judging-the-irc) does not enter `scientific_status`: a standalone `irc` reports `success` when the integration finishes, also at `--max-cycles`, and `all` judges the TS and endpoint optimizations.

## See Also

- {ref}`Exit codes <exit-codes>` — what each exit code means
- [Output layout](output-layout.md) — where `result.json` and `summary.json` are written
- [Troubleshooting](troubleshooting.md) — what to change after a failed or non-converged run
- [YAML Reference](yaml-reference.md) — configuration inputs whose values surface in these schemas
- [all](all.md), [path-search](path-search.md) — subcommands that write `summary.json` without `--out-json`
- [opt](opt.md), [sp](sp.md), [tsopt](tsopt.md), [freq](freq.md), [irc](irc.md), [scan](scan.md), [scan2d](scan2d.md), [scan3d](scan3d.md), [path-opt](path-opt.md), [dft](dft.md), [extract](extract.md), [trj2fig](trj2fig.md), [energy-diagram](energy-diagram.md) — subcommands that emit `result.json` only under `--out-json`
