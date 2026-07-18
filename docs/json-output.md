# JSON Output Reference

pdb2reaction emits machine-readable JSON for programmatic use by scripts and agents.

## `--out-json` flag

Every MLIP-based subcommand supports `--out-json / --no-out-json` (default: off).
When enabled, a `result.json` file is written to the output directory alongside the normal outputs.

```bash
pdb2reaction opt -i r.pdb -q -1 --out-json --out-dir result_opt
cat result_opt/result.json | python -m json.tool
```

The `all` and `path-search` commands write `summary.json` without an
`--out-json` flag once execution reaches their summary writer. Early
CLI/input validation can fail before the file exists.

### `summary.json` mirror

`write_result_json` serializes each per-stage payload once, atomically publishes
the `summary.json` mirror first, and publishes authoritative `result.json`
last. A successful return means both files contain identical bytes. If an I/O
failure interrupts publication, the writer raises instead of hiding mirror
divergence; MCP consumers additionally require the current run identity.

## Common envelope

The shared writer supplies the fields below; rows explicitly marked optional
are present only when the producer supplies the corresponding data:

| Field | Type | Description |
|-------|------|-------------|
| `schema_version` | string | Envelope schema version; current value lives at `pdb2reaction.core.utils.RESULT_JSON_SCHEMA_VERSION` — pin against that constant rather than the literal in this doc. A version bump signals a structural change. |
| `command` | string | Leaf envelopes use the subcommand name (e.g. `"opt"`); aggregate `all` / `path-search` summaries record the full invocation string |
| `pdb2reaction_version` | string | Package version |
| `status` | string | Value depends on the subcommand (see each section below): e.g. `converged` / `not_converged` / `stalled` (opt, tsopt), `completed` (irc, freq), `ok` / `partial` (bond-summary), `success` / `partial` / `failed` (aggregate workflows), and `error` on failure |
| `run_id` | string | Optional current MCP invocation UUID. It is injected only when the product-specific run environment is active; a conflicting producer value is rejected. |
| `elapsed_seconds` | float | Optional wall-clock time; omitted by producers that do not pass timing to the shared writer |
| `environment` | object | Hardware info (see below) |
| `mlip_backend` | string \| null | Optional backend identifier; null when a plot-only command did not evaluate a calculator |
| `mlip_model` | string \| null | Optional exact model/checkpoint identifier, kept separate from the backend |
| `mlip_precision` | string \| null | Effective public precision token (`fp32` or `fp64`); null for a custom calculator whose dtype is controlled by user code |

Leaf schemas also retain their local `backend` / `model` fields; consumers
should prefer `mlip_backend` / `mlip_model` / `mlip_precision` for a uniform
cross-command provenance contract.

### Rigid projection provenance

`freq`, `irc`, and `tsopt` results include a `rigid_projection` object; `opt` includes it when `--flatten` runs. `freq --dump` also writes the same object to `thermoanalysis.yaml`.

| Field | Type | Description |
|-------|------|-------------|
| `treatment` | string | Selected `--tr-projection` mode: `"constrained"` or `"legacy-active"` |
| `algorithm` | string | Projection-kernel identifier |
| `effective_rank` | int | Number of rigid directions removed from the active Hessian |
| `full_rigid_rank` | int | Rank of the full-system rigid basis before frozen-anchor constraints |
| `frozen_constraint_rank` | int | Rank removed by the frozen-anchor constraints |
| `active_atom_count` / `frozen_atom_count` | int | Active and frozen atom counts |
| `active_atoms` / `frozen_atoms` | int[] | 0-based atom indices used by the projection kernel |
| `hessian_space` | string | `"full"` or `"active"` input Hessian space |
| `hessian_source` / `source` | string | Hessian provenance. `freq`/`irc` use `hessian_source`; `opt`/`tsopt` use `source`. |
| `hessian_shape` / `raw_hessian_shape` | int[2] | Input Hessian shape. `freq`/`irc` use `hessian_shape`; `opt`/`tsopt` use `raw_hessian_shape`. |

`constrained` removes only full-system rigid motions that leave frozen anchors fixed. `legacy-active` is an isolated-active comparison treatment and is not a bitwise replay guarantee for near-linear or degenerate structures. See [Frozen Atoms](freeze-atoms.md#rigid-modes-with-frozen-boundaries).

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

Caught runtime exceptions write a best-effort `result.json` with
`"status": "error"` and an `"error_type"`. Usage/validation exits and failures
before the output directory is resolved may not create JSON, so a nonzero exit
code or missing expected JSON is also a failure signal. Check stderr/job logs
for the authoritative diagnostic.

For jobs that complete but did not converge, `result.json` is written with `"status": "not_converged"` and the final force/step values, allowing an AI agent to decide whether to retry with more cycles.

An optimizer may also report `"status": "stalled"`: the energy stopped decreasing over the configured window (an energy plateau) while the configured force/step convergence criteria remained unmet. A stall is a distinct, non-converged outcome — it is never reported as `converged`. In `tsopt` it also stops the flatten/retry loop rather than repeating a non-progressing search; `opt --flatten` still runs its flatten loop, which exists to displace along the remaining imaginary modes and leave the plateau. When present, a `stop_reason` string records the energy range, window, and the failed criteria. A stall may be retried (e.g. from a perturbed geometry or with tighter step control); it is not an alias for `max_cycles` exhaustion or a generic failure.

## Subcommand schemas

### `sp`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"ok"` |
| `stage` | string | `"sp"` |
| `input` | string | Prepared public input path |
| `backend` / `model` | string / string \| null | Local MLIP provenance, mirrored to the common `mlip_*` fields |
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
| `status` | string | `"converged"`, `"not_converged"`, or `"stalled"` (energy plateau; see above) |
| `stop_reason` | string | Present only for a non-converged stop (stalled/stopped); records the energy plateau range/window and failed criteria |
| `energy_hartree` | float | Final energy (Hartree) |
| `n_opt_cycles` | int | Optimization cycles completed |
| `opt_mode` | string | `"grad"`, `"hess"`, `"lbfgs"`, or `"rfo"` |
| `backend` | string | MLIP backend (`"uma"`, `"orb"`, `"mace"`, `"aimnet2"`, or `"custom"` with `--calc-file`) |
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
| `rigid_projection` | object | Optional; present when `--flatten` runs. See [projection provenance](#rigid-projection-provenance). |

### `tsopt`

All fields from `opt`, plus:

| Field | Type | Description |
|-------|------|-------------|
| `n_imaginary_modes` | int | Number of imaginary frequencies |
| `imaginary_frequencies_cm` | float[] | Imaginary frequencies (cm⁻¹, negative) |
| `opt_mode` | string | `"rsprfo"` (default), `"rsirfo"`, `"trim"`, or `"dimer"` |
| `reference_mode_file` | string\|null | Advanced path-mode file supplied through `--ref-mode`; normally generated and passed by `all` |
| `safeguards` | object | Hessian-TS rejection/recovery diagnostics, including rejected mode-loss trials, exact saddle checks, final target-mode identity/overlap and higher-order MEP re-anchoring flag, stop reason, and bounded path-mode restart attempts |
| `rigid_projection` | object | Rigid-mode and exact-Hessian provenance; see [projection provenance](#rigid-projection-provenance) |

The `files` object may include `imaginary_mode_files` (list of vib file paths).
For Cartesian Hessian modes, `status: "converged"` requires the final exact
PHVA result to contain exactly one significant imaginary mode. Supported
internal-coordinate Hessian modes instead require one negative root in the
exact optimizer-space Hessian. Higher-order saddles and `n_imag=0` structures
are reported as `not_converged`. An energy-plateau `stalled` outcome (see the
general note above) is likewise never reported as `converged` and stops further
flatten/retry work. Convergence details are
available for rsirfo mode; dimer mode also reports `status: "converged"`,
`"not_converged"`, or `"stalled"`, but provides `n_opt_cycles` only and omits
the per-cycle force/step convergence keys and Hessian-mode `safeguards` object.

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
| `rigid_projection` | object | Rigid-mode and Hessian provenance; also written to `thermoanalysis.yaml` with `--dump` |

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
| `energy_first_hartree` | float | Energy of the first stitched-path endpoint; standalone IRC does not assign chemical identity |
| `energy_ts_hartree` | float | TS energy |
| `energy_last_hartree` | float | Energy of the last stitched-path endpoint; standalone IRC does not assign chemical identity |
| `endpoint_energy_orientation` | string | `"finished_first_to_finished_last"` |
| `energy_reactant_hartree` / `energy_product_hartree` | float | Legacy aliases for first / last, respectively; do not infer chemical R/P identity from these names |
| `forward_converged` | bool \| null | Forward IRC converged? `null` when the integrator did not expose the flag |
| `backward_converged` | bool \| null | Backward IRC converged? `null` when the integrator did not expose the flag |
| `forward_energy_increased` | bool \| null | Final forward step raised the energy |
| `backward_energy_increased` | bool \| null | Final backward step raised the energy |
| `backend` | string | MLIP backend |
| `charge` | int | System charge |
| `spin` | int | Spin multiplicity |
| `model` | string | Model identifier |
| `never_stop` | bool | Whether opt-in energy-rise/plateau bypass mode was enabled |
| `never_stop_energy_bypasses` | int | Number of energy-rise/plateau stop events actually bypassed |
| `n_freeze_atoms` | int | Frozen atoms |
| `solvent` | string | Implicit solvent or `"none"` |
| `bond_changes` | object | Directed first→last `{formed: [...], broken: [...]}` of element-prefixed 1-based atom-pair strings (e.g. `"C7-O12"`); key is omitted when the comparison fails or `finished_first.xyz`/`finished_last.xyz` are absent. |
| `bond_changes_direction` | string | `"finished_first_to_finished_last"` when `bond_changes` is present |
| `step_length` | float | IRC step length (Bohr) |
| `max_cycles` | int | Maximum IRC steps |
| `input_file` | string | Input filename |
| `files` | object | Trajectory files (XYZ plus available PDB/CIF companions) |
| `rigid_projection` | object | Rigid-mode and initial-Hessian provenance; see [projection provenance](#rigid-projection-provenance) |

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
| `charge` | int \| null | System charge; null for plot-only `scan3d --csv` |
| `spin` | int \| null | Spin multiplicity; null for plot-only `scan3d --csv` |
| `backend` | string \| null | MLIP backend; null for plot-only `scan3d --csv` |
| `model` | string \| null | Model identifier; null for plot-only `scan3d --csv` |
| `solvent` | string \| null | Implicit solvent or `"none"`; null for plot-only `scan3d --csv` because imported energies have no calculator provenance |
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

`path-search` has no `--out-json` flag. Once execution reaches its summary
writer it writes `summary.json` with the shared envelope (`command`,
`pdb2reaction_version`, `environment`) plus the fields below; early CLI/input
validation can fail before the file exists.

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"success"` / `"partial"` |
| `n_segments` | int | Recursive MEP segment count |
| `segments` | object[] | Per-segment `index`, `tag`, `kind`, `barrier_kcal`, `delta_kcal`, `bond_changes` (list of `{title: [entries]}` dicts; bridge segments emit `""`). |
| `energy_diagrams` | object[] | Per-segment labeled energy profiles (kcal/mol) |
| `mlip_backend` | string | Backend identifier |
| `mlip_model` | string \| null | Model identifier, recorded separately from the backend |
| `mlip_precision` | string \| null | Effective `fp32` / `fp64` token; null for custom calculators |
| `charge` | int | System charge |
| `spin` | int | Spin multiplicity |

See also the extended [`summary.json` section](#summary-json-path-search-all) for the additional fields that `all` layers on top.

### `dft`

> **Note:** `dft` writes `result.json` only on a successful SCF (exit 0), with
> `status: "converged"`. A
> non-converged SCF returns exit code 3 and skips `result.json`; SCF status
> is encoded by the exit code when no result exists. The generic
> `not_converged` status does not apply to `dft`;
> an unhandled exception, however, still writes the standard `error` envelope
> (`result.json` + `summary.json`).

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"converged"` |
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
| `n_atoms_raw` | int | Atoms in selected residues before backbone/truncation filtering (not the whole input structure) |
| `n_atoms_extracted` | int | Selected atoms kept after truncation, before cap-H addition |
| `total_charge` | float | Computed total charge |
| `protein_charge` | float | Protein charge |
| `ligand_total_charge` | float | Ligand charge sum |
| `ion_total_charge` | float | Ion charge sum |
| `ion_charges` | list | `[[name, charge], ...]` |
| `unknown_residue_charges` | object | `{resname: charge}` |
| `n_link_hydrogens` | int | Cap hydrogens added at carbon-parent truncation bonds |
| `exclude_backbone` | bool | Whether backbone atoms were excluded |
| `include_h2o` | bool | Whether crystallographic waters were included |
| `ligand_charge_input` | string | User-supplied `--ligand-charge` mapping |
| `center` | string | Center residue |
| `radius` | float | Extraction radius (angstrom) |
| `input_files` | string[] | Original input PDB/mmCIF paths |
| `files` | object | Output PDB / cluster filenames |

### `trj2fig`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"ok"` |
| `n_frames` | int | Number of trajectory frames |
| `min_energy_hartree` | float | Minimum energy across frames |
| `max_energy_hartree` | float | Maximum energy across frames |
| `energy_source` | string | `"trajectory_comment"` or `"mlip_recomputed"` |
| `backend` | string or null | MLIP backend only when frame energies were recomputed; null in comment-energy mode |
| `charge` / `multiplicity` | int or null | Resolved recomputation state, otherwise null |
| `solvent` / `solvent_model` | string or null | Recomputed-calculator solvent settings, otherwise null |
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
| `mlip_backend` | string | Backend identifier |
| `mlip_model` | string \| null | Model identifier, recorded separately from the backend |
| `mlip_precision` | string \| null | Effective `fp32` / `fp64` token; null for custom calculators |
| `charge` | int | System charge |
| `spin` | int | Spin multiplicity |
| `environment` | object | Hardware info |

The `all` command additionally includes:

| Field | Type | Description |
|-------|------|-------------|
| `rate_limiting_step` | object | Highest barrier at the highest method complete across every reactive segment (`DFT//MLIP_Gibbs` > `DFT` > `MLIP_Gibbs` > `MLIP` > `MEP`), with explicit `method` and raw `mep_barrier_kcal` |
| `overall_reaction_energy_kcal` | float | Overall reaction energy |
| `overall_reaction_energy_method` | string | Method of the overall reaction energy (`MEP`, `MLIP`, `MLIP_Gibbs`, `DFT`, or `DFT//MLIP_Gibbs`) |
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
