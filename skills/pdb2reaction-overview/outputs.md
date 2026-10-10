# Reading pdb2reaction outputs

`pdb2reaction all` writes `summary.log` (for people) and `summary.json` (for scripts) at the top of `--out-dir`; the structures to report are `segments/seg_NN/{reactant,ts,product}.*`. The run succeeded when `summary.json["scientific_status"]` is `"success"`: the console prints `Scientific status: success` under the last `====== Pipeline summary ======`, and each TS prints `[tsopt] Converged (n_imag=1).` Whether the endpoints are the intended R and P is still your check.

Standalone commands write a smaller `result.json` (with `--out-json`); its keys are on each command page in [pdb2reaction-cli](../pdb2reaction-cli/SKILL.md), not here.

## Output tree

```text
result_all/
├─ summary.log                 # text summary; sections [1]–[5]
├─ summary.json                # machine-readable results (below)
├─ mep_trj.xyz, mep_trj.pdb    # MEP over all segments (Endpoint and Scan-list modes)
├─ energy_diagram_MEP.png, energy_diagram_*_all.png, irc_plot_all.png
├─ segments/
│  └─ seg_NN/                  # one reaction step
│     ├─ reactant.*, ts.*, product.*   # structures to report, in the input format (--tsopt)
│     ├─ energy_diagram_*.png
│     ├─ ts/                   # TS optimization; vib/imag_*_trj.xyz animates the imaginary modes
│     ├─ irc/                  # IRC trajectories and irc_plot.png
│     ├─ endpoint_opt/         # endpoint optimizations
│     ├─ freq/{R,TS,P}/        # --thermo
│     └─ dft/{R,TS,P}/         # --dft
└─ _work/                      # intermediate files, including the TS candidates (HEI)
   ├─ models/                  # extracted clusters (-c)
   ├─ scan/                    # staged scan (-s)
   └─ path_opt/                # MEP and hei_seg_NN.* (path_search/ with --refine-path)
```

- `summary.log` sections: [1] MEP overview, [2] MEP barrier, reaction energy, and bond changes per segment, [3] post-processing per segment (n_imag under `TS imaginary freq:`), [4] energy diagrams, [5] output tree.
- mmCIF or very large PDB input also gets `.cif` files with the original chain IDs and residue numbers.
- TS-only mode has no MEP files and no `_work/path_opt/`; R, TS, and P go to `segments/seg_01/`.

## summary.json: top-level keys

| Key | What to read |
|---|---|
| `schema_version` | Schema version (`"4.0"`); check it before parsing |
| `execution_status` | `completed` or `failed`: whether the requested steps ran |
| `scientific_status`, `scientific_status_reasons` | `success`, `partial`, or `failed`, and the reasons when it is not `success` |
| `pipeline_mode` | `path-opt`, `path-search`, or `tsopt-only` |
| `charge`, `spin` | Charge and multiplicity used |
| `mlip_backend`, `mlip_model`, `mlip_precision` | Backend, exact model, and effective `fp32`/`fp64` (null for a custom calculator) |
| `references` | Methods the run used, as `{method, citation, doi}` records |
| `segments`, `post_segments` | Per-segment results (next section) |
| `rate_limiting_step` | `{segment, barrier_kcal, method, mep_barrier_kcal}` for the highest local barrier |
| `overall_reaction_energy_kcal`, `overall_reaction_energy_method` | R → P energy from the best complete all-segment diagram, and its method |
| `n_segments`, `n_segments_reactive` | All segments, and the reactive ones |
| `energy_diagrams` | Diagram energies and image metadata |
| `key_output_files`, `current_output_paths` | Index of the files of this run; stale files in a reused directory are excluded |
| `config`, `environment`, `command`, `pdb2reaction_version` | Effective settings, hardware, the command line, and the version |

- `success` means every requested stage converged; with `--tsopt`, every TS also has n_imag = 1. n_imag ≥ 2 gives `partial`; n_imag = 0 stops that segment before IRC, and the later segments are not post-processed. How the IRC stopped does not enter the status.
- The `references` set also appears at the end of `summary.log` and of the final console output, immediately before elapsed time.
- `rate_limiting_step.barrier_kcal` uses the highest method available for every reactive segment: `DFT//MLIP_Gibbs` > `DFT` > `MLIP_Gibbs` > `MLIP` > `MEP`. Always report it with its `method`. It is the highest local barrier, not a kinetic rate-limiting-step assignment, and it is absent when no segment is reactive.

## Per-segment keys

`segments[]` holds one MEP-level record per segment:

- `index` (1-based) and `tag` (a label whose number can differ from `index`).
- `kind`: `"seg"` for a reactive segment, `"kink"` for a segment with no covalent bond change, `"bridge"` for a short connecting path (its `tag` ends in `_bridge`), `"tsopt"` in TS-only mode. Reactive segments are those with `kind` `"seg"` or `"tsopt"`.
- `barrier_kcal` and `delta_kcal`: barrier and reaction energy on the MEP before TS optimization (`null` for a kink); in TS-only mode, TS − R.
- `bond_changes`: see [Bond changes](#bond-changes).

`post_segments[]` appears with `--tsopt`. Match it to `segments` by `index`, not by list position.

- `tsopt`: `optimization_status`, `hessian_status`, `saddle_validation` (`first_order`, `higher_order`, `no_imaginary`, or `unavailable`), `n_imaginary_modes` (n_imag), `imaginary_frequencies_cm`, and `n_opt_cycles` / `max_cycles`.
- `mlip`: R, TS, and P energies with `barrier_kcal`, `delta_kcal`, and `energies_kcal`. `gibbs_mlip` (`--thermo`), `dft` (`--dft`), and `gibbs_dft_mlip` (both) have the same shape. When DFT fails for any state, `dft` is `{"scientific_status": "failed", "failed_states": [...]}` and no DFT diagram is written.
- `mep_barrier_kcal`, `mep_delta_kcal`: the MEP values of the same segment.
- `irc`: per-direction stop diagnostics; there is no IRC pass/fail verdict.
- `endpoint_assignment`: how the IRC ends were matched to R and P.
- `endpoint_opt`: `reactant` and `product`, each with `optimization_status`, `n_opt_cycles`, `max_cycles`, and any `stop_reason`. In Endpoint and Scan-list modes, `connectivity_validated` (with `connectivity.match_matrix`) tells whether the optimized R and P kept the bond topology of the MEP ends. `scientific_status` `success` does not check this, and `bond_changes` ([R/TS/P paths](#rtsp-paths)) cannot show what the TS connects; treat `false` as a TS that connects other states.
- `thermo_symmetry`: point group and rotational symmetry per state.

These block names never change with the backend: `mlip`, `gibbs_mlip`, `dft`, and `gibbs_dft_mlip`. Energy-diagram filenames use `MLIP` for every backend, while top-level `mlip_backend` / `mlip_model` / `mlip_precision` record the exact provenance.

## R/TS/P paths

- Report `segments/seg_NN/{reactant,ts,product}.*`: the optimized TS, and the IRC ends after endpoint optimization.
- `segments/seg_NN/` exists only for reactive segments, and NN is `segments[].index`, so after `--refine-path` the list can start at `seg_02` or skip numbers. List `segments/` or read `post_segments[].index` instead of assuming `seg_01`, take the barrier TS from the segment in `rate_limiting_step.segment` (not `post_segments[0]`), and copy R, TS, and P of one segment together.
- In Endpoint and Scan-list modes, the IRC ends are matched to the MEP's left and right states by bond pattern, then RMSD.
- In TS-only mode, the higher-energy IRC end is named the reactant (the left end on an exact tie). This names the direction, not the chemical direction of the reaction; `endpoint_assignment` records `chemical_direction_known: false`. The barrier from P is `barrier_kcal − delta_kcal`.
- The raw IRC ends before optimization are in `segments/seg_NN/structures/{reactant,product}_irc.*`; use them only to debug a difference between the IRC and the endpoint optimization.
- `bond_changes` comes from the MEP ends in Endpoint and Scan-list modes, and from the optimized R and P in TS-only mode.

## Reading summary.json with Python

```python
import json

d = json.load(open("result_all/summary.json"))
if d.get("schema_version") != "4.0":
    raise RuntimeError(f"unsupported summary schema: {d.get('schema_version')!r}")
if d.get("scientific_status") != "success":
    raise RuntimeError(f"requested workflow is incomplete: {d.get('scientific_status_reasons', [])}")

# Barriers after TS optimization (post_segments is matched by index)
for ps in d.get("post_segments", []):
    mlip = ps.get("mlip") or {}
    if mlip.get("barrier_kcal") is None:
        continue
    line = (f"seg_{ps['index']:02d}: ΔE‡ = {mlip['barrier_kcal']:.1f}, "
            f"ΔE = {mlip['delta_kcal']:.1f} kcal/mol")
    gibbs = ps.get("gibbs_mlip") or {}
    if gibbs.get("barrier_kcal") is not None:
        line += f", ΔG‡ = {gibbs['barrier_kcal']:.1f} kcal/mol"
    line += f", n_imag = {(ps.get('tsopt') or {}).get('n_imaginary_modes')}"
    print(line)

# MEP barrier before TS optimization: keep the "MEP" label when you print it
for seg in d["segments"]:
    if seg.get("kind") in ("seg", "tsopt") and seg.get("barrier_kcal") is not None:
        print(f"seg_{seg['index']:02d}: MEP barrier = {seg['barrier_kcal']:.1f} kcal/mol")

# Highest local barrier, labelled by its method
rls = d.get("rate_limiting_step")
if rls:
    print(f"highest barrier: seg_{rls['segment']:02d}, "
          f"{rls['barrier_kcal']:.1f} kcal/mol ({rls['method']})")
```

## Bond changes

`segments[i]["bond_changes"]` is a list of single-key dicts; the key names the change and `(k)` counts the entries:

```json
[
  {"Bond formed (1)": ["C508-C567 : 3.166 Å --> 1.675 Å"]},
  {"Bond broken (1)": ["S507-C508 : 1.798 Å --> 3.459 Å"]}
]
```

- "Bond formed" lists bonds present in P but not in R; "Bond broken" lists bonds present in R but not in P. A section with no change is `{"Bond formed": ["None"]}`, without `(k)`.
- A bond is counted at up to 1.20 × the sum of covalent radii (margin 0.05). `bond-summary` in [utilities.md](../pdb2reaction-cli/utilities.md) uses the same rule.
- Standalone `irc` writes a flat `{"formed": [...], "broken": [...]}` dict in its `result.json`, with `bond_changes_direction`.
- An unexpected or long list is not a verdict. Look at the endpoint structures, then check the segment with TS optimization, n_imag, and IRC.

## Energy diagrams

| File | Written with |
|---|---|
| `energy_diagram_MEP.png` (top) | Endpoint and Scan-list modes |
| `seg_NN/energy_diagram_MLIP.png` | `--tsopt` |
| `seg_NN/energy_diagram_G_MLIP.png` | `--thermo` |
| `seg_NN/energy_diagram_DFT.png` | `--dft` |
| `seg_NN/energy_diagram_G_DFT_plus_MLIP.png` | `--dft` and `--thermo` |
| `energy_diagram_*_all.png` (top) | The same diagram over all segments |

Energies are in kcal/mol relative to R. A PNG is written only when the energies are finite and the export succeeds; `summary.json["energy_diagrams"]` keeps the energies either way, so check `image_written` there before you use a PNG. In TS-only mode the diagrams are in `segments/seg_01/`.

Read a barrier from `barrier_kcal` of its level (`mlip`, `gibbs_mlip`, `dft`, `gibbs_dft_mlip`), which is TS − R, or from the TS-labelled entry of `energy_diagrams`. Do not compute `max(energies) − E(R)`: after thermal or DFT corrections, P or a later state can lie above the TS.

To combine energies from several runs into one diagram, use `energy-diagram` ([utilities.md](../pdb2reaction-cli/utilities.md)):

```bash
pdb2reaction energy-diagram -i "[0.0, 21.5, -0.7, 2.2, -18.2]" --label-x "['R','TS1','IM','TS2','P']" -o my_diagram.png
```

## Failed runs

When `execution_status` is `failed` or `scientific_status` is not `success`:

1. Read `scientific_status_reasons` and `summary.log`; an early stop is shown as `Pipeline stop`.
2. Look in the stage directories under `segments/seg_NN/`: `ts/` and `irc/` hold a `result.json`, and an endpoint optimization that raised an error writes `endpoint_opt/failure.json`.
3. Optimizer logs are not written by default; rerun with `--dump` to keep the optimizer trajectories.

Partial outputs are kept. `_work/path_opt/seg_NNN_<tag>/` (`_work/path_search/` with `--refine-path`) exists for every segment that finished the MEP. `segments/seg_NN/` is created when post-processing starts, so its presence is not a sign of success; use the statuses in `summary.json`.

A barrier or reaction energy of hundreds to thousands of kcal/mol usually means that one endpoint optimization diverged, not a summary bug. Compare the Hartree energy on the second line of `segments/seg_NN/structures/{reactant,product}_irc.xyz` with that of `{reactant,product}.xyz` in the same directory. If the raw IRC end is sensible and the optimized one is not, rerun only that end with `opt` from its `*_irc.xyz` file with the same charge, spin, and frozen atoms; `summary.json["freeze_atoms"]` is 0-based, so add 1 for `--freeze-atoms`.

Read the exit code with the status: 0 is `success` or `partial` (tell them apart by `scientific_status`); 1 is non-convergence, no usable result, a runtime exception, or an output failure; 2 is invalid input, arguments, or configuration (including invalid YAML), so fix the command instead of resubmitting it; 130 is an interrupt. A run that fails while the options or inputs are checked can stop before `summary.json` exists; treat a nonzero exit code as a failure and read stderr.

## Next step

- [all.md](../pdb2reaction-cli/all.md) and its mode pages: the command that writes `summary.json`, and resuming a failed segment.
- [tsopt.md](../pdb2reaction-cli/tsopt.md), [irc.md](../pdb2reaction-cli/irc.md), [freq.md](../pdb2reaction-cli/freq.md), [dft.md](../pdb2reaction-cli/dft.md): `result.json` of each standalone stage.
- [ts-strategy.md](ts-strategy.md): what to do when n_imag or the endpoints are wrong.
- Docs: [json-output](../../docs/json-output.md), [output-layout](../../docs/output-layout.md), [all](../../docs/all.md#reading-the-run-status).
