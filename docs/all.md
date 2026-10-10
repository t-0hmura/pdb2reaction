# `all` (end-to-end workflow)

`all` runs the whole workflow in one command: it extracts the active-site model and builds the minimum energy path (MEP). When asked, it also optimizes the transition state (TS) of each reaction step and runs the intrinsic reaction coordinate (IRC), frequency, and DFT calculations on it.

Without `--tsopt`, the run ends with TS candidates: the highest-energy image (HEI) of each MEP segment. The default backend is **UMA**, Meta's pretrained [machine-learning interatomic potential (MLIP)](backends.md).

---

## What it is for

What you pass selects the mode:

* **Path and energy diagram from R and P**: give two or more structures in reaction order (reactant, intermediates, product); `all` finds the MEP between each neighbouring pair and draws the energy diagram.
* **Path from a reactant alone**: give one structure and the bonds to form or break with `-s`; a staged scan makes the intermediates, and the MEP search runs through them.
* **Check one TS candidate (TS-only mode)**: give one structure with `--tsopt` and no `-s`; `all` optimizes the TS and runs IRC from it. The TS is confirmed when n_imag = 1 and the IRC ends at the intended R and P.

---

## Examples

The examples use the GPP C6-methyltransferase BezA ([Tsutsumi et al., *Angew. Chem. Int. Ed.* 2022, 61, e202111217](https://doi.org/10.1002/anie.202111217)); the full scripts are in [`examples/`](https://github.com/t-0hmura/pdb2reaction/tree/main/examples). `1.R.pdb` (reactant), `2.IM.pdb` (intermediate), and `3.P.pdb` (product) are full structures with every hydrogen; your own structures need hydrogens too. Examples 1–3 are walked through, with how to check the results, in [Quickstart: `all`](quickstart-all.md), [Quickstart: `--scan-lists`](quickstart-scan.md), and [Quickstart: TS-only mode](quickstart-tsopt.md).

### 1. MEP with TS optimization, thermochemistry, and DFT

`-c` names the extraction centers, and `-l` gives the charges of the non-standard residues.

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo --dft --out-dir ./result_mep
```

Every requested stage finished when the console prints `[tsopt] Converged (n_imag=1).` for each TS and `Scientific status: success` under the last `====== Pipeline summary ======`; `result_mep/summary.json` holds the same values. Then check the endpoints as in [Reading the run status](#reading-the-run-status). The optimized structures are in `result_mep/segments/seg_NN/`.

### 2. Path from the reactant by a staged scan

Stage 1 brings the methyl carbon of SAM (CS1) to C7 of GPP (1.60 Å), and stage 2 moves H11 of GPP onto OE2 of Glu186 (0.90 Å).

```bash
pdb2reaction all -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -s '[("CS1 SAM 320","GPP 321 C7",1.60)]' '[("GPP 321 H11","GLU 186 OE2",0.90)]' \
    --tsopt --thermo --out-dir ./result_scan
```

The targets inside one literal move together in one stage. Literals given in a row run as successive stages, each starting from the end of the one before, and the stage ends become the inputs of the MEP search. Give `-s` once and list every literal after it. To decide how to split a reaction, see {ref}`Decide how to split the reaction <mechanism-split>`. In a PDB with an empty chain field, write an atom as its residue name, residue number, and atom name in any order (`"CS1 SAM 320"`); with chains, write `A:SAM:320:CS1`. All accepted forms are in {ref}`Scan-list spec <scan-list-spec>`.

### 3. Check a TS candidate (TS-only mode)

One input with `--tsopt` and no `-s` skips the MEP search.

```bash
pdb2reaction all -i TS_candidate.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo --dft
```

The optimized R, TS, and P are written to `result_all/segments/seg_01/`.

### 4. Resume post-processing from a segment

To redo the post-processing from segment N, repeat the original command with the same inputs, extraction, path, and calculator options and the same `--out-dir`, and add `--resume-segment N`. Post-processing options such as `--tsopt-max-cycles` may change.

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo --dft \
    --resume-segment 1 --out-dir ./result_mep
```

The segments before N are kept; the post-processing from segment N onward, the summary, and the diagrams are written again.

---

## How it works

```text
Full structure(s) (PDB / mmCIF / XYZ / GJF)
  ├─ (with -c) active-site extraction: extract
  │   └─ active-site model(s)
  ├─ (one structure with -s) staged scan: scan
  │   └─ stage ends as intermediates
  ├─ MEP search: path-opt (default) or path-search (--refine-path)
  │   └─ mep_trj.xyz and energy_diagram_MEP.png
  └─ (with --tsopt) TS optimization and IRC: tsopt → irc
      ├─ (with --thermo) frequencies and thermochemistry: freq
      └─ (with --dft) DFT single points: dft
```

1. **Preparing the input**: for a PDB with alternate locations (altloc), `all` keeps one label per residue, the one with the highest mean occupancy; for a PDB with blank element columns, it fills them in. With `-c`, it cuts out the active-site model around the given residues and caps the cut bonds with hydrogens.
2. **Building the path**: the input structures are optimized first (`--preopt`). With `-s`, the staged scan makes the intermediates. `path-opt` then finds the MEP between each neighbouring pair by GSM (growing string method) or DMF (direct max flux); with `--refine-path`, the recursive `path-search` refines the path and splits it into steps where bonds change. The HEI of each step is its TS candidate.
3. **Optimizing the TS** (`--tsopt`): each HEI is optimized by RS-P-RFO (restricted-step partitioned rational function optimization) by default, and the final Hessian gives n_imag.
4. **Following the IRC**: from the TS, the IRC is traced in both directions with EulerPC (Euler predictor–corrector), and both ends are optimized to minima. These become the R and P of the segment.
5. **Thermochemistry and DFT**: `--thermo` runs `freq` on R, TS, and P for the Gibbs energy, and `--dft` adds DFT single points on the same structures. Each adds its own energy diagram.

`all` continues from the TS to IRC only when the TS optimization converged, its final Hessian was computed, and n_imag ≥ 1:

| TS result | What `all` does |
| --- | --- |
| Converged, n_imag = 1 | Runs IRC and optimizes both IRC ends. |
| Converged, n_imag ≥ 2 | Runs IRC with a warning along the imaginary mode that best matches the MEP direction (the lowest one when none matches). The result is `partial`. |
| Converged, n_imag = 0 | Stops before IRC. |
| Not converged (cycle limit or `--stop-plateau`), the final Hessian skipped with `--skip-final-freq`, or a failed Hessian | Stops before IRC. |

When `all` stops before IRC, the result is not `success`; the TS files stay in `segments/seg_NN/ts/`, and the later segments are not post-processed. The full table of how a TS optimization can end is in [`tsopt` → Reading the TS result](tsopt.md#reading-the-ts-result).

If one endpoint optimization does not converge, the result is `partial` and `segments/seg_NN/endpoint_opt/` is kept for inspection. If an endpoint optimization fails with an error, the error is written to `segments/seg_NN/endpoint_opt/failure.json`, and that segment stops before the frequency and DFT stages; the TS and IRC structures are kept.

---

## Reading the run status

A successful TS optimization gives one imaginary mode along the reaction coordinate (n_imag = 1). Even if the IRC does not converge, the result is usable when the endpoint optimizations reach the intended R and P.

Read the outcome in three places:

* **Console**: each TS optimization ends with `[tsopt] Converged (n_imag=1).` when it converged with one imaginary mode. The `====== Pipeline summary ======` block prints `Execution status:` and `Scientific status:`. When the result is not `success`, `RESULT WARNING:` lines give the reasons.
* **`summary.log`**: the header shows `Pipeline mode` (`MEP`, `Scan`, or `TS-only`) and both statuses. Section [1] is the MEP overview; [2] lists the barrier ΔE‡, the reaction energy ΔE, and the bond changes of each segment on the MEP; [3] gives the post-processing of each segment, with n_imag under `TS imaginary freq:`; [4] tabulates the energy diagrams; [5] shows the output tree.
* **`summary.json`**: `scientific_status` holds `success`, `partial`, or `failed`, and `scientific_status_reasons` holds the [reasons](json-output.md#execution-and-requested-stage-completion). n_imag of each TS is `post_segments[].tsopt.n_imaginary_modes`.
  * **Barriers in `summary.json`**: with `--tsopt`, the barrier of each segment is `post_segments[].mlip.barrier_kcal` (MLIP energy of the optimized TS minus R). With `--thermo` and `--dft`, the same `barrier_kcal` is also under `gibbs_mlip`, `dft`, and `gibbs_dft_mlip`. `segments[].barrier_kcal` is the barrier on the MEP before TS optimization, or TS − R in TS-only mode.

`success` means that every requested stage converged; with `--tsopt`, it also means that every TS has n_imag = 1. Whether the endpoints are the intended R and P is for you to check: compare the bond changes in section [2] of `summary.log` and the structures `segments/seg_NN/reactant.*` and `product.*` with the R and P you intended. If n_imag ≠ 1 or the endpoints are not the intended ones, see {ref}`When the TS search fails <ts-search-fails>`.

---

## Output files

`all` writes these files to `--out-dir`:

```text
result_all/
├─ summary.log                  # Text summary
├─ summary.json                 # Machine-readable results (always written; all has no --out-json)
├─ mep_trj.xyz                  # MEP trajectory over all segments
├─ mep_trj.pdb                  # Same trajectory as PDB
├─ mep_trj.cif                  # Same trajectory as mmCIF (mmCIF or very large PDB input)
├─ mep_w_ref.pdb                # MEP merged into the full input (--write-ref-merge)
├─ energy_diagram_MEP.png       # MEP energy profile over all segments
├─ energy_diagram_*_all.png     # R → TS → P diagrams over all segments (--tsopt, --thermo, --dft)
├─ irc_plot_all.png             # IRC profiles over all segments (--tsopt)
├─ segments/
│  └─ seg_NN/                   # One reaction step: seg_01, seg_02, ...
│     ├─ reactant.*             # Optimized R, TS, and P in the input format (--tsopt)
│     ├─ ts.*
│     ├─ product.*
│     ├─ energy_diagram_*.png   # R → TS → P diagrams of this step
│     ├─ ts/                    # TS optimization; vib/imag_*_trj.xyz animates the imaginary modes
│     ├─ irc/                   # IRC trajectories and irc_plot.png
│     ├─ endpoint_opt/          # Endpoint optimizations (kept with --dump, or when an endpoint did not converge or failed)
│     ├─ freq/{R,TS,P}/         # Frequencies and thermochemistry (--thermo)
│     └─ dft/{R,TS,P}/          # DFT single points (--dft)
└─ _work/                       # Intermediate files, including the TS candidates (HEI)
   ├─ models/                   # Extracted models, model_<input>.pdb (with -c)
   ├─ scan/                     # Staged scan (with -s)
   └─ path_opt/                 # MEP search and hei_seg_NN.* (path_search/ with --refine-path)
```

* **Structures to report**: cite `segments/seg_NN/reactant.*`, `ts.*`, and `product.*`. The subdirectories of `seg_NN/` hold the files of each stage.
* **TS-only mode**: there is no MEP search, so the MEP files and `_work/path_opt/` are absent; R, TS, and P go to `segments/seg_01/`.

The energy diagrams are named by method:

| File | Written when | Content |
| --- | --- | --- |
| `energy_diagram_MEP.png` | The MEP search finishes | MEP energy profile over all segments |
| `energy_diagram_MLIP.png` | `--tsopt` | R → TS → P, MLIP energy |
| `energy_diagram_G_MLIP.png` | `--thermo` | R → TS → P, MLIP Gibbs energy |
| `energy_diagram_DFT.png` | `--dft` | R → TS → P, DFT energy on the MLIP geometries |
| `energy_diagram_G_DFT_plus_MLIP.png` | `--dft` and `--thermo` | R → TS → P, DFT energy plus the MLIP thermal correction |
| `energy_diagram_*_all.png` | Same as the diagram without `_all` | The same diagram over all segments, at the top of the output directory |
| `irc_plot.png` (in `seg_NN/irc/`), `irc_plot_all.png` | `--tsopt` | IRC energy profile of one segment, and of all segments |

Energies in the diagrams are in kcal/mol relative to the first state (the reactant).

---

## Main options

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path(s) | (required) | Two or more structures in reaction order, or one structure with `-s` or `--tsopt` (`.pdb`, `.cif`, `.xyz`, `.gjf`). Give several files after one `-i`, or repeat `-i` |
| `-c, --center` | text | `None` | Extraction centers, normally the substrate and catalytic residues: residue names (`'SAM,GPP'`), residue IDs (`'A:123,B:456'`), or chain-qualified names (`'A:SAM'`, `'A:SAM:123'`). Omit to use the full input |
| `-l, --ligand-charge` | text | `None` | Charges of non-standard residues (e.g. `'SAM:1,GPP:-3'`), or their total charge as one number (total ligand charge). PDB/mmCIF input only |
| `-q, --charge` | integer | `None` | Total charge. With `-c`, it comes from the extracted model; without `-c`, required unless `-l` is given or the input is `.gjf`. An explicit value overrides the derived one with a warning |
| `-m, --multiplicity` | integer | `1` | Spin multiplicity (2S+1) |
| `-b, --backend` | text | `uma` | Calculator backend (`uma`, `orb`, `mace`, `aimnet2`, `dft`) |
| `-r, --radius` | float | `2.6` | Extraction cutoff (Å) around the center atoms. `0` keeps only the `-c` and `--selected-resn` residues |
| `--selected-resn` | text | `""` | Residues to include without radius expansion, in the same forms as `-c` |
| `-s, --scan-lists` | text | `None` | Staged scan targets for one input, one literal per stage (e.g. `'[("A:SAM:320:CS1","A:GPP:321:C7",1.60)]'`; format in {ref}`Scan-list spec <scan-list-spec>`) |
| `--tsopt/--no-tsopt` | flag | `False` | Optimize the TS of each segment and run IRC |
| `--thermo/--no-thermo` | flag | `False` | Frequencies and thermochemistry on R, TS, and P (needs `--tsopt`) |
| `--dft/--no-dft` | flag | `False` | DFT single points on R, TS, and P (needs `--tsopt`) |
| `--refine-path/--no-refine-path` | flag | `False` | Run the recursive `path-search` instead of one `path-opt` per pair |
| `--mep-mode` | `gsm` / `dmf` | `gsm` | MEP method: GSM or DMF |
| `--opt-mode` | `grad` / `hess` | `grad` | Optimizer for the single-structure optimizations and the scan: `grad` = L-BFGS, `hess` = RFO |
| `--opt-mode-post` | `grad` / `hess` | `hess` (the `--opt-mode` value when `--opt-mode` is given on the command line) | Optimizer for the TS and the endpoints after IRC: `grad` = Dimer for the TS and L-BFGS for the endpoints, `hess` = RS-P-RFO for the TS and RFO for the endpoints |
| `--preopt/--no-preopt` | flag | `True` | Optimize the input structures before the scan and the MEP search |
| `--flatten/--no-flatten` | flag | `False` | Remove extra imaginary modes left after the TS optimization |
| `--stop-plateau/--no-stop-plateau` | flag | `False` | Stop an optimization when the energy stops changing before convergence; the run is reported as stalled, not converged |
| `--tsopt-max-cycles` | integer | `100000` | Cycle limit of the TS optimization |
| `--resume-segment` | integer | `None` | Redo the post-processing from segment N, reusing the MEP in `--out-dir` (example 4) |
| `--dry-run/--no-dry-run` | flag | `False` | Check the options and print the plan without running a calculation. With `-c`, extraction runs in a temporary directory to check the charge |
| `-o, --out-dir` | path | `./result_all/` | Output directory |

For every option, run `pdb2reaction all --help-advanced` or see the [generated CLI reference](reference/commands/all.md).

> **Note:** In YAML (`--config`), you can set what the options above do not cover. See [YAML Reference](yaml-reference.md) for the sections and keys.

---

## Notes

* **`--dft` and `-b dft`**: they cannot be used together, and the run stops with an error at startup. To add DFT single points after a `-b dft` run, run `pdb2reaction sp -b dft` or `pdb2reaction dft` as a separate job.
* **Cost of `--dft`**: memory use depends on the structure, basis, functional, precision, and software stack. Try a representative structure on the target node and watch the peak memory. For a large model, finish the MLIP run first and run the DFT single points as a separate job.
* **R and P in TS-only mode**: the higher-energy IRC end is named the reactant (on an exact tie, the left end). The names, the file names, the barrier, and the reaction energy follow this energy order, not a known chemical direction; the barrier from P is `barrier_kcal − delta_kcal`. `summary.json` records the rule under `endpoint_assignment`, with `chemical_direction_known: false`.
* **`summary.log` in TS-only mode**: section [1] is the TS and IRC overview, and [2] comes from the optimized TS and endpoints.
* **Thermochemistry file**: with `--thermo`, `thermoanalysis.yaml` is kept even under [`--no-dump`](reference/commands/all.md), because `all` reads the thermochemistry from it.
* **Extraction radius**: `-r 0` disables radius-based expansion, so the model starts from the residues selected by `-c` and `--selected-resn`. Structural safeguards can still add a disulfide partner or the backbone of an adjacent residue. A zero radius is evaluated internally as 0.001 Å.
* **Without `-c`**: extraction is skipped, and the full input structures go to the MEP search, `tsopt`, `freq`, and `dft`. One structure still needs `-s` or `--tsopt`.
* **Input formats**: with `-c`, the input must be PDB or mmCIF; without `-c`, XYZ and GJF are accepted too. All structures of one run must have the same atoms in the same order.
* **Charge and multiplicity**: with `-c`, the total charge is the sum over the extracted model: built-in values for amino acids, ions, and water, `-l` for the other residues, and 0 for residues not listed in `-l`. Without `-c`, it comes from `-l` applied to the input, or from the `.gjf` header. The multiplicity is `-m`, otherwise the `.gjf` header, otherwise 1. See {ref}`Charge specification <charge-specification>`.
* **Separately prepared structures**: when the input structures were prepared independently, their differences outside the reaction coordinate enter the barrier. Compare the structures before reading the barrier. For two mechanisms of the same composition, use one common atom set and atom order for both paths.
* **`--write-ref-merge`**: writes the path merged back into the original full input, for inspection: `mep_w_ref*` in the output directory and `hei_w_ref_seg_NN.pdb` in `_work/path_search/`. It needs `--refine-path`, `-c`, and PDB or mmCIF input.
* **`--resume-segment`**: it needs `--tsopt`, `--thermo`, or `--dft`, and cannot be combined with `--dry-run`. The run stops with an error when the saved inputs and MEP do not match the command.

### Comparing a mutant with the wild type

Within one path, every structure has the same atoms in the same order. A mutant and the wild type (WT) differ in residues and often in atom count, so their total energies cannot be subtracted. Compare the barriers computed within each system instead:

`ΔΔG‡ = (G_TS − G_R)_mutant − (G_TS − G_R)_WT`

* Select the same residue positions and the same boundary and cap rules for both models, so that the mutation is the only designed difference. Two independent radius-based extractions can differ, because a boundary residue may enter one model and not the other; compare the two selections.
* Use the same protonation rules, charge assignment, backend and model, precision, restraints, and thermochemistry settings. If the mutation changes a protonation state or a formal charge, the total charges differ; do not force the same `-q` on both.

The two runs use the same options except for the input and the output directory. Give R and P of each system (Endpoint mode), so that R is the chemical reactant; `G_TS − G_R` is `post_segments[].gibbs_mlip.barrier_kcal`:

```bash
pdb2reaction all -i wt_R.pdb wt_P.pdb -c 'SAM,GPP,MG' -l 'GPP:-3,SAM:1' --tsopt --thermo -o result_wt
pdb2reaction all -i mutant_R.pdb mutant_P.pdb -c 'SAM,GPP,MG' -l 'GPP:-3,SAM:1' --tsopt --thermo -o result_mutant
```

---

## See also

* [extract](extract.md) — extraction of the active-site model
* [scan](scan.md) — staged scans of distances, angles, and dihedrals
* [path-opt](path-opt.md) — one MEP between two structures (GSM / DMF)
* [path-search](path-search.md) — recursive MEP search that splits the path into steps
* [tsopt](tsopt.md) — TS optimization
* [irc](irc.md) — IRC from a TS
* [freq](freq.md) — vibrational analysis and thermochemistry
* [dft](dft.md) — DFT single points
* [Refine an MLIP TS with DFT](dft-backend.md) — `-b dft` and `--dft`
* [Tips for studying reaction mechanisms](mechanism-tips.md) — splitting the reaction, checking the TS, and what to try when it fails
* [Troubleshooting](troubleshooting.md) — what to do when a run fails
* [Getting Started](getting-started.md) — the shortest run and what to read next
