# `scan` (restrained coordinate scan)

`scan` drives chosen distances, angles, or dihedrals of one structure step by step with harmonic restraints, relaxing every other degree of freedom at each step, and so builds a candidate reaction path from a single structure. The coordinates in one literal (or one YAML stage) move together as one **stage**; several literals run as stages in sequence, each starting from the relaxed end of the previous one.

---

## What it is for

* **A path from one structure**: drive the reacting bonds of a reactant to get intermediate- and product-like structures for [`path-search`](path-search.md).
* **Testing the order of events**: drive bond formation and proton transfer in one stage or in separate stages, and compare the energy profiles.
* **Running the scan step of `all` on its own**: repeat the scan that [`all`](all.md) runs for `-s`, with other step sizes or restraints.

The default backend is **UMA** (Meta); `-b/--backend` also selects **ORB**, **MACE**, **AIMNet2**, or DFT (`dft`). For an energy grid over two or three independent coordinates, use [`scan2d`](scan2d.md) or [`scan3d`](scan3d.md).

---

## Examples

The examples use `input.pdb`, the cluster model cut from the bundled enzyme structure with [extract](extract.md), and take its charge from `-l 'SAM:1,GPP:-3'`.

```bash
pdb2reaction extract -i examples/1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' -o input.pdb
```

Its PDB has an empty chain column, so an atom is written as residue name, residue number, and atom name, in any order, separated by commas or spaces (`'SAM,320,CS1'` or `'CS1 SAM 320'`).

### 1. From a YAML spec

Write the stages in a file and add `--out-json` for a summary.

```yaml
# scan.yaml
stages:
  - [["SAM,320,CS1", "GPP,321,C7", 1.60]]
  - [["GPP,321,H11", "GLU,186,OE2", 0.90]]
```

```bash
pdb2reaction scan -i input.pdb -l 'SAM:1,GPP:-3' -m 1 -s scan.yaml --out-json -o ./result_scan
```

The console prints the covalent-bond changes of each stage and ends with `====== Scan summary ======`. `result_scan/result.json` gives `scientific_status`, with the reasons in `scientific_status_reasons` when it is not `success`.

### 2. Inline literal

A short single-stage scan can be given on the command line.

```bash
pdb2reaction scan -i input.pdb -l 'SAM:1,GPP:-3' -m 1 -s '[("SAM,320,CS1","GPP,321,C7",1.60)]'
```

### 3. Two coordinates in one stage

Coordinates in the same literal move together (a concerted step).

```bash
pdb2reaction scan -i input.pdb -l 'SAM:1,GPP:-3' -m 1 \
    -s '[("CS1 SAM 320","GPP 321 C7",1.60),("GPP 321 H11","GLU 186 OE2",0.90)]' -o ./result_concerted
```

### 4. Two stages in sequence

Give several literals after one `-s`; stage 2 starts from the relaxed result of stage 1.

```bash
pdb2reaction scan -i input.pdb -l 'SAM:1,GPP:-3' -m 1 \
    -s '[("SAM,320,CS1","GPP,321,C7",1.60)]' '[("GPP,321,H11","GLU,186,OE2",0.90)]' -o ./result_staged
```

### 5. Bidirectional scan

A [4-tuple](#bidirectional-scan-4-tuple) scans one distance in both directions from the input geometry.

```bash
pdb2reaction scan -i input.pdb -l 'SAM:1,GPP:-3' -m 1 -s '[("SAM,320,CS1","GPP,321,C7",1.60,3.00)]'
```

### 6. Dump trajectories

Add `--dump` to keep the optimizer trajectory of every step.

```bash
pdb2reaction scan -i input.pdb -l 'SAM:1,GPP:-3' -m 1 -s scan.yaml --dump -o ./result_scan_dump
```

---

## How it works

1. **Reading the structure**: the {ref}`charge <charge-specification>` comes from `-q` or `-l`. With `--preopt`, the structure is first optimized without restraints; if that does not converge, the input geometry is used.
2. **Splitting each stage into steps**: for every coordinate, `scan` takes the change Δ = target − current and divides the stage into N = ceil(max(|Δ| / h)) steps, where h is `--max-step-size` (Å) for distances, `--max-angle-step-size` for angles, and `--max-dihedral-step-size` for dihedrals (degrees). Each coordinate moves by Δ / N per step, so all coordinates of a stage arrive together.
3. **Restrained relaxation**: at each step, a harmonic restraint E = ½ k (q − q_target)² holds every scanned coordinate q at its step target (k from `--restraint-k`), and the rest of the structure is relaxed with L-BFGS (`--opt-mode grad`, default) or RFO (`--opt-mode hess`). The energy written for each step is computed with the restraints removed.
4. **End of the stage**: with `--endopt`, the last structure of the stage is optimized once more without restraints. `scan` then compares the first and last structures of the stage for covalent-bond changes and writes the stage result.
5. **Next stage**: the next stage starts from this result. After the last stage, the trajectories of all stages are joined into one file.

### Bidirectional scan (4-tuple)

A range `(i, j, low, high)` instead of a target `(i, j, target)` scans in both directions from the input geometry. It expands into two stages:

1. **Pass 1**: drive `i`–`j` from the current distance toward `low`.
2. **Pass 2**: restore the input geometry and drive `i`–`j` toward `high`.

The joined trajectory runs `low → input geometry → high`, a continuous path through the starting structure. Angle ranges `(i, j, k, low, high)` and dihedral ranges `(i, j, k, l, low, high)` are scanned the same way.

(section-bond)=
### Bond-change detection

Let T be the sum of the covalent radii of two atoms scaled by `bond_factor` (default `1.20`). The atoms count as bonded when their distance is at most T − `margin_fraction` × T (default `0.05`). A pair is reported as formed or broken only when its distance changed by at least `delta_fraction` × T (default `0.05`). `path-search` uses the same rules; the keys are in the YAML [`bond`](yaml-reference.md#bond) section.

---

(scan-checking-result)=
## Reading the result

| Where | What to check |
| --- | --- |
| Console, each stage | `[stage k] Covalent-bond changes (start vs final): Yes` with the formed and broken bonds listed, or `No` with `(no covalent changes detected)` |
| Console, end of run | `====== Scan summary ======`: targets, number of steps, and bond changes of each stage |
| `result.json` (`--out-json`) | `scientific_status`: `success` when every step of every stage converged (and the `--endopt` optimization, when requested) with a finite energy; `partial` when only some stages did; `failed` when none did |
| `result.json` (`--out-json`) | `stages[].converged`, `stages[].bond_changes.changed`, `stages[].final_energy_hartree`, and the energy of every step in `stages[].energies_hartree` |

A `partial` run exits with 0 and a `failed` run with 1; for stages that did not converge, see {ref}`max_cycles and plateau stops <troubleshooting-max-cycles>`. A converged scan with the intended bond changes gives a candidate path; the highest-energy step is a TS candidate for [`tsopt`](tsopt.md), which you can {ref}`extract <trajectory-one-frame>` from `scan_trj.xyz`.

---

(scan-direction-and-barrier-sign)=
## Barrier sign

`scan` records energies but does not report a barrier. If you read a barrier off a scan, the forward barrier is always computed from the reactant:

| You ran | Difference from the starting structure | Forward barrier |
| --- | --- | --- |
| A scan from the reactant | `E(TS) − E(reactant)` | The same difference |
| A scan from the product | `E(TS) − E(product)`, the **reverse** barrier | `E(TS) − E(reactant)`, **not** the difference from the starting structure; E(reactant) comes from an optimized reactant, for example an IRC endpoint optimized with [`opt`](opt.md) |

No option changes this. Before quoting a barrier, check which endpoint the scan started from, especially when the starting structure was a crystallographic product complex.

---

## Output files

`scan` writes these files to `--out-dir`:

```text
result_scan/
├─ preopt/
│  └─ result.xyz                    # Pre-optimized structure (--preopt)
├─ stage_01/                        # One directory per stage (stage_NN)
│  ├─ result.xyz                    # Final geometry of the stage
│  ├─ scan_trj.xyz                  # Structure and energy of every step in the stage
│  └─ scan_s0001_optimization_trj.xyz  # Optimizer trajectory of each step (--dump)
├─ scan_trj.xyz                     # All stages joined
└─ result.json                      # Summary (--out-json); summary.json has the same content
```

For PDB/mmCIF input, each `result.xyz` is also written as `result.pdb` and each `scan_trj.xyz` as `scan.pdb`, in the same directory; for Gaussian input, the structures also as `result.gjf`. {ref}`mmCIF input <mmcif-input>`, and PDB input too large for the PDB columns, also get `.cif` files that keep the original identifiers.

* **Stage results**: `stage_NN/result.*` is the structure at the end of stage NN. For [`path-search`](path-search.md), give the starting structure (`preopt/result.*` with `--preopt`) followed by the `stage_NN/result.*` files in stage order, as `all` does.
* **Energy profile**: the comment line of each frame in `scan_trj.xyz` holds the energy without restraints (Hartree); plot it with [`trj2fig`](trj2fig.md).

---

## Main options

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path | (required) | Input structure (`.pdb`, `.cif`, `.mmcif`, `.xyz`, ...) |
| `-q, --charge` | integer | `None` | Total charge. Required unless `-l` is given or the input is `.gjf` |
| `-m, --multiplicity` | integer | `1` | Spin multiplicity (2S+1) |
| `-l, --ligand-charge` | text | `None` | Total ligand charge (for example `-1`) or a charge per residue name (for example `'GPP:-3,SAM:1'`), used when `-q` is omitted (PDB/mmCIF input or `--ref-pdb`) |
| `-s, --scan-lists` | text | (required) | A YAML/JSON spec file, or one or more inline literals (one per stage): distance targets `(i,j,target)`, or ranges for a distance `(i,j,low,high)`, an angle `(i,j,k,low,high)`, or a dihedral `(i,j,k,l,low,high)` |
| `-o, --out-dir` | path | `./result_scan/` | Output directory |
| `--one-based/--zero-based` | flag | `--one-based` | Read atom indices in `-s` as 1-based or 0-based |
| `--max-step-size` | float | `0.2` | Largest change of a distance per step (Å) |
| `--max-angle-step-size` | float | `5.0` | Largest change of an angle per step (degrees) |
| `--max-dihedral-step-size` | float | `10.0` | Largest change of a dihedral per step (degrees) |
| `--restraint-k` | float | `300.0` | Restraint strength k (eV/Å² for distances, eV/rad² for angles); alias `--bias-k` |
| `--preopt/--no-preopt` | flag | `False` | Optimize the input structure without restraints before the scan |
| `--endopt/--no-endopt` | flag | `False` | Optimize the result of each stage without restraints |
| `--dump/--no-dump` | flag | `False` | Write the optimizer trajectory of every step |
| `--opt-mode` | `grad` / `hess` | `grad` | Relaxation: L-BFGS / RFO (on `tsopt` the same words select other optimizers; see {ref}`--opt-mode by command <opt-mode-semantics>`) |
| `--freeze-links/--no-freeze-links` | flag | `True` | Freeze the parent atoms of cap hydrogens at the cluster boundary |
| `--out-json/--no-out-json` | flag | `False` | Write a summary to `result.json` ([JSON Output Reference](json-output.md)) |

See the [generated CLI reference](reference/commands/scan.md) for every option.

---

## Notes

* **`--preopt` depends on the caller**: `scan` alone pre-optimizes only with `--preopt`. Inside `all`, it follows `all --preopt` (on by default), and `all --scan-preopt/--no-scan-preopt` overrides it.
* **Tuples from `all -s`**: write them in the `scan` form ({ref}`Scan-list spec <scan-list-spec>`).
* **Targets and ranges are not mixed inline**: one inline literal, and all literals of one run, hold either targets `(i,j,target)` or ranges. To combine them, list them under `stages:` in a YAML/JSON spec.
* **Stage numbers with ranges**: a range becomes two stages, toward `low` and then toward `high`. Inline, all ranges of one literal share these two stages. In a YAML `stages:` list, a stage of targets only stays one stage; a stage that holds a range is split into one stage per target and two per range.
* **Target distances must be positive**, and one coordinate may appear only once per stage.
* **Check the spec without computing**: `--dry-run` reads the input, the charge and spin, and `-s`, prints the number of stages, and exits without any optimization.
* **Cycle limit**: `--relax-max-cycles` (default `100000`) limits each relaxation; when given, it overrides YAML `opt.max_cycles`.
* **Restraint strength in YAML**: in `--config`, {ref}`bias.k <bias-section>` applies when `--restraint-k` is not given.

---

## See also

* {ref}`Scan-list spec <scan-list-spec>` — YAML/JSON spec files, inline literals, and atom selectors
* [scan2d](scan2d.md) — energy map over two coordinates
* [scan3d](scan3d.md) — energy grid over three coordinates
* [path-search](path-search.md) — minimum energy path (MEP) search from the scan results
* [all](all.md) — the full workflow, including a scan from one structure with `-s`
* [Troubleshooting](troubleshooting.md) — what to do when a run fails
