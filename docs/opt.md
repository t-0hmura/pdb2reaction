# `opt` (geometry optimization)

`opt` optimizes one structure to a local minimum.

---

## What it is for

* **Preparing R, P, and intermediates**: relax the reactant, product, and intermediate structures before a path search or a frequency calculation, and confirm each minimum (n_imag = 0) with [`freq`](freq.md).
* **Relaxing with fixed distances**: keep chosen atom pairs at a set distance while everything else relaxes.
* **Turning IRC endpoints into R and P**: optimize the endpoints of an [`irc`](irc.md) run to the minima they lead to.

The default backend is **UMA** (Meta); `-b/--backend` also selects **ORB**, **MACE**, **AIMNet2**, or **DFT**.

---

## Examples

### 1. Basic minimization

Give the charge and the spin multiplicity explicitly, and write a summary with `--out-json`.

```bash
pdb2reaction opt -i input.pdb -q 0 -m 1 --out-json --out-dir ./result_opt
```

The run converged when the console prints `[opt] Converged!` and `result_opt/result.json` has `"optimization_status": "converged"`.

### 2. Tighter threshold with trajectory

Use the `gau_tight` criteria and keep the optimization trajectory.

```bash
pdb2reaction opt -i input.pdb -q 0 -m 1 --thresh gau_tight --dump \
    --out-dir ./result_opt_tight
```

### 3. Distance restraint

Pull atoms 1 and 5 toward 2.0 Å with a weak harmonic restraint (20 eV·Å⁻²).

```bash
pdb2reaction opt -i input.pdb -q 0 -m 1 \
    --distance-restraint '[(1,5,2.0)]' --restraint-k 20.0 --out-dir ./result_opt_rest
```

### 4. RFO

Switch to RFO, which starts from an exact Hessian, with `--opt-mode hess`.

```bash
pdb2reaction opt -i input.pdb -q 0 -m 1 --opt-mode hess --out-dir ./result_opt_hess
```

---

## How it works

1. **Reading the structure and freezing the boundary**: the {ref}`charge <charge-specification>` comes from `-q` or `-l`. With `--freeze-links` (on by default), the parent atoms of the {ref}`cap hydrogens <link-hydrogen-and-frozen-atoms>` of a cut-out cluster are frozen; `--freeze-atoms` freezes more atoms.
2. **Choosing the optimizer** (`--opt-mode`): `grad` (alias `lbfgs`) runs **L-BFGS**, which uses gradients only. `hess` (alias `rfo`) runs **RFO**, which starts from an exact Hessian, updates it with [TS-BFGS](glossary.md#optimization-algorithms) (YAML `rfo.hessian_update`), and recomputes it every 500 cycles. In `tsopt`, the same tokens {ref}`select other methods <opt-mode-semantics>`.
3. **Adding distance restraints** (`--distance-restraint`): each `(i, j, target)` adds a harmonic term with force constant `--restraint-k` (eV·Å⁻²) that pulls atoms i and j toward `target` in Å; `(i, j)` keeps their starting distance. Indices are 1-based unless `--zero-based` is given.
4. **Minimizing**: the optimizer moves the structure until the convergence criteria are met or `--max-cycles` is reached. The default `--thresh gau` asks for a max force below 4.5 × 10⁻⁴ and an RMS force below 3.0 × 10⁻⁴ hartree/bohr, and a max step below 1.8 × 10⁻³ and an RMS step below 1.2 × 10⁻³ bohr, the same as Gaussian's default.
5. **Removing imaginary modes (only with `--flatten`)**: after the optimization, `opt` computes the Hessian, displaces the structure by 0.10 Å along every imaginary mode (ν < −5.00 cm⁻¹), and optimizes again, for up to 50 rounds or until no imaginary mode is left. With `--flatten`, the console prints n_imag in the line `[Imaginary modes] n=…` after each round, and `[flatten] WARNING: Remaining imaginary modes after the flatten loop: N` when modes are left after the last round.

---

## Checking convergence

How the run ended is printed on the console and recorded in `result.json` (`--out-json`):

| How it ended | `optimization_status` | Console line | `scientific_status` / exit code |
| --- | --- | --- | --- |
| Converged | `converged` | `[opt] Converged!` | `success` / 0 |
| Reached `--max-cycles` without converging | `not_converged` | `[opt] Reached max cycles (N/M).` | `failed` / 1 |
| Stopped on an energy plateau (`--stop-plateau`) | `stalled` | `[opt] Stalled (energy plateau; not converged)` | `failed` / 1 |

Each of these lines is followed by `[opt] Total cycles: N`. For what to change after `not_converged` or `stalled`, see {ref}`max_cycles and plateau stops <troubleshooting-max-cycles>`.

Convergence gives a stationary point, not necessarily a minimum. `opt` computes no final Hessian unless `--flatten` is on, so run [`freq`](freq.md) on the final geometry and check that n_imag = 0.

---

## Output files

When the run finishes, `--out-dir` contains:

```text
result_opt/
├─ final_geometry.xyz      # Final geometry (always written)
├─ final_geometry.pdb      # Same, for PDB/mmCIF input (.gjf for Gaussian input)
├─ optimization_trj.xyz    # Optimization trajectory (--dump)
├─ optimization.pdb        # Same trajectory as PDB (--dump, PDB/mmCIF input)
├─ restart_NNN.yaml        # Optimizer state (--dump with YAML opt.dump_restart)
└─ result.json             # Summary (--out-json)
```

{ref}`mmCIF input <mmcif-input>`, and PDB input too large for the PDB columns, also get `.cif` files that keep the original identifiers.

* **Final geometry**: `final_geometry.*` is the optimized structure to pass to [`freq`](freq.md) or to a path search.
* **Summary**: with `--out-json`, [`result.json`](json-output.md) records `optimization_status`, the final energy `energy_hartree` (without the restraint energy), and the number of cycles `n_opt_cycles`.
* **Console**: the cycle table and the elapsed time; with `-v 3`, also the `geom`, `calc`, `opt`, and `lbfgs` / `rfo` settings actually used.

---

## Main options

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path | (required) | One structure (`.pdb`, `.cif`, `.mmcif`, `.xyz`, `.gjf`). For a trajectory, extract one frame to `.xyz` first (see {ref}`Extract one frame from a trajectory <trajectory-one-frame>`) |
| `-q, --charge` | integer | `None` | Total charge. Required unless `-l` is given or the input is `.gjf` |
| `-l, --ligand-charge` | text | `None` | Total ligand charge (for example `-1`) or a charge per residue name (for example `'GPP:-3,SAM:1'`), used when `-q` is omitted (PDB/mmCIF input or `--ref-pdb`) |
| `-m, --multiplicity` | integer | `1` | Spin multiplicity (2S+1) |
| `--ref-pdb` | path | `None` | PDB/mmCIF topology for an `.xyz` / `.gjf` input; the coordinates come from `-i` (for example IRC endpoints, see [irc](irc.md)) |
| `-b, --backend` | text | `uma` | Backend (`uma`, `orb`, `mace`, `aimnet2`, `dft`) |
| `--opt-mode` | `grad` / `hess` | `grad` | Optimizer: L-BFGS / RFO (`lbfgs` and `rfo` are aliases) |
| `--coord-type` | `cart` / `redund` / `dlc` / `tric` | `cart` | Optimization coordinates: Cartesian / redundant internal / delocalized internal (DLC) / translation-rotation internal (TRIC) |
| `--thresh` | preset | `gau` | Convergence criteria (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`) |
| `--max-cycles` | integer | `100000` | Maximum number of optimization cycles, shared with the `--flatten` rounds |
| `--dump/--no-dump` | flag | `False` | Write the optimization trajectory `optimization_trj.xyz` |
| `--distance-restraint` | text | `None` | Harmonic distance restraints, inline (`'[(i,j,target_Å),...]'`) or as a YAML/JSON file that lists the same entries under `constraints:`; `(i,j)` keeps the starting distance. Atoms can also be given as [atom selectors](cli-conventions.md#atom-selectors) such as `'SAM,320,CS1'` |
| `--restraint-k` | float | `300` | Force constant of the distance restraints (eV·Å⁻²) |
| `--one-based/--zero-based` | flag | `--one-based` | Count `--distance-restraint` indices from 1 or from 0 |
| `--freeze-links/--no-freeze-links` | flag | `True` | Freeze the parent atoms of cap hydrogens (PDB/mmCIF input or `--ref-pdb`) |
| `--freeze-atoms` | text | `None` | Atoms to freeze (1-based, comma-separated, for example `'1,3,5'`) |
| `--flatten/--no-flatten` | flag | `False` | Remove imaginary modes after the optimization |
| `--reject-uphill/--no-reject-uphill` | flag | `False` | With `hess`, reject RFO steps that raise the energy by more than 1e-4 hartree and shrink the trust radius |
| `--stop-plateau/--no-stop-plateau` | flag | `False` | Stop when the energy stops changing (range below 1e-4 hartree over 50 cycles) and report `stalled` |
| `-o, --out-dir` | path | `./result_opt/` | Output directory |

See the [generated CLI reference](reference/commands/opt.md) for every option.

> **Note:** in YAML (`--config`), `geom.freeze_atoms` adds frozen atoms (1-based), merged with `--freeze-links` and `--freeze-atoms`. Every key is listed under [`geom`](yaml-reference.md#geom), [`opt`](yaml-reference.md#opt), [`lbfgs`](yaml-reference.md#lbfgs), and [`rfo`](yaml-reference.md#rfo) in the YAML Reference.

---

## Notes

* **Plateau stop**: `--stop-plateau` saves cycles when force noise keeps the force criteria out of reach, but a flat energy is no evidence of a stationary point. `--max-cycles` remains the real limit. `--stop-plateau-thresh` and `--stop-plateau-window` set the energy range and the number of cycles.
* **Rigid motions with frozen atoms**: the RFO curvature checks in Cartesian coordinates and `--flatten` treat rigid motions as [`freq`](freq.md#rigid-modes-with-frozen-boundaries) does. L-BFGS is not affected.
* **Frozen atoms and restraints in general**: how to choose frozen atoms and restraints for a cluster model is described in {ref}`Frozen atoms and distance restraints <freeze-atoms-and-restraints>`.
* **Optimizer state dumps**: with `--dump`, set YAML `opt.dump_restart` to a positive integer N to write `restart_NNN.yaml` every N cycles. pdb2reaction does not read this file back, so rerun `opt` from the final geometry to continue a stopped calculation.
* **Model and precision**: `--backend-model` selects the model of the backend and `--precision` (`fp32` / `fp64`) its precision; see the generated reference.

---

## See also

* [freq](freq.md) — check that the optimized structure is a minimum (n_imag = 0)
* [tsopt](tsopt.md) — optimize a TS (saddle point) instead of a minimum
* [irc](irc.md) — trace the reaction path from a TS to the endpoints to optimize
* [extract](extract.md) — cut out the active-site model before optimizing
* [all](all.md) — the full workflow, which also optimizes the IRC endpoints
* [Troubleshooting](troubleshooting.md) — when a run fails
* [YAML Reference](yaml-reference.md) — every `opt`, `lbfgs`, and `rfo` setting
* [Glossary](glossary.md) — L-BFGS, RFO, and other terms
* {ref}`Exit codes <exit-codes>` — what each exit status means
