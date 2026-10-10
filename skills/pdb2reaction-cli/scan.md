# scan, scan2d, scan3d

## When to use

- `scan` drives distances, angles, or dihedrals step by step under harmonic
  restraints and relaxes everything else, in stages that run in sequence. Use
  it to build a candidate path from one structure or to pick frames as later
  path endpoints. When you also want the MEP, TS, and IRC, use `all -s`
  ([all-scan-list.md](all-scan-list.md)); run `scan` alone when the scan
  itself is the task.
- `scan2d` relaxes a grid over two coordinates. It shows how they couple, but
  the grid alone does not prove a concerted or stepwise mechanism.
- `scan3d` relaxes a grid over three coordinates. Use it only when three are
  needed; fewer coordinates need fewer optimizations but answer a different
  question. A grid neither models a change of electronic state nor
  establishes a mechanism.

## Writing -s

`-s` is a Python literal: single quotes outside, double quotes inside. A tuple
names its atoms by 1-based index (`(1, 5, 1.4)`) or by atom spec. An atom spec
has residue name, residue number, and atom name in any order, separated by
spaces, commas, colons, slashes, backticks, or backslashes (`"CS1 SAM 320"`,
`"SAM 320 CS1"`). When residue names or numbers repeat, use the positional
`CHAIN:RESNAME:RESSEQ[ICODE]:ATOM`, for example
`("A:SAM:320:CS1", "B:GPP:321:C7", 1.60)`.

- `(i, j, target)` drives a distance to a target (`scan` only).
- `(i, j, low, high)`, `(i, j, k, low, high)`, and `(i, j, k, l, low, high)`
  are distance, angle, and dihedral ranges, in Å and degrees. In `scan` a
  range becomes two stages, toward `low` and then toward `high`. `scan2d` takes
  exactly two ranges and `scan3d` exactly three, in one literal or under
  `pairs:` in a YAML/JSON file; each range is one grid axis.
- `scan` reads a 4-tuple as a distance range; `all -s` reads it as an angle
  target and takes no ranges. Within one `scan` run, use either targets or
  ranges.

In `scan`, the tuples of one literal move together as one stage
(`'[(a, b, 1.6), (c, d, 3.0)]'` drives two bonds at once). Several literals
after one `-s` run as stages in sequence, each starting from the final
geometry of the previous stage.

## Minimal run

### scan

```bash
pdb2reaction scan -i 1.R.pdb -l 'SAM:1,GPP:-3' \
    -s '[("CS1 SAM 320","C7 GPP 321",1.60)]' \
    -b uma -o result_scan

# Two stages in sequence
pdb2reaction scan -i 1.R.pdb -l 'SAM:1,GPP:-3' \
    -s '[("CS1 SAM 320","C7 GPP 321",1.60)]' \
       '[("H11 GPP 321","OE2 GLU 186",0.90)]' \
    -b uma -o result_scan_staged
```

### scan2d

```bash
pdb2reaction scan2d -i 1.R.pdb -l 'SAM:1,GPP:-3' \
    -s '[("CS1 SAM 320","C7 GPP 321",1.60,3.10), ("H11 GPP 321","OE2 GLU 186",0.90,2.40)]' \
    -b uma -o result_scan2d
```

### scan3d

```bash
pdb2reaction scan3d -i 1.R.pdb -l 'SAM:1,GPP:-3' \
    -s '[("CS1 SAM 320","C7 GPP 321",1.50,3.00), ("H11 GPP 321","OE2 GLU 186",0.90,2.50), ("SD SAM 320","CS1 SAM 320",1.80,3.00)]' \
    -b uma -o result_scan3d
```

To redraw a finished `scan3d` grid without computing energies, give its
`surface.csv` with `--csv`, for example with a new `--zmin`/`--zmax` and a
separate `-o`; `-i`, `-s`, and `-q` are then not needed.

## Judge success

**scan.** Each stage prints `[stage k] Covalent-bond changes (start vs final):`
followed by `Yes` and the formed and broken bonds, or `No`; the run ends with
`====== Scan summary ======`. With `--out-json`, `result.json` separates
`execution_status` (`completed`/`failed`) from `scientific_status`
(`success`/`partial`/`failed`); `partial` exits with 0 and `failed` with 1.
A point that does not converge does not stop its stage: the next point starts
from its geometry, and `stages[i]["converged"]` reports only the last point
(the end optimization when `--endopt/--no-endopt` is on). A stage is usable
only when every point converged: read `usable` in its `stage_outcomes` entry,
and `scientific_status`, which is `success` only when every stage is usable. A
point stopped on an energy plateau (YAML `opt.energy_plateau`) is not
converged. Then check the target, `final_energy_hartree`, and the trajectory.
Inside `all`, read
`_work/scan/result.json`; `all` continues past a `partial` scan whose stages'
last points converged, printing `[all] WARNING: Scan has incomplete
intermediate steps; continuing from valid terminal seeds.`, and stops on a
`failed` one.

```text
result_scan/
├─ preopt/result.xyz       # with --preopt
├─ stage_NN/result.xyz     # final attempted geometry of stage NN
├─ stage_NN/scan_trj.xyz   # steps of stage NN; empty when the target equals the start
├─ stage_NN/scan_*.xyz     # optimizer steps, with --dump
├─ scan_trj.xyz            # all stages joined
└─ scan.pdb                # PDB/mmCIF input or --ref-pdb; scan.cif for mmCIF or oversized PDB
```

**scan2d and scan3d.** A point is usable when its relaxation converged with a
finite energy. In `result.json` (`--out-json`), `scientific_status` is
`success` when every point is usable, `partial` when some are (exit 0), and
`failed` when none is (exit 1); `n_points_attempted` and `n_points_usable`
give the counts, and `grid_points[]` maps each grid index to its values,
energy, convergence, and structure file. `execution_status: completed` only
means the run finished. `surface.csv` holds the per-point energies and
`bias_converged`.

- `grid/point_iDDD_jDDD[_kDDD].{xyz,pdb,cif,gjf}` (DDD = target × 100) is
  the final attempted geometry of each point; when rounded tags collide,
  later names append `_grid_III_JJJ[_KKK]` with zero-based indices.
- `grid/preopt_…` is the starting snapshot (optimized only with `--preopt`);
  `grid/inner_path_d1_NNN[_d2_MMM]_trj.xyz` are written with `--dump`.
- Plots: `scan2d_map.png` and `scan2d_landscape.html`, or
  `scan3d_density.html`. A `--csv` redraw reads the CSV without copying it
  into the output directory.

## Pitfalls and recovery

- Put all stage literals after one `-s`, or repeat `-s` once per stage;
  mixing the two forms is rejected.
- Stage k+1 starts from the final geometry of stage k, so a diverged stage
  spoils every later stage.
- An explicit `--relax-max-cycles` (default 100000 per point;
  `--scan-relax-max-cycles` in `all`) overrides YAML `opt.max_cycles`. Lower
  it for exploratory scans so that points that do not converge do not use up
  the walltime.
- `--max-step-size` (0.20 Å by default) is the largest change of a driven
  distance per scan point and also bounds every relaxation step at that point:
  the L-BFGS `max_step` and the RFO trust radius become the smaller of their
  own value and this size, so a small `--max-step-size` also shortens the
  relaxation steps.
- The three commands have no `--show-config`; use `--dry-run` to check the
  spec without optimizing.
- Grid cost is the product of the axis lengths: a 10 × 10 grid runs 100
  restrained optimizations and a 5 × 5 × 5 grid 125. Use `scan` or `scan2d`
  when fewer coordinates answer the question.
- With too few usable points (under three or all on one line for `scan2d`,
  under four or all in one plane for `scan3d`) only the plots are skipped,
  with `[plot] NOTE: Plots skipped` or `Volume plot skipped` and exit 0. With
  no usable point the run prints `[plot] No finite data for plotting.` and
  exits with 1.
- A `--csv` table needs `d1_A`, `d2_A`, `d3_A`, and `energy_hartree` or
  `energy_kcal`; pre-optimization rows, unconverged rows, and non-finite
  energies are left out.

## Next step

- Plot `scan_trj.xyz`: [`trj2fig`](utilities.md#trj2fig).
- Optimize a high-energy frame or grid point as a TS candidate:
  [tsopt.md](tsopt.md).
- Staged scans inside the full workflow, which avoid a grid when the stages
  are decoupled: [all-scan-list.md](all-scan-list.md).
- Flags and defaults: `--help-advanced` and [SKILL.md](SKILL.md#where-flags-and-defaults-live).
