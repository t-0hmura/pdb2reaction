# Tips for studying reaction mechanisms

The TS and the path that pdb2reaction finds are candidates for the mechanism you propose. This page takes you from a hypothesis to a calculation, a checked TS, a comparison of candidate mechanisms, and a barrier, with the next moves to try when the TS search fails.

## Quick guide

| Goal or symptom | Next move | Section |
| --- | --- | --- |
| Start from a proposed mechanism | List the bonds that form and break and the H atoms that move, then pick the input mode | {ref}`Start from a hypothesis <mechanism-hypothesis>` |
| Set up a concerted or stepwise calculation | Group or split the `-s` literals | {ref}`Decide how to split the reaction <mechanism-split>` |
| See whether you got a TS | Check n_imag and the IRC endpoints | {ref}`Check the TS <mechanism-check-ts>` |
| n_imag ≥ 2, and the extra modes lie outside the reacting site | Add `--flatten` | {ref}`When the TS search fails <ts-search-fails>` |
| n_imag ≥ 2, and two modes both move the reacting bonds | Try the reaction as separate stages | {ref}`When the TS search fails <ts-search-fails>`, {ref}`Decide how to split the reaction <mechanism-split>` |
| n_imag = 0, or the candidate slid toward R or P | Add `--refine-path`, or start from another structure | {ref}`When the TS search fails <ts-search-fails>` |
| Compare candidate mechanisms, or the IRC endpoints are not the intended R and P | Swap the stage order, compare stepwise and concerted runs, revisit the model | {ref}`Compare candidate mechanisms <mechanism-compare>` |
| Read the barrier | Count it from the minimum just before that stage | {ref}`Read the barrier <mechanism-barrier>` |
| The optimization stops at max cycles | Switch the optimizer (`tsopt --opt-mode` / `all --opt-mode-post`), reduce the step size, or start from another structure | {ref}`Troubleshooting: TS optimization <troubleshooting-ts>` |

(mechanism-hypothesis)=
## Start from a hypothesis

Before you run anything, write down what the proposed mechanism does: the bonds that form, the bonds that break, and every H atom that moves, including each H relayed by a water or a base. This list becomes the `-s` coordinates and the points you check at the IRC endpoints.

Pick the [input mode](getting-started.md#choosing-an-input-mode) from the structures you have:

- **R and P (and any intermediates)**: list them in [`all`](quickstart-all.md).
- **R only**: build the path from R with a [scan](quickstart-scan.md).
- **A TS candidate only**: use [TS-only mode](quickstart-tsopt.md).

(mechanism-split)=
## Decide how to split the reaction

Each `-s` literal (one bracketed list after `-s`) is one stage. The coordinates inside one literal move together in the same stage (concerted); several literals run one after another (stepwise), and the restraints of a stage are released when the next stage starts. The number of steps in a stage is set by the coordinate that needs the most steps when its change (target minus start) is divided by the step cap. The step cap is `--scan-max-step-size` for distances (0.20 Å in `all`), 5° for angles, and 10° for dihedrals.

The two commands below move the same four coordinates of the bundled example: the C–C bond that forms, the C–S bond that breaks, and the H that moves from GPP to Glu186, written as two distances (C7–H11 and OE2–H11).

Run all four coordinates as one concerted stage:

```bash
pdb2reaction all -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -s '[("SAM,320,CS1","GPP,321,C7",1.50),("SAM,320,CS1","SAM,320,SD",3.30),("GPP,321,C7","GPP,321,H11",2.90),("GLU,186,OE2","GPP,321,H11",1.00)]' \
    --tsopt --thermo -o ./result_concerted
```

Run the C–C formation with the C–S cleavage first, then the H transfer as a second stage:

```bash
pdb2reaction all -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -s '[("SAM,320,CS1","GPP,321,C7",1.50),("SAM,320,CS1","SAM,320,SD",3.30)]' \
       '[("GPP,321,C7","GPP,321,H11",2.90),("GLU,186,OE2","GPP,321,H11",1.00)]' \
    --tsopt --thermo -o ./result_stepwise
```

The bundled `1.R.pdb` has an empty chain column, so its selectors use the residue name, the residue number, and the atom name (`"SAM,320,CS1"`); with chain IDs, write `"A:SAM:320:CS1"`. All accepted forms are in {ref}`Scan-list spec <scan-list-spec>`.

- **Write `-s` once**: list every literal after a single `-s`. This form works in both `all` and `scan`; `all` rejects a repeated `-s`.
- **Put every moving coordinate in its stage**: include the bond that breaks and each H that moves, not only the bond that forms. Do not drive two coordinates and expect the rest to follow.
- **Without a scan**: when you can prepare R and P (and intermediates, if any), list them in `-i` for an MEP search. Each neighbouring pair becomes one segment; with `--refine-path`, the path is split into segments where bonds change.

(mechanism-check-ts)=
## Check the TS

A successful TS optimization gives one imaginary mode along the reaction coordinate. That gives you a TS candidate, and the IRC confirms it by ending at the intended R and P. An IRC that stops before it converges is still usable when the endpoint optimizations reach the intended R and P.

- **How the run ended**: `tsopt` prints one closing line, such as `[tsopt] Converged (n_imag=1).`; add `--out-json` to also get `result.json`. `all --tsopt` prints the same line for each TS and ends with `Scientific status: success` under `====== Pipeline summary ======` when every requested stage converged; section [3] (`Per-segment post-processing`) of `summary.log` lists n_imag and the IRC outputs of each segment. The verdicts are explained in [tsopt](tsopt.md) and [all → Reading the run status](all.md#reading-the-run-status).
- **Max cycles and plateau stops**: a TS optimization (RS-P-RFO, RS-I-RFO, TRIM, or Dimer) that reaches max cycles without converging does not compute the Hessian, so no n_imag is reported. A run stopped on an {ref}`energy plateau <troubleshooting-max-cycles>` (`--stop-plateau`) always computes the Hessian and reports n_imag.
- **Endpoints**: an exit code of 0 alone does not show that you have the TS you wanted. After the endpoint optimizations, check that the covalent bonds and the positions of the moving H atoms at both IRC ends match the intended R and P.
- **Diagnostic IRC**: when n_imag ≥ 2 but the optimization converged numerically (not a plateau stop) and the final PHVA (partial Hessian vibrational analysis) finished, `all` still runs IRC. It follows the imaginary mode that matches the reaction direction, or the lowest imaginary mode when none can be chosen, and the log says `this is not first-order TS certification`.
- **Modes**: open `vib/imag_*_trj.xyz` in a viewer and see which atoms move in each mode.

(ts-search-fails)=
## When the TS search fails

### Extra imaginary modes

- **Extra modes outside the reacting site** (such as a rotating side chain or water): rerun with `--flatten` (`tsopt`, `opt`, `all`). When n_imag ≥ 2 after the optimization, it displaces the structure along the extra modes and optimizes again, for up to {ref}`50 rounds <flatten-precedence-caveat>` by default. The closing line of `tsopt` suggests this move too.
- **After `--flatten`**: check both IRC ends again. n_imag can reach 1 on a TS candidate of a different reaction.
- **Without `--flatten`**: the mode file of an extra mode, `vib/imag_*_trj.xyz` (also `vib/imag_*.pdb` for PDB input), holds 20 frames, and frames 6 and 16 are the largest displacements in the two directions. Save each of the two as its own file in a viewer and start TS-only mode (`all -i <frame>.pdb -l ... --tsopt`) or `tsopt` from each, with the same charge and spin.
- **Extra modes that also move the reacting bonds**: two stages may overlap in one candidate. Run them as {ref}`separate stages <mechanism-split>`.

The other remedies (precision, coordinate type) are in [`tsopt` → Wrong imaginary-mode count after optimization](tsopt.md#wrong-imaginary-mode-count-after-optimization).

### No imaginary mode, or the candidate slid away

- **Look at the HEI** (highest-energy image): when the forming and the breaking bond are both long at the HEI, or the HEI sits after the bond has broken, the TS optimization tends to lose the reaction mode. Check the bond lengths at the HEI.
- **`--refine-path`**: in `all`, it replaces the single-pass `path-opt` with a recursive `path-search` that refines the path and picks the HEI again. It is off by default; look at the coarse MEP first. The closing line of `tsopt` suggests this move too.
- **Path settings**: the number of MEP images (`--max-nodes`, default 20), the method (`--mep-mode dmf`), and the GSM node placement all move the HEI. `--gsm-param energy` puts more nodes where the energy is high; the default `equi` spaces them evenly.
- **Another starting structure**: restart TS-only mode from the HEI of the unsplit MEP (`hei.xyz` from `path-opt`), the HEI of each segment (`hei_seg_NN.xyz` from `path-search`), or a frame near the top of a scan. In `all`, the HEI files are `hei_seg_NN.xyz` under `_work/path_opt/`, or `_work/path_search/` with `--refine-path`.

(mechanism-compare)=
## Compare candidate mechanisms

- **Stage order**: the stages set the order of the structures passed to the MEP search (R → end of stage 1 → end of stage 2 …). Another order gives another path, so run the proposed order and its reverse, or every plausible order for three or more stages, and compare the energy diagrams.
- **Stepwise or concerted**: run it stepwise and check whether the intermediate bond state survives an unbiased optimization once the restraints are released (`--scan-endopt`). If it relaxes back, the concerted path fits better. If a concerted run splits into two segments with `--refine-path` and passes the same intermediate, the stepwise path fits better.
- **Model**: when the IRC endpoints are not the intended ones, or a residue that takes part (H acceptor, water, partner chain) is not wholly in the model, [enlarge the model](model-setup.md).
- **DFT**: to check a candidate TS with DFT, see [Refine an MLIP TS with DFT](dft-backend.md).

(mechanism-barrier)=
## Read the barrier

- The top of a scan is a TS candidate. Read the barrier from the energy after the TS optimization.
- Count the barrier from the minimum just before that stage: `all` reports each segment's barrier as E(TS) − E(R) of that segment. In TS-only mode, R is the higher-energy end of the IRC. When you compare candidate mechanisms, use the same reference R for all of them.
- After `--tsopt`, each segment's barrier from the optimized TS is in section [3] (`Per-segment post-processing`) of `summary.log` and in `post_segments[].mlip.barrier_kcal` of `summary.json`. In an MEP search, section [2] (`Segment-level MEP summary`) gives the MEP barrier before the TS optimization. Section [4] (`Energy diagrams`) tabulates the energy diagrams. See [all → Reading the run status](all.md#reading-the-run-status).

## Notes

- The diagnostic IRC is a clue to where a mode leads, not a TS check.
- `--flatten` removes extra imaginary modes; it cannot create a reaction mode that is missing.
- On a coarse path, `--refine-path` can split the reaction into unneeded stages, each with its own MEP, TS optimization, IRC, and frequency calculation.
- Changing the cutoff that counts n_imag does not bring a structure closer to a TS; how n_imag is counted is in [tsopt](tsopt.md).

## See also

- [Getting Started](getting-started.md) — choosing an input mode
- [tsopt](tsopt.md) — reading the TS result
- [all](all.md) — the full workflow and its outputs
- [scan](scan.md) — restrained scans and stage literals
- [path-opt](path-opt.md) — the MEP between two structures and its HEI
- [path-search](path-search.md) — recursive MEP search with segments
- [irc](irc.md) — following the reaction mode to R and P
- [freq](freq.md) — counting imaginary modes
- [Building the cluster model](model-setup.md) — checking and enlarging the active-site model
- [Refine an MLIP TS with DFT](dft-backend.md) — DFT checks of a candidate TS
- [Troubleshooting](troubleshooting.md) — errors and convergence problems
- [Common options and selectors](cli-conventions.md) — atom selectors for `-s`
