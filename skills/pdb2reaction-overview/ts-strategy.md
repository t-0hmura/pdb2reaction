# TS strategy

## When to read

Read this page when you study a reaction mechanism with pdb2reaction: from a hypothesis to how the reaction is split into stages, a checked TS, what to try when the TS does not come out, multistep paths, comparisons between candidates, and the barrier. Before any run, write down the bonds that form, the bonds that break, and every H atom that moves; this list becomes the `-s` coordinates and the points you check at the IRC endpoints.

## Precision (fp32 or fp64)

- `--precision fp32|fp64` is accepted by UMA, ORB, and MACE; AIMNet2 accepts fp32 only. Left unset, each backend keeps its default: `uma` fp32, `orb` fp64, `mace` fp64, `aimnet2` fp32.
- For UMA, fp64 can change TS optimization and Hessian behaviour, and its speed cost depends on the hardware and model; before production, compare both settings on your system and time one Hessian in each precision on the production GPU.
- fp64 is slow on consumer GPUs, where UMA in fp32 gives the best balance; on HPC GPUs, ORB in fp64 is cost-effective. ORB and MACE in fp32 leave extra imaginary modes more often.
- Precision is not proof of a saddle; check n_imag and the IRC endpoints with either setting.
- `--deterministic` requests run-to-run repeatability on the same software and hardware; it does not increase
  numerical precision, remove physical imaginary modes, or guarantee identical
  output across PyTorch/backend versions, hardware, or custom calculators.

## Two routes to a TS candidate

- **MEP between R and P** (Endpoint mode, or `path-opt`/`path-search`): use it when you have both ends. The highest-energy image (HEI) is the TS candidate; `path-search` splits a multi-step path where bonds change and writes `hei_seg_NN.*` per segment.
- **Restrained scan from R** (`all -s`, or `scan`): use it when you have only R, or want to drive a chosen distance. Each stage holds the listed distances with a harmonic restraint (k = 300 by default) and optimizes everything else, so scan frames are biased structures, not stationary points. In `scan`, `--preopt` optimizes the start structure without restraints and `--endopt` the end of each stage; both are off by default.
- Pick the scan atoms by measuring the structure, not from atom numbers in a paper: a forming bond is long in R and close to a bond length in P. Before scanning a proton transfer, measure the donor–acceptor heavy-atom distance in R; if it is well beyond a hydrogen-bond contact (roughly 3.2 Å) and no water or residue in the model bridges it, revise the hypothesis instead of driving it. If the TS or product puts the moving H on an atom you did not name, the coordinate was wrong; rebuild it rather than accepting the new acceptor.
- Scan relaxations and endpoint pre-optimization use `gau` unless you set `--thresh`, while TS and post-IRC endpoint optimizations in `all` use `--thresh-post baker`. When the scan seeds the TS candidate or its profile looks wrong (for example a barrier near zero), rerun it with `--thresh baker` and a smaller step (`--scan-max-step-size` in `all`, `--max-step-size` in `scan`) before changing the model.
- `opt` does not scan, but it can hold distances with `--distance-restraint` and `--restraint-k`.
- Send either candidate to `tsopt` → `irc`. The final partial Hessian vibrational analysis (PHVA) of `tsopt` counts n_imag; add `freq` for all modes or thermochemistry.

## Wrong n_imag after tsopt

A successful TS optimization gives one imaginary mode along the reaction coordinate. When n_imag ≥ 2, or the one mode does not move the reacting atoms, look first at every mode (`vib/imag_*_trj.xyz`), the starting structure, and why the optimizer stopped. Then:

- **A small extra mode in a model with a frozen boundary** (a few to about 20 cm⁻¹): first re-optimize with a tighter preset (`tsopt --thresh gau_tight` or `gau_vtight`, or `all --thresh-post`) and count n_imag again. If the small mode vanishes or changes sign, it was residual curvature; only a persistent one calls for `--flatten`.
- **Extra modes outside the reacting site**: rerun with `--flatten` (`tsopt`, `opt`, `all`; off by default). It displaces the structure along the extra modes and optimizes again, up to 50 rounds.
- **n_imag still ≥ 2 after `--flatten`**: read `flatten_skip_reason` in the tsopt `result.json`. `target mode is not negative` (or `target mode sign never determined`) means that the reaction mode followed along a reference direction (the MEP tangent in `all`, or `--ref-mode`) is no longer imaginary, so flattening was refused to avoid making a saddle of another reaction: start again from the HEI or another candidate instead of retrying `--flatten`. `max-cycles budget exhausted before flattening` (or `during flattening`) means the flatten rounds share `--max-cycles`.
- **Precision**: compare fp32 and fp64 where the backend allows it (not AIMNet2). Neither removes a real second direction of negative curvature.
- **Coordinates**: `--coord-type` takes `cart` (default), `redund`, `dlc`, or `tric`; `path-opt`, `path-search`, and `all` take only `cart` or `dlc`. Delocalized internal coordinates can change the conditioning; try both on the problem structure, since neither is faster or more reliable everywhere.
- **Poor HEI**: rerun the parent `all` command with `--refine-path` so the recursive `path-search` refines the MEP before TS optimization. This is deliberately off by default: a poor or noisy path can be split into unneeded segments, each with its own MEP, TS optimization, IRC, and frequency cost.
- `--ref-mode` is not a normal standalone remedy. It takes atom-order-matched Cartesian 3N vectors (`.npz`, `.npy`, or text) that guide which Hessian mode is followed; it does not replace the Hessian, and Dimer does not use it. `all` already supplies the MEP tangent to Hessian-based TS optimizers (`--no-tsopt-from-mep-tan` turns this off). Pass it by hand only for a path from outside.

## Product-start scans give the reverse barrier

A scan or path that starts from P reports E(TS) − E(P), the reverse barrier. The forward barrier is E(TS) − E(R). Check which end the run started from before quoting a barrier. In TS-only mode, R is the higher-energy IRC end ([outputs.md](outputs.md#rtsp-paths)).

Read the barrier from the energy after TS optimization, not from the top of the scan, and count it from the minimum just before that stage.

## Multistep paths

For each TS of a multistep path, report the local barrier (counted as above) and the TS height above one common reference: the lowest optimized R of the same model and protonation state, which after IRC and endpoint optimization is often several kcal/mol below the structure you started from. Do not add local barriers, do not count a drop that comes only from switching the reference R, and do not join steps from different models, snapshots, or protonation states into one profile.

When only one step gives a validated TS, or the optimized IRC ends show only part of the intended change (for example the heavy-atom transfer without the proton transfer), treat that optimized end as an intermediate, not a failure. Start the next step from it rather than from a scan stage end: a new Scan-list run from `segments/seg_NN/product.*` (or `reactant.*`) without `-c`, with the same charge and spin, driving only the remaining coordinates. Before joining steps computed separately, compare the optimized P of one step with the optimized R of the next: covalent bonds, the owner of every H, each residue's protonation, then all-atom RMSD and energy. Identical bonding can still differ by a few tenths of an Å and several kcal/mol; connect that conformational gap with a path between the two minima, or report it as unverified. Then draw all steps on one diagram from the common R with `energy-diagram` ([outputs.md](outputs.md#energy-diagrams)).

## Staged vs concerted scans

Each literal after `-s` is one stage. The tuples inside one literal move together in that stage (concerted); several literals run one after another, each in its own `stage_NN/` (staged). Write `-s` once and list every literal after it: `all` rejects a repeated `-s`, and `scan` accepts repeats but rejects mixing the two forms.

```bash
# concerted: one stage, both distances move together
pdb2reaction scan -i R.pdb -l 'SAM:1,GPP:-3' -s '[("SAM,320,CS1","GPP,321,C7",1.60),("GLU,186,OE2","GPP,321,H11",1.00)]' -o scan_concerted
# staged: stage_01, then stage_02
pdb2reaction scan -i R.pdb -l 'SAM:1,GPP:-3' \
    -s '[("SAM,320,CS1","GPP,321,C7",1.60)]' \
       '[("GLU,186,OE2","GPP,321,H11",1.00)]' -o scan_staged
```

- **Choose** from the mechanism you propose: staging a truly concerted event can create an artificial intermediate, and coupling unrelated coordinates can hide a stepwise route. With both endpoints and an unknown mechanism, `path-search` proposes bond-change segments; they are candidates to check with TS optimization, frequencies, and IRC.
- **Hold what should not move yet**: a stage restrains only the coordinates it lists, so an X–H that should react later can move in an earlier stage. Add that distance to each earlier stage with its target set to the value measured on the structure the stage starts from; a tuple whose target equals its start does not move, so its restraint holds the distance until the next stage starts. For the first stage, `all` pre-optimizes the input inside the scan, so measure on that result (`_work/scan/preopt/result.*`, or your own `opt` result) and start `all --no-scan-preopt` from it. The hold is harmonic: check the held distance in `scan_trj.xyz` and again at the MEP ends.
- **A stage end is not an intermediate**: it is a restrained structure, and the MEP pre-optimizes each stage end without restraints (`--preopt`, on by default), which often breaks a bond the stage formed or re-forms one it broke. Compare the driven distances in `_work/scan/stage_NN/result.*` with the same distances at the MEP ends, or rerun with `--scan-endopt` (off by default) and check whether the intermediate survives once the restraints are released; if it relaxes back, the concerted path fits better.
- **Count the steps from the TSs**, not from the `-s` stages: one literal with several tuples often refines into separate TSs, and two stages can collapse into one segment. For each TS, check which bonds and H atoms its mode moves and where they cross in the IRC frames. A proton still on its donor at the TS and at the raw IRC end (`segments/seg_NN/structures/*_irc.*`) that moves only during endpoint optimization completes downhill after that TS; it shows neither a concerted event nor a separate step.
- In `scan`, a range `(i, j, low, high)` becomes two stages, one toward each end.

## When the TS does not come out

- **Judge**: you have a TS when n_imag = 1, the mode moves the reacting atoms, and the optimized IRC ends match the intended R and P in covalent bonds and in the positions of the moving H atoms. An exit code of 0 is not that evidence.
- **Extra modes outside the reacting site**: `--flatten` (see *Wrong n_imag after tsopt*). After `--flatten`, check both IRC ends again; n_imag can reach 1 on a TS candidate of a different reaction.
- **Start from an extra mode**: `vib/imag_*_trj.xyz` holds 20 frames, and frames 6 and 16 are the largest displacements in the two directions. Save each as its own file and start TS-only mode or `tsopt` from each, with the same charge and spin.
- **No imaginary mode, the candidate slid toward R or P, or the step's MEP barrier is low**: the HEI may not yet be a TS. Check its imaginary mode and the forming and breaking bond lengths, then read the MEP energy profile. Two maxima mean the pair spans two steps: give the intermediate or use `--refine-path` rather than raising `--max-nodes`. For a single HEI that does not look like a TS, rerun with `--refine-path`, raise `--max-nodes` (20 by default), or try `--mep-mode dmf` or `--gsm-param energy` (more GSM nodes where the energy is high, when the equidistant path skips the region near the HEI); or run a fine scan between the two ends of that step and take a TS-like frame. Other starting structures are the HEI of the whole MEP (`hei.xyz` from `path-opt`), a segment's `hei_seg_NN.xyz` (`_work/path_opt/`, or `_work/path_search/` with `--refine-path`), or a scan frame near the top. Once a TS has n_imag = 1 and its IRC ends are the intended R and P, the node count of the MEP that seeded it does not change the result; in a batch of snapshots, allow one such path retry per snapshot, then move to the next snapshot.
- **Two modes that both move the reacting bonds, or endpoints that are not the intended R and P**: swap the stage order and compare, compare staged and concerted runs, put every moving coordinate (the breaking bond, each moving H) in its stage, or prepare R and P and run an endpoint MEP.
- **Model**: when a residue that takes part (H acceptor, water, partner chain) is missing or the charge is doubtful, revisit the model with [Check the boundary and the charge](../pdb2reaction-model-setup/SKILL.md#check-the-boundary-and-the-charge).

More moves and examples: [docs/mechanism-tips.md](../../docs/mechanism-tips.md).

## Controlled comparisons

Do not subtract the raw total energies of WT and mutant models. Compute each barrier separately and compare them:

`ΔΔG‡ = (G_TS − G_R)_mutant − (G_TS − G_R)_WT`

Use the same backend, model, precision, and thermochemistry settings for both. How to keep the atoms, boundary, protonation, and caps consistent within one path and across variants: [pdb2reaction-model-setup](../pdb2reaction-model-setup/SKILL.md#same-atoms-across-states-and-variants).

- **Mechanism candidates**: cut one cluster and run every hypothesis on it, so that all barriers share the same R. Build each candidate explicitly (its own `-s` stages or endpoints, concerted and stepwise); an MEP or TS search that ends on the saddle of one mechanism, or a model with a reacting fragment removed, is not evidence against another. Rank the candidates by barriers whose TS has n_imag = 1 and the intended IRC ends, after `--tsopt --thermo` on every candidate that could be the lowest, not by MEP-level barriers (`segments[].barrier_kcal`, or `rate_limiting_step` with `method` `MEP`): the verified barrier can differ by many kcal/mol and reverse the order.
- **MD snapshots**: one structure gives one static barrier. Run the same protocol and flags on several snapshots that sample the reactive arrangement ([building their models](../pdb2reaction-model-setup/SKILL.md#build-the-cluster)), and compare barriers, not total energies. Before counting, deduplicate TSs by energy and geometry: several scan orders or settings from one snapshot often reach the same TS, and two inputs built from the same endpoints are one candidate. Report the TSs, the snapshots with a validated TS, and the fully connected paths as separate counts, with the lowest and median barrier rather than a single value.
- **Repeat runs**: runs are not bitwise repeatable by default. On the same input and GPU model, optimizer paths can separate within a few tens of cycles and end at different structures, or one TS run can converge while the other stops. Do not attribute a difference between two single runs to the setting you changed, and do not discard a candidate after one failed run; repeat the run, or compare with `--deterministic` on the same hardware and software.

## Next step

- [outputs.md](outputs.md): where the barriers, n_imag, and R/TS/P structures are.
- [pdb2reaction-cli](../pdb2reaction-cli/SKILL.md): `scan`, `path-opt`/`path-search`, `tsopt`, `irc`, and `freq`.
- Docs: [mechanism-tips](../../docs/mechanism-tips.md), [tsopt](../../docs/tsopt.md#wrong-imaginary-mode-count-after-optimization), [scan](../../docs/scan.md), [backends](../../docs/backends.md) (precision and determinism).
