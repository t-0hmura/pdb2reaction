# TS strategy

## When to read

Read this page when you study a reaction mechanism with pdb2reaction: from a hypothesis to how the reaction is split into stages, a checked TS, what to try when the TS does not come out, comparisons between candidates, and the barrier. Before any run, write down the bonds that form, the bonds that break, and every H atom that moves; this list becomes the `-s` coordinates and the points you check at the IRC endpoints.

## Precision (fp32 or fp64)

- `--precision fp32|fp64` works on every MLIP backend. Left unset, each backend keeps its default: `uma` fp32, `orb` fp64, `mace` fp64, `aimnet2` fp32.
- AIMNet2 accepts fp32 only.
- For UMA, fp64 can change TS optimization and Hessian behaviour, and its speed cost depends on the hardware and model; compare both settings on your system before production.
- Precision is not proof of a saddle; check n_imag and the IRC endpoints with either setting.
- `--deterministic` requests run-to-run repeatability on the same software and hardware; it does not increase
  numerical precision, remove physical imaginary modes, or guarantee identical
  output across PyTorch/backend versions, hardware, or custom calculators.

## Two routes to a TS candidate

- **MEP between R and P** (endpoint MEP mode, or `path-opt`/`path-search`): use it when you have both ends. The highest-energy image (HEI) is the TS candidate; `path-search` splits a multi-step path where bonds change and writes `hei_seg_NN.*` per segment.
- **Restrained scan from R** (`all -s`, or `scan`): use it when you have only R, or want to drive a chosen distance. Each stage holds the listed distances with a harmonic restraint (k = 300 by default) and optimizes everything else, so scan frames are biased structures, not stationary points. In `scan`, `--preopt` optimizes the start structure without restraints and `--endopt` the end of each stage; both are off by default.
- `opt` does not scan, but it can hold distances with `--distance-restraint` and `--restraint-k`.
- Send either candidate to `tsopt` → `irc`. The final partial Hessian vibrational analysis (PHVA) of `tsopt` counts n_imag; add `freq` for all modes or thermochemistry.

## Wrong n_imag after tsopt

A successful TS optimization gives one imaginary mode along the reaction coordinate. When n_imag ≥ 2, or the one mode does not move the reacting atoms, look first at every mode (`vib/imag_*_trj.xyz`), the starting structure, and why the optimizer stopped. Then:

- **Extra modes outside the reacting site**: rerun with `--flatten` (`tsopt`, `opt`, `all`; off by default). It displaces the structure along the extra modes and optimizes again, up to 50 rounds.
- **Precision**: compare fp32 and fp64 where the backend allows it (not AIMNet2). Neither removes a real second direction of negative curvature.
- **Coordinates**: `--coord-type` takes `cart` (default), `redund`, `dlc`, or `tric`; `path-opt`, `path-search`, and `all` take only `cart` or `dlc`. Delocalized internal coordinates can change the conditioning; try both on the problem structure, since neither is faster or more reliable everywhere.
- **Poor HEI**: rerun the parent `all` command with `--refine-path` so the recursive `path-search` refines the MEP before TS optimization. This is deliberately off by default: a poor or noisy path can be split into unneeded segments, each with its own MEP, TS optimization, IRC, and frequency cost.
- `--ref-mode` is not a normal standalone remedy. It takes atom-order-matched Cartesian 3N vectors (`.npz`, `.npy`, or text) that guide which Hessian mode is followed; it does not replace the Hessian, and Dimer does not use it. `all` already supplies the MEP tangent to Hessian-based TS optimizers (`--no-tsopt-from-mep-tan` turns this off). Pass it by hand only for a path from outside.

## Product-start scans give the reverse barrier

A scan or path that starts from P reports E(TS) − E(P), the reverse barrier. The forward barrier is E(TS) − E(R). Check which end the run started from before quoting a barrier. In TS-only mode, R is the higher-energy IRC end ([outputs.md](outputs.md#rtsp-paths)).

Read the barrier from the energy after TS optimization, not from the top of the scan, and count it from the minimum just before that stage.

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

- Choose from the mechanism you propose: staging a truly concerted event can create an artificial intermediate, and coupling unrelated coordinates can hide a stepwise route.
- With both endpoints and an unknown mechanism, `path-search` proposes bond-change segments; they are candidates to check with TS optimization, frequencies, and IRC.
- In `scan`, a range `(i, j, low, high)` becomes two stages, one toward each end.

## When the TS does not come out

- **Judge**: you have a TS when n_imag = 1, the mode moves the reacting atoms, and the optimized IRC ends match the intended R and P in covalent bonds and in the positions of the moving H atoms. An exit code of 0 is not that evidence.
- **Extra modes outside the reacting site**: `--flatten` (see *Wrong n_imag after tsopt*). After `--flatten`, check both IRC ends again; n_imag can reach 1 on a TS candidate of a different reaction.
- **Start from an extra mode**: `vib/imag_*_trj.xyz` holds 20 frames, and frames 6 and 16 are the largest displacements in the two directions. Save each as its own file and start TS-only mode or `tsopt` from each, with the same charge and spin.
- **No imaginary mode, or the candidate slid toward R or P**: check the bond lengths at the HEI, try `--refine-path`, or start from another structure: the HEI of the whole MEP (`hei.xyz` from `path-opt`), a segment's `hei_seg_NN.xyz` (`_work/path_opt/`, or `_work/path_search/` with `--refine-path`), or a scan frame near the top.
- **Two modes that both move the reacting bonds, or endpoints that are not the intended R and P**: swap the stage order and compare, compare staged and concerted runs, put every moving coordinate (the breaking bond, each moving H) in its stage, or prepare R and P and run an endpoint MEP.
- **Model**: when a residue that takes part (H acceptor, water, partner chain) is missing or the charge is doubtful, revisit the model with [pdb2reaction-model-setup](../pdb2reaction-model-setup/SKILL.md).

More moves and examples: [docs/mechanism-tips.md](../../docs/mechanism-tips.md).

## Controlled comparisons

Do not subtract the raw total energies of WT and mutant models. Compute each barrier separately and compare them:

`ΔΔG‡ = (G_TS − G_R)_mutant − (G_TS − G_R)_WT`

Use the same backend, model, precision, and thermochemistry settings for both. How to keep the atoms, boundary, protonation, and caps consistent within one path and across variants: [pdb2reaction-model-setup](../pdb2reaction-model-setup/SKILL.md#same-atoms-across-states-and-variants).

## Next step

- [outputs.md](outputs.md): where the barriers, n_imag, and R/TS/P structures are.
- [pdb2reaction-cli](../pdb2reaction-cli/SKILL.md): `scan`, `path-opt`/`path-search`, `tsopt`, `irc`, and `freq`.
- Docs: [mechanism-tips](../../docs/mechanism-tips.md), [tsopt](../../docs/tsopt.md#wrong-imaginary-mode-count-after-optimization), [scan](../../docs/scan.md), [backends](../../docs/backends.md) (precision and determinism).
