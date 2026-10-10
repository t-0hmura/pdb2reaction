# `pdb2reaction all`: Multi-structure MEP search

Give two or more structures in reaction order; `all` finds the MEP between each
neighbouring pair and, with `--tsopt`, optimizes each TS candidate and runs
IRC. It succeeded when the console prints `[tsopt] Converged (n_imag=1).` for
each TS and `Scientific status: success` under the last
`====== Pipeline summary ======`.

## When to use

You have two or more structures in reaction order (reactant, optional
intermediates, product) with the same atoms in the same order, typically R and
P (sometimes IM) from a published QM or QM/MM study. By default `all` runs one
single-pass `path-opt` per neighbouring pair and does not look for hidden
intermediates. With `--refine-path`, the recursive `path-search` runs once over
the whole ordered series and splits it where bonds change, so you need to give
only the steps you already know.

## Minimal run

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo -o result_mep
```

Add `--dft` (and `--func-basis 'wb97m-v/def2-tzvpd'`) for DFT single points on
R, TS, and P. For a known multistep mechanism, give each intermediate:

```bash
pdb2reaction all -i 1.R.pdb 2.IM.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo -o result_mep_3pt
```

## Same atoms in the same order

Every input needs the same number of atoms, the same element sequence, and the
same residue assignments. `extract` checks that several inputs match, but it
does not map or repair a mismatched series. If the inputs came from different
programs or were renumbered, first reorder them to one common atom order, then
cut all of them in one extraction so they share the same cluster boundary:

```bash
pdb2reaction extract -i 1.R_raw.pdb -i 3.P_raw.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -o 1.R.pdb -o 3.P.pdb
```

For states and variants, see
[Same atoms across states and variants](../pdb2reaction-model-setup/SKILL.md#same-atoms-across-states-and-variants).

## Judge success

Read the console, `summary.json`, and the endpoints as in
[all.md](all.md#judge-success). For this mode also check:

- **MEP**: open `mep_trj.xyz` and `energy_diagram_MEP.png`; the TS candidates are `_work/path_opt/hei_seg_NN.*` (`_work/path_search/` with `--refine-path`).
- **Segments**: `summary.json["segments"]` lists `index`, `kind`, `barrier_kcal`, `delta_kcal`, and `bond_changes` for each segment. The default gives one segment per neighbouring pair. With `--refine-path`, compare `n_segments_reactive` with the number of inputs minus one; `n_segments` also counts `"kink"` and `"bridge"` segments.
- **R/TS/P**: with `--tsopt`, `segments/seg_NN/{reactant,ts,product}.*` are written after IRC and the endpoint optimizations (RFO by default, L-BFGS with `--opt-mode-post grad`).

## Pitfalls and recovery

- **More bond changes than the inputs imply.** The reaction encoded by the endpoints and the optimized candidate path may differ. Inspect the structures and the bond report first. If justified, compare `peak` and `minima` refinement or supply a verified intermediate; changing the refinement rule is not itself a fix. `all` has no `--refine-mode`; with `--refine-path`, set `search.refine_mode` in the `--config` YAML.
- **More reactive segments than input pairs** (`--refine-path` only). This is a candidate decomposition, not proof that the hidden intermediates are real; validate each IM and its TS/IRC.
- **n_imag ≥ 2.** A higher-order saddle or a numerical or seed problem, not a validated TS. Inspect the modes, improve the MEP seed, or try `--flatten`; use Dimer (`--opt-mode-post grad`) or a second backend as a cross-check rather than calling it a first-order saddle. See [Wrong n_imag after tsopt](../pdb2reaction-overview/ts-strategy.md#wrong-n_imag-after-tsopt).
- **Different atoms or order across inputs.** Re-extract with one controlled selection and compare the ordered (element, chain, residue, atom name) lists; equal PDB line counts are not enough.
- **GSM or DMF.** The better choice depends on the system and environment. GSM is the default; `--mep-mode dmf` needs cyipopt ([Core package](../pdb2reaction-install/backends.md#core-package)). Inspect and validate either MEP.
- **`-s` with several inputs.** It stops with an error; `-s` takes exactly one structure ([all-scan-list.md](all-scan-list.md)).

## Next step

- [all.md](all.md): mode choice, success criteria, resume, output tree.
- [path.md](path.md): what `path-opt` and `path-search` do.
- [bond-summary](utilities.md#bond-summary): what bond-change detection reports.
- [Reading outputs](../pdb2reaction-overview/outputs.md#per-segment-keys): multi-segment results.
