# `pdb2reaction all` — endpoint-MEP mode

## When to use

You have **two or more reaction-ordered structures** (reactant, optional
intermediate(s), product), all with the **same atom count and atom
ordering**. The default pipeline optimizes one MEP per adjacent input pair.
Add `--refine-path` when recursive bond-change segmentation and hidden
intermediate discovery are desired.

This is the most common mode for a published-mechanism reproduction
where you have R and P (and sometimes IM) coordinates from a prior QM
or QM/MM study.

## Synopsis

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    [--dft --dft-func-basis 'wb97m-v/def2-tzvpd'] \
    -o result_mep
```

For a known multistep mechanism, supply each intermediate explicitly:

```bash
pdb2reaction all -i 1.R.pdb 2.IM.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    -o result_mep_3pt
```

With `--refine-path`, the recursive bond-change segmentation in
`path-search` runs **once over the full ordered series of inputs**, so you
don't have to provide every elementary step — just the "obvious" ones from
the literature. By default `all` runs single-pass `path-opt` between adjacent
pairs only and does not auto-discover hidden intermediates.

## Mode-specific flags

| Flag | Default | Meaning |
|---|---|---|
| `--mep-mode` | `gsm` | `gsm` (Growing String) or `dmf` (Direct Max Flux) |
| `--max-nodes` | 20 | Maximum string nodes per segment (final string ≤ `max-nodes + 2`) |

> `--refine-mode {peak,minima}` is exposed only on the standalone
> `pdb2reaction path-search` subcommand, not on `all`. To override it
> from `all`, add a small YAML config (`--config foo.yaml`) with a
> `search.refine_mode` key.

`--scan-lists` is **not** allowed in this mode — it triggers
`all-scan-list.md` instead.

## Atom-count consistency requirement

All `-i` inputs must have:

- the same number of atoms,
- the same element sequence (atom ordering),
- the same residue assignments.

`extract` validates multi-input identity/order but does not atom-map or repair
an already mismatched series. If inputs came from different programs or were
renumbered, first map/reorder them to a common full-structure signature; then
use one multi-input extraction to apply and verify the same cluster boundary:

```bash
pdb2reaction extract -i 1.R_raw.pdb -i 3.P_raw.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -o 1.R.pdb -o 3.P.pdb
```

## Output

Same as `all.md`. Specifically for endpoint-MEP mode:

| Path | When | Content |
|---|---|---|
| `<out_dir>/mep_trj.xyz`, `energy_diagram_MEP.png`; PDB/CIF/GJF companions when topology/template exists | successful MEP | stitched MEP (promoted from `<work_path>` = `_work/path_opt/` by default, `_work/path_search/` with `--refine-path`) |
| `<out_dir>/_work/path_opt/mep_seg_NN_trj.xyz`, `mep_seg_NN.{pdb,gjf}` | successful MEP (`_work/path_search/` with `--refine-path`) | per-segment MEP frames (no bare `.xyz`); scratch |
| `<out_dir>/_work/path_search/seg_NNN_<tag>/` | `--refine-path` | recursive-splitter scratch (3-digit, `_mep` / `_maxdepth` / `_bridge`); internal |
| `<out_dir>/segments/seg_NN/` | a reactive segment enters requested post-processing | per-segment outputs; may be partial after a failed stage, so inspect JSON status |
| `<out_dir>/segments/seg_NN/{reactant,ts,product}.{pdb,xyz}` | successful `--tsopt` + IRC/endpoint processing | canonical R/TS/P (after IRC + RFO endpoint opt; LBFGS via `--opt-mode-post grad`) |
| `<out_dir>/summary.json["segments"]` | pipeline reaches its summary writer | list of `{index, barrier_kcal, delta_kcal, bond_changes, ...}`; early validation can produce no summary |

## Distinctive failure modes

| Symptom in `summary.json` | Likely cause | Fix |
|---|---|---|
| `status == "partial"` with `bond-summary` reporting extra changes between inputs vs the optimized MEP | Bond-change detector found extra changes; the reaction encoded by the endpoints and the optimized candidate path may differ. | Inspect the structures and bond report first. If justified, compare `peak`/`minima` refinement or supply a verified intermediate; changing the refinement rule is not itself a correctness fix. |
| `post_segments[i]["ts_imag"]["n_imag"] > 1` | Higher-order saddle or numerical/seed problem; not a validated TS | Inspect the modes, improve the MEP seed and/or try `--flatten`; use Dimer or a second backend as a cross-check rather than relabeling it a first-order saddle. |
| Different atoms/order across `-i` inputs | Inconsistent extractions | Re-extract with one controlled residue/atom selection and compare the ordered `(element, chain, residue, atom-name)` signatures; equal PDB line counts alone are insufficient. |

## Caveats

- The appropriate choice between `--mep-mode gsm` and `--mep-mode dmf`
  is system- and environment-dependent. GSM is the software default; DMF is
  optional and needs its separate environment. Inspect and validate either MEP.
- With `--refine-path`, path search may propose **more reactive segments**
  than adjacent input pairs. Compare `summary.json["n_segments_reactive"]`
  with `len(inputs) - 1`; top-level `n_segments` also counts bridge segments.
  Any extra reactive segmentation is a candidate decomposition, not proof
  that hidden intermediates are chemically correct. The default
  single-pass `path-opt` yields one segment per adjacent input pair.

## See also

- `all.md` — base orientation (output tree, summary.json schema).
- `path-search.md` — recursive MEP search internals.
- `bond-summary.md` — what bond-change detection looks like.
- `pdb2reaction-workflows-output/SKILL.md` — interpreting multi-segment
  results.
