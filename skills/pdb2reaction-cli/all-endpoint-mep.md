# `pdb2reaction all` — endpoint-MEP mode

## When to use

You have **two or more reaction-ordered structures** (reactant, optional
intermediate(s), product), all with the **same atom count and atom
ordering**. The pipeline interpolates an MEP between adjacent
structures and segments multi-step paths automatically.

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

The recursive bond-change segmentation in `path-search` will still run
inside each adjacent pair, so you don't have to provide every
elementary step — just the "obvious" ones from the literature.

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

If the inputs come from different programs or were re-numbered, run
them through `extract` once to canonicalize ordering:

```bash
pdb2reaction extract -i 1.R_raw.pdb -i 3.P_raw.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -o 1.R.pdb -o 3.P.pdb
```

## Output

Same as `all.md`. Specifically for endpoint-MEP mode:

| Path | When | Content |
|---|---|---|
| `<out_dir>/mep.pdb`, `mep_trj.xyz`, `mep_w_ref.pdb`, `energy_diagram_MEP.png` | always | stitched MEP (promoted from `_work/path_search/`) |
| `<out_dir>/_work/path_search/mep_seg_NN_trj.xyz`, `mep_seg_NN.{pdb,gjf}` | always | per-segment MEP frames (no bare `.xyz`); scratch |
| `<out_dir>/_work/path_search/seg_NNN_<tag>/` | always | recursive-splitter scratch (3-digit, `_mep` / `_maxdepth` / `_bridge`); internal |
| `<out_dir>/segments/seg_NN/` | always | per-segment deliverables + post-processing (ts/irc/freq/dft) + energy diagrams |
| `<out_dir>/segments/seg_NN/{reactant,ts,product}.pdb` | always | canonical R/TS/P (after IRC + LBFGS endpoint opt) |
| `<out_dir>/summary.json["segments"]` | always | list of `{index, barrier_kcal, delta_kcal, bond_changes, ...}` |

## Distinctive failure modes

| Symptom in `summary.json` | Likely cause | Fix |
|---|---|---|
| `path_search.status == "partial"` with `bond-summary` reporting extra changes between inputs vs the optimized MEP | Bond-change detector found extra changes; the reaction in the inputs and the reaction the optimizer found don't match. | Check which bonds changed via `bond-summary -i 1.R.pdb 3.P.pdb`; rerun the standalone `pdb2reaction path-search` with `--refine-mode minima`, or supply the missing intermediate(s) explicitly via additional `-i` inputs. |
| `tsopt.n_imaginary_modes > 1` for a segment | Multi-imaginary-mode TS — common with Orb on tricky systems | Re-run that segment with `pdb2reaction tsopt --opt-mode dimer -b uma` (default is already RS-I-RFO; the dimer is the lighter alternative when full-Hessian RS-I-RFO is unstable). |
| Different atom counts across `-i` inputs | Inconsistent extractions | Re-extract per the snippet above, verify with `wc -l 1.R.pdb 3.P.pdb`. |

## Caveats

- The "right" choice between `--mep-mode gsm` and `--mep-mode dmf`
  depends on system size and how much you trust the initial MEP guess;
  GSM is the safer default.
- Path search may discover **more** segments than you have inputs:
  if `summary.json["n_segments"] > len(inputs) - 1`, that's the
  recursive bond-change segmentation finding intermediates the inputs
  didn't contain. Often this is the *correct* answer.

## See also

- `all.md` — base orientation (output tree, summary.json schema).
- `path-search.md` — recursive MEP search internals.
- `bond-summary.md` — what bond-change detection looks like.
- `pdb2reaction-workflows-output/SKILL.md` — interpreting multi-segment
  results.