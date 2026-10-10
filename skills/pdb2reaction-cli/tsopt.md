# `pdb2reaction tsopt`

## When to use

`tsopt` optimizes a TS candidate to a first-order saddle point, then computes
the Hessian at the final geometry and counts n_imag. A successful TS
optimization gives one imaginary mode along the reaction coordinate. Run it on
the HEI from `path-opt` / `path-search`, on the top of a `scan`, or on a
candidate you built. The default optimizer is RS-P-RFO; the Hessian-guided
Dimer is the alternative.

## Minimal run

```bash
pdb2reaction tsopt -i hei.xyz -q 0 -m 1 -b uma --out-json -o result_tsopt
```

Success: the console prints `[tsopt] Converged (n_imag=1).`.

Dimer, when RS-P-RFO does not converge:

```bash
pdb2reaction tsopt -i hei.xyz -q 0 -m 1 \
    --opt-mode dimer -b uma -o result_tsopt_dimer
```

Another RFO-family optimizer on a difficult candidate:

```bash
pdb2reaction tsopt -i hei.pdb -l 'SAM:1,GPP:-3' \
    --opt-mode rsirfo --max-cycles 200 -b mace \
    -o result_tsopt_rsirfo
```

## Judge success

The `[tsopt]` verdict line tells how the run ended:

| Verdict line | Meaning |
|---|---|
| `[tsopt] Converged (n_imag=1).` | First-order saddle. Check that the mode moves the reacting atoms, then run IRC |
| `[tsopt] WARNING: Higher-order stationary point (n_imag=N, …)` | Converged with extra imaginary modes; not yet a TS |
| `[tsopt] No imaginary mode detected. …` | Converged toward a minimum |
| `[tsopt] ERROR: Not converged (plateau stop, n_imag=N).` | Stopped on an energy plateau (`--stop-plateau`). Not converged, but the Hessian is always computed and n_imag reported |
| `[tsopt] ERROR: Not converged.` | Reached `--max-cycles`. No Hessian is computed, so n_imag is not reported |
| `[tsopt] Converged; terminal PHVA is unavailable.` | `--skip-final-freq` or a failed final Hessian; the saddle order is unchecked |

A mode counts as imaginary below −5 cm⁻¹, the tsopt cutoff (YAML
`freq.zero_cutoff_cm`). The warning `[tsopt] WARNING: the leading imaginary
mode is … cm^-1, below 50 cm^-1` changes neither the verdict nor n_imag; judge
such a TS by its mode and IRC ends like any other. With `--out-json`, a
first-order TS needs all of:
`optimization_status` is `converged`, `hessian_status` is `completed`,
`saddle_validation` is `first_order` (`n_imaginary_modes` is 1), and the
displacement in `vib/imag_*_trj.xyz` follows the reacting atoms.
`scientific_status` does not look at n_imag.

```python
import json
d = json.load(open("result_tsopt/result.json"))
print(d["optimization_status"])      # "converged" / "not_converged" / "stalled"
print(d["saddle_validation"])        # "first_order" / "higher_order" / "no_imaginary" / "unavailable"
print(d["n_imaginary_modes"], d["imaginary_frequencies_cm"])
print(d["energy_hartree"], d["files"]["final_geometry_xyz"])
print(d["reaction_mode_index"], d["reaction_mode_frequency_cm"],
      d["reaction_mode_source"], d["reaction_mode_overlap"])
```

`reaction_mode_index` and `reaction_mode_frequency_cm` name the imaginary mode
that `all` follows into IRC. With a reference direction (the MEP tangent in
`all` unless `--no-tsopt-from-mep-tan`, or `--ref-mode`), it is the imaginary
mode closest to that direction, `reaction_mode_source` is
`"mep-reference-overlap"`, and `reaction_mode_overlap` gives the overlap.
Otherwise, or when the final
Hessian was computed again, it is the lowest imaginary mode and the source is
`"lowest-imaginary"`, as always for Dimer and for `tsopt` without
`--ref-mode`. Neither value shows that the mode is the reaction: check the
`vib/imag_*_trj.xyz` of that frequency as above.

With `higher_order`, `optimization_status` still reports only the optimizer.
Watch every `vib/imag_*_trj.xyz` and do not accept the structure as a
first-order saddle. `all` may continue a warning-labelled diagnostic IRC from
it; that is not a first-order TS. Even a first-order TS candidate does not by
itself establish the intended elementary reaction: confirm it with `irc`.

## Choosing --opt-mode

| Mode | Algorithm | When |
|---|---|---|
| `hess` / `rsprfo` (default) | RS-P-RFO, full Hessian | Default; start here |
| `grad` / `dimer` | Hessian-guided Dimer | RS-P-RFO does not converge, or recomputing the full Hessian is too costly on a large cluster. It follows the lowest mode by dimer rotation and refreshes it from an exact Hessian only at intervals |
| `rsirfo` | RS-I-RFO | Another try on a difficult candidate |
| `trim` | TRIM | Another try on a difficult candidate |

`rsirfo` and `trim` are their own values, not aliases of `grad` or `hess`. In
`all`, `--opt-mode-post grad` selects the same Dimer.

## Pitfalls and recovery

- **Not converged at `--max-cycles`.** The default of 100000 is a safety
  bound; more cycles rarely help. Inspect the trajectory (`--dump`) and the
  candidate, then switch `--opt-mode` or start from a better candidate.
- **Stopped by the scheduler.** A run ended by walltime, a node failure, or
  cancellation is unfinished, not a convergence failure: it prints no
  `[tsopt]` verdict line and reports no n_imag. Do not count the candidate as
  one that does not converge or change the method for it; rerun with more
  walltime, or resume `all` ([all.md](all.md#resume-a-failed-segment)).
- **RS-P-RFO stops with a `ValueError`.** `RS-P-RFO exhausted its micro
  cycles outside the trust radius.`, `RS-P-RFO alpha update is not finite and
  positive.`, and `RS-P-RFO combined step exceeds the trust radius.` are
  numerical stops of the step solver: the run exits with 1 and a traceback,
  without a final Hessian or n_imag, and more cycles do not help. They do not
  show that the candidate is bad. Rerun once; if the stop repeats, switch
  `--opt-mode` to `rsirfo` or `dimer`, or start from another candidate.
- **n_imag ≥ 2.** Diagnose constraints, precision, and the character of each
  mode, then re-run with `--flatten` or get a better candidate. See
  [Wrong n_imag after tsopt](../pdb2reaction-overview/ts-strategy.md#wrong-n_imag-after-tsopt).
- **n_imag = 0.** The candidate is not near a saddle. From `all`,
  `--refine-path` gives a finer HEI, but recursive segmentation can multiply
  the MEP, TS, IRC, and freq cost, so it is off by default; inspect the coarse
  MEP before enabling it.
- **Precision.** Curvature is precision-sensitive on every backend. Keep the
  backend's precision (UMA fp32, ORB and MACE fp64 by default), inspect the
  actual mode, and check the final Hessian and IRC rather than ranking a saddle
  by optimizer convergence alone. Use standalone `freq` for more analysis.
- **`--ref-mode`.** This advanced option gives a Cartesian 3N reaction direction
  for choosing the initial root and tracking overlap. `all` supplies it from
  the MEP; with `all --no-tsopt-from-mep-tan`, its TS optimization chooses
  from the Hessian modes of the starting structure. Ordinary standalone `tsopt` runs should omit
  it. Do not add it because a run is difficult; it is not a convergence switch.
  Supply it only when an external path gives a deliberate non-zero direction.
- **Frozen atoms.** With frozen atoms, the Hessian covers the movable atoms
  and only rigid motions that keep every frozen atom in place are removed. This
  is separate from `--ref-mode`. See
  [PHVA treatment](extract.md#phva-treatment)
  and the [JSON record](../../docs/json-output.md#rigid-projection-provenance).
- **Uphill steps.** `tsopt` always forces `reject_uphill=False`, whatever the
  optimizer mode or YAML, because uphill steps can be part of following a
  saddle mode. `--reject-uphill` belongs to `opt` and to endpoint optimization
  in `all`.
- **Analytical Hessian with UMA workers.** See [Cross-cutting pitfalls](SKILL.md#cross-cutting-pitfalls).

## Outputs

- `final_geometry.xyz`: always written. `--convert-files` (on by default) adds
  `.pdb` for PDB/mmCIF input, `.cif` with the original IDs for mmCIF or very
  large PDB input, and `.gjf` for Gaussian input.
- `vib/imag_<freq>cm-1_trj.xyz`: animation of each imaginary mode, with
  `.pdb` / `.cif` when the input carries topology.
- `--dump`: `optimization_trj.xyz` (RS-P-RFO, RS-I-RFO, TRIM) or
  `optimization_all_trj.xyz` (Dimer).
- `--dump-hess <file>.npy`: the final Hessian, written only when it was
  computed, for `--read-hess` in `freq`, `tsopt`, or `irc`.
- `result.json`: with `--out-json`.

## Next step

- [irc.md](irc.md): confirm that the TS connects the intended R and P.
- [freq.md](freq.md): thermochemistry or more mode analysis.
- [path.md](path.md): make TS candidates.
- Backend notes: [UMA](../pdb2reaction-install/backends.md#uma),
  [MACE](../pdb2reaction-install/backends.md#mace-separate-environment).
- Defaults: `import pdb2reaction.core.defaults as d; print(d.RSIRFO_KW, d.DIMER_KW, d.HESSIAN_DIMER_KW)`
