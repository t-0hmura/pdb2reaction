# `pdb2reaction all`: TS-only mode

Give one TS candidate with `--tsopt` and no `-s`; `all` optimizes the TS, runs
IRC, and optimizes both IRC ends (`--thermo` and `--dft` add R/TS/P frequencies
and DFT). It succeeded when the console prints
`[tsopt] Converged (n_imag=1).` and `Scientific status: success` under the last
`====== Pipeline summary ======`; then check that the IRC ends are the intended
R and P.

## When to use, and when not

Use it when you already have a TS candidate (from another QM code, an earlier
run such as `result_all/_work/path_opt/hei_seg_01.pdb`, or a manual guess) and
want only the validation stages, without an MEP search.

Without a TS candidate, use the [Multi-structure MEP search](all-endpoint-mep.md)
or [Single structure + scan](all-scan-list.md) mode, or `path-search`
([path.md](path.md)). A candidate of unknown connectivity can also be tested
with standalone `tsopt`, `freq`, and `irc`; inspect both IRC ends. If verified
R and P exist and the seed proves wrong, build an MEP between them.

## Minimal run

```bash
pdb2reaction all -i ts_candidate.xyz \
    -q -1 -m 1 -b uma \
    --tsopt --thermo \
    -o result_ts_only
```

With a residue-labelled PDB/mmCIF, `-l` derives the charge:

```bash
pdb2reaction all -i ts_candidate.pdb \
    -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    -o result_ts_only
```

Add `--dft` (and `--func-basis 'wb97m-v/def2-tzvpd'`) for DFT single points on
R, TS, and P. A PDB/mmCIF candidate is cut into a cluster only when `-c` is
given; otherwise it is used as is.

## How the mode is chosen

`all` runs TS-only mode for exactly one `-i` input with `--tsopt` and no
`-s`; `summary.log` shows `Pipeline mode` as `TS-only`. One input without
`-s` or `--tsopt` stops with `BadParameter`. One input with both `-s` and
`--tsopt` runs the scan mode instead.

## Judge success

- **TS**: a successful TS optimization gives one imaginary mode along the reaction coordinate. `post_segments[0].tsopt.n_imaginary_modes` should be 1 and `.imaginary_frequencies_cm` gives its wavenumber; play `segments/seg_01/ts/vib/imag_*_trj.xyz` to see that the mode moves the bonds that form or break. If the TS optimization stops unconverged, the run stops before IRC and keeps the TS files in `segments/seg_01/ts/`; n_imag is computed after a `--stop-plateau` stop but not at the cycle limit.
- **Status**: `scientific_status` is `success` only when every requested stage converged and n_imag = 1; otherwise read `scientific_status_reasons`.
- **Endpoints**: open `segments/seg_01/irc/finished_irc_trj.xyz` and `segments/seg_01/reactant.*` and `product.*`, and read `segments[0].bond_changes`. Even if the IRC does not converge, the result is usable when the endpoint optimizations reach the intended R and P.
- **R and P names**: with no MEP, the higher-energy IRC end is named the reactant (on an exact tie, the left end). The names and the barrier follow this energy order, not a known chemical direction; `post_segments[0].endpoint_assignment` records the rule with `chemical_direction_known: false`. The barrier from P is `barrier_kcal − delta_kcal`. Compare both ends with the intended states before reporting a forward barrier.
- **Energies**: `post_segments[0].mlip.barrier_kcal` and `.delta_kcal` (same values in `segments[0]`); `gibbs_mlip` (`--thermo`) and `dft` (`--dft`) carry the same keys.

```python
import json
d = json.load(open("result_ts_only/summary.json"))
seg, post = d["segments"][0], d["post_segments"][0]
print(d["scientific_status"], d.get("scientific_status_reasons"))
print(post["tsopt"]["n_imaginary_modes"], post["tsopt"]["imaginary_frequencies_cm"])
print(seg["barrier_kcal"], seg["delta_kcal"], seg["bond_changes"])
print(post["endpoint_assignment"], post["mlip"]["energies_au"])
```

## Pitfalls and recovery

- **`scientific_status` is `failed` and `post_segments[0]` lacks `tsopt` or `mlip`.** TS optimization or a later check did not finish. Read `summary.log` and `segments/seg_01/ts/` before a targeted retry: a better seed, Dimer (`--opt-mode-post grad`), other coordinates (`--coord-type`), or `--flatten` for extra modes.
- **n_imag = 0.** The run stops before IRC and is not `success`. The geometry reached a minimum, or a near-zero mode was classified differently; inspect the frequencies and displacements and start from a better seed, such as the HEI of a validated MEP.
- **n_imag ≥ 2.** The result is `partial`; IRC follows one mode only as a diagnostic, and the structure is not a first-order saddle. Inspect every displacement, check the frozen atoms, PHVA, and precision, then retry from a better seed or with `--flatten`. See [Wrong n_imag after tsopt](../pdb2reaction-overview/ts-strategy.md#wrong-n_imag-after-tsopt).
- **`bond_changes` is empty, or an end is not the intended state.** The ends may differ only in conformation or proton position, be the same basin, or the bond cutoff may miss the event; inspect both ends and the mode. An empty report alone does not judge the TS.
- **XYZ candidate.** The charge must come from `-q`, or from `--ref-pdb cluster.pdb` with `-l 'RES:Q'`. `-m` defaults to 1; set it for open-shell systems.
- **TS still not found.** See [When the TS does not come out](../pdb2reaction-overview/ts-strategy.md#when-the-ts-does-not-come-out).

## Run the stages yourself

For finer control:

```bash
TOTAL_CHARGE=-1  # replace with the verified cluster charge
pdb2reaction tsopt -i ts.xyz -q "$TOTAL_CHARGE" -m 1 -o result_tsopt -b uma
pdb2reaction irc   -i result_tsopt/final_geometry.xyz -q "$TOTAL_CHARGE" -m 1 -o result_irc -b uma
pdb2reaction freq  -i result_tsopt/final_geometry.xyz -q "$TOTAL_CHARGE" -m 1 -o result_freq -b uma
```

## Outputs

`summary.json` and `summary.log` sit at the top of `--out-dir`. Cite
`segments/seg_01/reactant.*`, `ts.*`, and `product.*` (in the input format);
`seg_01/structures/` holds working copies.
`seg_01/` also has `ts/` (`final_geometry.*`, `vib/imag_*_trj.xyz`,
`optimization_trj.xyz` with `--dump`), `irc/{forward,backward,finished}_irc_trj.xyz`,
`freq/{R,TS,P}/` with `frequencies_cm-1.txt` and `thermoanalysis.yaml`
(`--thermo`), `dft/{R,TS,P}/result.yaml` (`--dft`), and the energy diagrams.
There are no MEP files and no `_work/path_opt/`.

## Next step

- [all.md](all.md): mode choice, success criteria, resume.
- [tsopt.md](tsopt.md), [irc.md](irc.md), [freq.md](freq.md), [dft.md](dft.md): each stage on its own.
- [Reading outputs](../pdb2reaction-overview/outputs.md#bond-changes): IRC ends and bond changes.
