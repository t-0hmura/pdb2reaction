# `pdb2reaction all` — TS-only mode

## When to use

You already have a **TS candidate** (typically from another QM code, an
older `pdb2reaction` run, or a manual guess) and want to run only the
validation stages without the upstream extract/path search: `tsopt → irc`,
plus R/TS/P `freq` when `--thermo` and DFT when `--dft`.

## Synopsis

```bash
pdb2reaction all -i ts_candidate.xyz \
    -q -1 -m 1 -b uma \
    --tsopt --thermo \
    [--dft --dft-func-basis 'wb97m-v/def2-tzvpd'] \
    -o result_ts_only
```

Or with a residue-labeled PDB/mmCIF from which `-l` can derive the charge:

```bash
pdb2reaction all -i ts_candidate.pdb \
    -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    -o result_ts_only
```

## How it differs from the other two modes

`pdb2reaction all` falls into TS-only mode when **all three** hold:

- exactly **one** `-i` input is given,
- **no** `--scan-lists` is provided,
- `--tsopt` is passed.

Without `--scan-lists` or `--tsopt` the CLI raises `BadParameter`
(`all.md` covers the orchestrator's input gate).

For finer control, run the underlying subcommands directly:

```bash
TOTAL_CHARGE=-1  # replace with the verified cluster charge
pdb2reaction tsopt -i ts.xyz -q "$TOTAL_CHARGE" -m 1 -o result_tsopt -b uma
pdb2reaction irc   -i result_tsopt/final_geometry.xyz -q "$TOTAL_CHARGE" -m 1 -o result_irc -b uma
pdb2reaction freq  -i result_tsopt/final_geometry.xyz -q "$TOTAL_CHARGE" -m 1 -o result_freq -b uma
```

`extract` and MEP search are skipped entirely. The output tree for a
completed run is:

| Path | When | Content |
|---|---|---|
| `<out_dir>/summary.json` | pipeline reaches its summary writer | machine-readable result; early CLI/input validation can fail before this file exists |
| `<out_dir>/summary.log` | pipeline reaches its summary writer | human-readable text + dir tree; early CLI/input validation can fail before this file exists |
| `<out_dir>/segments/seg_01/{reactant,ts,product}.{pdb,xyz}` | successful TSOPT + IRC/endpoint processing | canonical R/TS/P (extension follows the input/topology available) |
| `<out_dir>/segments/seg_01/structures/` | successful TSOPT + IRC/endpoint processing | `.xyz` working copies of R/TS/P, `.pdb` when topology is available, and `.cif` when the bridge retained public IDs |
| `<out_dir>/segments/seg_01/ts/final_geometry.xyz`; `.pdb`/`.cif` companions; `optimization_trj.xyz` with `--dump` | TS optimizer reaches output; companions require `--convert-files` and topology/bridge metadata | final TS attempt and optional trajectory |
| `<out_dir>/segments/seg_01/irc/{forward,backward,finished}_irc_trj.xyz` | IRC completes | IRC trajectories |
| `<out_dir>/segments/seg_01/freq/{R,TS,P}/{frequencies_cm-1.txt, thermoanalysis.yaml}` | `--thermo` and each freq stage succeeds | per-state frequencies + thermochemistry |
| `<out_dir>/segments/seg_01/dft/{R,TS,P}/result.yaml` | `--dft` and the corresponding state calculation succeeds | per-state DFT (`result.json` only when `dft` is run standalone with `--out-json`) |

## Output keys

```python
import json
d = json.load(open("result_ts_only/summary.json"))
seg  = d["segments"][0]            # MEP-style block (barrier/delta/bond_changes)
post = d["post_segments"][0]       # post-processing block (ts_imag / mlip / dft …)
print(seg["barrier_kcal"], seg["delta_kcal"])
print(seg["bond_changes"])         # what bonds broke / formed along the IRC
print(post["ts_imag"]["n_imag"])   # should be 1 for a true TS
print(post["mlip"]["energies_au"])
```

If `post["ts_imag"]["n_imag"] != 1`, the geometry is **not a true first-order
saddle**; see "Distinctive failure modes" below.

`all` requires the TS stage `result.json` to report `status=converged` and
`n_imaginary_modes=1` before starting IRC. A rejected result leaves the TS
attempt for diagnosis and stops before endpoint post-processing.

## Distinctive failure modes

| Symptom | Likely cause | Fix |
|---|---|---|
| top-level `status == "partial"` and the segment lacks `mlip`/`ts_imag` post data | TS optimization or subsequent validation did not complete | Inspect `summary.log` and the `segments/seg_01/ts/` artifacts before choosing a targeted retry such as a better MEP seed, another optimizer family, coordinates, or `--flatten` for surplus modes. |
| `post["ts_imag"]["n_imag"] == 0` | Geometry reached a local minimum or a near-zero mode was classified differently | Inspect frequencies/displacements and obtain a better saddle seed, commonly from a validated/refined endpoint path. |
| `post["ts_imag"]["n_imag"] >= 2` | Higher-order saddle or numerical/constraint artifact | Not a valid first-order TS. Inspect every displacement, verify freeze/PHVA and precision, then retry from a better seed and/or with `--flatten`. |
| `segments[0]["bond_changes"]` is empty (no cutoff-defined covalent bonds change) | Endpoints may differ only conformationally/proton-position-wise, may be the same basin, or the geometric cutoff may miss the event | Inspect both endpoints and the imaginary displacement; an empty covalent bond-change report alone does not classify the TS as physical or non-physical. |

## When *not* to use TS-only mode

- You do not yet have a TS candidate. Run `path-search` (or the
  full `all` in endpoint-MEP / scan-list mode) instead.
- A TS candidate with unknown connectivity can still be tested by standalone
  `tsopt`/`freq`/`irc`; inspect both IRC endpoints. If verified R/P endpoints
  are available and the seed proves wrong, build/refine an endpoint MEP with
  `path-opt` or `path-search`.

## Caveats

- For an XYZ TS candidate, you must supply `-q` and `-m` explicitly
  (XYZ has no header). Use `--ref-pdb cluster.pdb` if you want
  `-l 'RES:Q'` to work.
- The IRC step here is the **canonical validation** that the TS connects the
  expected R and P. Always inspect `seg_01/{reactant,product}.xyz`; `.pdb` or
  `.cif` companions are available only when conversion topology exists.
- TS-only mode has no supplied R/P references with which to orient the two
  IRC ends. The current pipeline names the higher-energy optimized endpoint
  `reactant` and the lower-energy one `product`; those names are an energy
  convention, not chemical identity. Compare both structures with the
  intended states before reporting a forward barrier.

## See also

- `all.md` — base orientation.
- `tsopt.md`, `irc.md`, `freq.md`, `dft.md` — the underlying
  subcommands (which you can also run standalone if you want
  fine-grained control).
- `pdb2reaction-workflows-output/SKILL.md` — IRC interpretation
  and bond-change conventions.
