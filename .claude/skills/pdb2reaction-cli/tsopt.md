# `pdb2reaction tsopt`

## Purpose

Transition-state optimization. The default algorithm is RS-I-RFO
(full-Hessian; `--opt-mode hess`/`rsirfo`). Hessian-Guided Dimer
(`--opt-mode grad`/`dimer`) is the **alternative** TS optimizer when
RS-I-RFO struggles (e.g. unstable full-Hessian eigenstructure or very
large clusters where full-Hessian recomputation is prohibitive); the
dimer still uses an initial Hessian to set the search direction, then
updates the lowest mode via dimer rotation rather than recomputing the
full Hessian each cycle. Use after `path-search` or `scan` to refine a
HEI to a true first-order saddle, or as a standalone validator on an
externally-generated TS guess.

## Synopsis

```bash
pdb2reaction tsopt -i ts_guess.{pdb,xyz,gjf} \
    [-q 0 -m 1] [-l 'RES:Q,...'] \
    [--opt-mode grad|hess|dimer|rsirfo] \
    [--max-cycles 10000] \
    [-b uma|orb|mace|aimnet2] [-o ./result_tsopt/]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | TS candidate; `.pdb` / `.xyz` / `.gjf` |
| `-q` / `-l` / `-m` | — | — | Charge / spin (common conventions) |
| `--opt-mode` | str | `hess` | `grad`/`dimer` (Hessian-Guided Dimer) or `hess`/`rsirfo` (RS-I-RFO) |
| `--max-cycles` | int | 10000 | Optimization step cap |
| `--hessian-calc-mode` | str | (live default) | `Analytical` or `FiniteDifference` (default: `FiniteDifference`); selects how the initial Hessian is computed |
| `-b, --backend` | str | `uma` | MLIP backend |
| `--solvent` | str | none | xTB-ALPB solvent |
| `-o, --out-dir` | path | `./result_tsopt/` | Output directory |
| `--config` / `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

## Examples

### Default RS-I-RFO

```bash
pdb2reaction tsopt -i hei.xyz -q 0 -m 1 -b uma --out-json -o result_tsopt
```

`--out-json` is required for the `result.json` shown below; omit it if you only want the geometry.

### Dimer mode (alternative when RS-I-RFO struggles)

```bash
pdb2reaction tsopt -i hei.xyz -q 0 -m 1 \
    --opt-mode dimer -b uma -o result_tsopt_dimer
```

### Tighter convergence on an ill-conditioned saddle

```bash
pdb2reaction tsopt -i hei.xyz -l 'SAM:1,GPP:-3' \
    --opt-mode rsirfo --max-cycles 200 -b mace \
    -o result_tsopt_rsirfo
```

## Output

```
result_tsopt/
├── result.json                # only when --out-json is passed
├── final_geometry.xyz              # converged TS (always)
├── final_geometry.pdb              # PDB companion (with --convert-files)
├── final_geometry.gjf              # Gaussian companion (when input is .gjf)
├── optimization_trj.xyz            # full optimization trajectory (RS-I-RFO; with --dump)
├── optimization_all_trj.xyz        # full optimization trajectory (HessianDimer mode; with --dump)
└── vib/                            # imaginary-mode vibrations
    └── imag_*.{pdb,xyz}            # mode displacement visualization
```

`result.json` (only when `--out-json` is passed) keys:

```python
import json
d = json.load(open("result_tsopt/result.json"))
print(d["status"])                      # "converged" / "not_converged" (Dimer always "converged"; "error" on failure)
print(d["energy_hartree"])
print(d["n_imaginary_modes"])           # should be 1 for a real TS
print(d["imaginary_frequencies_cm"])    # list of cm⁻¹
print(d["files"]["final_geometry_xyz"]) # path under out_dir
```

## `--opt-mode` choice

| Mode | Algorithm | When |
|---|---|---|
| `hess` / `rsirfo` (default) | RS-I-RFO with full Hessian | Default. Robust for tricky / multi-imaginary-mode candidates; slower per cycle but converges in fewer cycles |
| `grad` / `dimer` | Hessian-Guided Dimer | Alternative when RS-I-RFO fails to converge or full-Hessian recomputation is prohibitive (e.g. large clusters > 600 atoms with UMA-m). Dimer uses the initial Hessian to set the search direction and then rotates the dimer pair rather than recomputing the full Hessian each cycle. |

Try `rsirfo` first; switch to `dimer` only if RS-I-RFO does not
converge or the full-Hessian cost becomes prohibitive.

## Validation: imaginary modes

A real TS has exactly one imaginary frequency that corresponds to the
reaction coordinate.

```python
import json
d = json.load(open("result_tsopt/result.json"))
if d["n_imaginary_modes"] == 1:
    print("OK: single imaginary mode at", d["imaginary_frequencies_cm"][0], "cm-1")
elif d["n_imaginary_modes"] == 0:
    print("BAD: collapsed to a minimum during refinement")
elif d["n_imaginary_modes"] > 1:
    print("AMBIGUOUS: multiple imaginary modes; inspect vib/imag_*.pdb")
```

For multi-imaginary cases, visualize the modes (`pymol vib/imag_*.pdb`)
to decide whether the extra modes are spurious (translation/rotation
of frozen residues) or real chemical second-order saddle points.

## Caveats

- A converged `tsopt` is **not** a complete validation; always follow
  with `irc.md` to confirm the TS connects the expected R and P.
- `--max-cycles` defaults to 10000 as a safety upper bound;
  well-conditioned cases converge in well under 200 cycles, so hitting
  >200 cycles usually means the TS guess is too far off — re-run
  `path-search` rather than raising the cap.
- Backend choice matters here more than for minima: UMA / MACE are
  usually safer than Orb for TS curvature.

## See also

- `path-search.md` — produces TS candidates for `tsopt`.
- `irc.md`, `freq.md` — downstream validation.
- [`pdb2reaction-install-backends/uma.md`](../pdb2reaction-install-backends/uma.md) / `mace.md` — TS-accurate
  backends.
- Defaults: `import pdb2reaction.defaults as d; print(d.RSIRFO_KW, d.DIMER_KW, d.HESSIAN_DIMER_KW)`
