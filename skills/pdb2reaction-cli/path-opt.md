# `pdb2reaction path-opt`

## Purpose

MEP optimization for **one** segment between two endpoints. It exposes the
same GSM / optional-DMF path engines that the recursive `path-search` workflow
uses, but is a separate command: it does not recursively split the path or
classify bond-change segments. Use it standalone when one endpoint pair is the
intended candidate unit; the result still needs TS/frequency/IRC validation
before calling it elementary.

## Synopsis

```bash
pdb2reaction path-opt -i reactant.pdb product.pdb \
    [--mep-mode gsm|dmf] [--max-nodes 20] \
    [-l 'RES:Q,...'] [-b uma|orb|mace|aimnet2] \
    [-o ./result_path_opt/]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path(s) | required (= 2) | Reactant and product with identical atom identities and ordering |
| `--mep-mode` | str | `gsm` | `gsm` (Growing String) or `dmf` (Direct Max Flux) |
| `--max-nodes` | int | 20 | Max internal nodes (final string ≤ `max-nodes + 2`) |
| `--thresh` | str | `gau` | Endpoint preoptimization convergence preset |
| `--thresh-gsm` | str | `gau_loose` | GSM string-optimizer convergence preset |
| `--thresh-dmf` | str/float | `tight` | DMF IPOPT dual-infeasibility tolerance: `tight`, `middle`, `loose`, or a positive float |
| `--preopt / --no-preopt` | flag | `--preopt` | Optimize each endpoint before constructing the string |
| `--fix-ends / --no-fix-ends` | flag | `--fix-ends` | Keep GSM endpoint images fixed; accepted but unused with DMF |
| `-q, --charge` / `-l` / `-m` | — | — | Charge / spin (common conventions) |
| `-b, --backend` | str | `uma` | MLIP backend |
| `--solvent` | str | none | xTB-ALPB solvent name |
| `-o, --out-dir` | path | `./result_path_opt/` | Output directory |

## Examples

### Default GSM single segment

```bash
pdb2reaction path-opt -i R.xyz P.xyz -q 0 -m 1 -b uma -o result_path_opt
```

### Optional DMF engine

```bash
pdb2reaction path-opt -i R.pdb P.pdb -l 'GPP:-3' --mep-mode dmf -b mace \
    -o result_path_opt_dmf
```

DMF is optional and cannot be combined with `--solvent`: its ASE path has no
xTB solvent wrapper, so pdb2reaction rejects that PES mismatch. Use GSM for a
solvent-corrected path.

## Output

| Path | When | Content |
|---|---|---|
| `<out_dir>/result.json` | `--out-json` | machine-readable result |
| `<out_dir>/final_geometries_trj.xyz` | completed path optimizer run | final string trajectory; inspect JSON `status` before calling it converged |
| `<out_dir>/final_geometries.{pdb,cif,gjf}` | `--convert-files` (default on) AND topology/template present; CIF needs a bridged input/reference | PDB / CIF / GJF companions |
| `<out_dir>/hei.{xyz,pdb,cif,gjf}` | path result reaches HEI writing; non-XYZ companions need a topology/template plus `--convert-files` | highest-energy image (TS candidate) |
| `<out_dir>/dmf_initial_trj.xyz`, `dmf_ipopt.out` | `--mep-mode dmf` | DMF-mode artifacts |

`result.json` reports the final string energies and status. The status is
`converged` / `not_converged` when the selected engine exposes a convergence
flag, or `completed` when that engine completed but did not expose one.
`completed` must not be relabeled as convergence.

## When to use vs path-search

- **`path-search`** if you want recursive bond-change segmentation
  for a possibly multi-step mechanism.
- **`path-opt`** if one endpoint pair is the intended candidate unit and you
  want its MEP without recursive segmentation. TS/frequency/IRC must still
  establish whether it is one elementary step.

## Caveats

- Convergence is sensitive to endpoint geometries. Endpoint preoptimization is
  already enabled by default; if `not_converged`, inspect its result and the
  string before deciding whether separate `opt`, a different endpoint, or a
  path-engine setting is appropriate.
- If GSM stalls after endpoint/preoptimization checks, DMF is an optional
  alternative where its separate dependencies are installed. Changing the
  engine or raising `--max-nodes` increases cost and does not repair an
  inconsistent endpoint pair; inspect the resulting path rather than treating
  engine completion as validation.

## See also

- `path-search.md` — the recursive driver around this command.
- `opt.md` — pre-relax endpoints before path-opt.
- Defaults: `import pdb2reaction.core.defaults as d; print(d.GS_KW, d.DMF_KW, d.STOPT_KW)`
