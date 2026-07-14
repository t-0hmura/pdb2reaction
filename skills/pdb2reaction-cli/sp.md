# `pdb2reaction sp`

## Purpose

Evaluate one MLIP energy and force array, with an optional Hessian. Use
this for backend/precision checks and stationary-point diagnostics without
starting an optimizer. Use `dft` instead for a PySCF/GPU4PySCF single point.

## Synopsis

```bash
pdb2reaction sp -i geom.{pdb,cif,mmcif,xyz,gjf} \
    [-q 0 -m 1] [-l 'RES:Q,...'] \
    [--hess --hessian-calc-mode Analytical] \
    [-b uma|orb|mace|aimnet2] [-o ./result_sp/]
```

## Key flags

| flag | default | description |
|---|---|---|
| `-i, --input` | required | Input PDB, mmCIF, XYZ, or GJF |
| `-q` / `-l` / `-m` | resolved by common rules | Total charge / per-residue charge mapping / multiplicity |
| `--hess / --no-hess` | off | Also write a Hessian; without frozen atoms it is `(3N, 3N)`, while YAML `geom.freeze_atoms` selects the active partial block |
| `--hessian-calc-mode` | automatic for `sp --hess` | Automatic mode uses Analytical for UMA and FiniteDifference for ORB/MACE/AIMNet2; an explicit `Analytical` request is supported by all four backends |
| `--precision` | backend-dependent | Unset means UMA/AIMNet2 fp32 and ORB/MACE fp64; AIMNet2 rejects fp64 |
| `--workers` / `--workers-per-node` | `1` / `1` | UMA predictor parallelism; other built-in backends filter these keys |
| `-b, --backend` | `uma` | MLIP backend |
| `-o, --out-dir` | `./result_sp/` | Output directory |
| `--out-json` | off | Write identical `result.json` and `summary.json` payloads |
| `--show-config` / `--dry-run` | off | Print-and-continue / validate-and-exit |

With UMA, `--workers` greater than 1 cannot honor an explicit Analytical
Hessian request and raises `BackendError`; select FiniteDifference or use one
worker. No Hessian restriction applies when `--hess` is absent.

## Examples

```bash
# Energy and forces only
pdb2reaction sp -i geom.pdb -l 'SAM:1,GPP:-3' -b orb -o result_sp

# Trusted ORB Hessian: the unset precision resolves to fp64
pdb2reaction sp -i geom.xyz -q -2 -m 1 -b orb \
    --hess --hessian-calc-mode Analytical --out-json -o result_sp_hess

# Numerical cross-check
pdb2reaction sp -i geom.xyz -q -2 -m 1 -b mace \
    --hess --hessian-calc-mode FiniteDifference -o result_sp_fd
```

## Output

| Path | When | Content |
|---|---|---|
| `<out_dir>/forces.npy` | successful calculation | `(N, 3)` forces in Hartree/Bohr |
| `<out_dir>/hessian.npy` | `--hess` | mass-unweighted Hessian in Hartree/Bohr² (full or YAML-freeze active block) |
| `<out_dir>/result.json`, `summary.json` | `--out-json` | status, energy, provenance, charge/spin, and file paths |
| stdout | successful single-point calculation | energy in Hartree, maximum force, and elapsed time |

The JSON payload uses `status: "ok"`, `energy_au`, `forces_path`, and
`hessian_path`. `hessian_path` is null when `--hess` was not requested.
`sp` does not write `summary.log`.

## Validation

- Compare Hessians only at the same geometry, atom ordering, freeze set,
  backend model, precision, charge, and multiplicity.
- An fp64 request is not supported by AIMNet2; it fails instead of silently
  changing precision.
- `--show-config` does not stop execution. Use `--dry-run` for a preview.

## See also

- `freq.md` — frequencies, PHVA, and thermochemistry.
- `tsopt.md` — transition-state optimization and saddle validation.
- [`pdb2reaction-install-backends/SKILL.md`](../pdb2reaction-install-backends/SKILL.md) — backend installation and precision behavior.
