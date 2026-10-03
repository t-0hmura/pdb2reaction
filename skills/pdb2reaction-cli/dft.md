# `pdb2reaction dft`

## When to use

Run one DFT single point (energy, plus Mulliken, meta-Löwdin, and IAO charges
and spin densities) with GPU4PySCF (`--dft-engine gpu`, the default) or PySCF
(`--dft-engine cpu`), typically on the R, TS, and P from `tsopt` or `irc` to
get DFT//MLIP energies. `-b dft` is different: it makes DFT the calculator of
another command (`sp`, `opt`, `tsopt`, `irc`, `freq`, scan and path commands,
`all`), and those iterative runs reuse the last converged density.

## Minimal run

```bash
pdb2reaction dft -i seg_01/ts.pdb -l 'SAM:1,GPP:-3' \
    --func-basis 'wb97m-v/def2-tzvpd' --dft-engine gpu --out-json
```

The default method is `wb97m-v/def2-svp`. XYZ input needs `-q` and `-m`, or
`--ref-pdb` so that `-l` works. `--solvent` adds PySCF implicit solvent
(`--solvent-model smd` by default, or `pcm`); `--dft-nprocs` and
`--dft-memory` override the detected CPU thread count and host RAM limit.
On a machine without a usable GPU:

```bash
pdb2reaction dft -i ts.xyz -q 0 -m 1 --func-basis 'wb97m-v/def2-svp' \
    --dft-engine cpu -o result_dft_cpu
```

## Judge success

The console prints `E_total (Hartree): …`, the run exits with 0, and
`result_dft/result.yaml` has `energy.converged: true` together with the grid
level and the per-atom charges and spin densities. `input_geometry.xyz` is the
geometry sent to PySCF. With `--out-json`, `result.json` (and its copy
`summary.json`) is also written:

```python
import json
d = json.load(open("result_dft/result.json"))
print(d["energy_hartree"], d["converged"])  # exit code 1 if not converged
print(d["xc_functional"], d["basis_set"])   # e.g. "wb97m-v", "def2-tzvpd"
print(d["engine"])  # "gpu4pyscf(rks_lowmem)", "gpu4pyscf", or "pyscf(cpu)"
print(d["used_gpu"], d["used_lowmem"])  # lowmem is False for open shell, CPU, or --no-dft-low-memory
```

## Pitfalls and recovery

- An unconverged SCF prints `WARNING: SCF did not converge to the requested
  tolerance.`, writes `converged: false`, and exits with 1.
- `OSError: libcusolver.so.11 not found`: capture `pip check` and compare with
  the clean-environment library-loading test in
  [backends.md](../pdb2reaction-install-backends/backends.md#cuda-and-pytorch);
  do not guess a library path.
- `cupy ... invalid device ordinal`: keep the scheduler's
  `CUDA_VISIBLE_DEVICES` and select a valid local ordinal (usually device 0 in
  a one-GPU allocation). Do not unset it.
- `RuntimeError: CUDA out of memory`: rerun the same method with
  `--dft-engine cpu` or on a larger GPU. A smaller grid or basis changes the
  method, so run it only as a new, labeled calculation.
- When GPU4PySCF cannot run, `dft` stops with an error that suggests
  `--dft-engine cpu`; it does not switch by itself. The PyPI wheel is
  x86_64-only, so on aarch64 use `--dft-engine cpu` or build `gpu4pyscf` from
  source.
- `--func-basis` follows PySCF names. Test a basis name directly:
  `python -c "from pyscf import gto; print(len(gto.basis.load('def2-tzvpd', 'C')))"`.

## DFT//MLIP on the TS candidate

After `all --tsopt`, the structures to feed `dft` are
`<out_dir>/segments/seg_NN/{reactant,ts,product}.pdb` (plus `.cif` for mmCIF
input). With `result_mep` as the output directory of `all`:

```bash
LIGAND_CHARGE='SAM:1,GPP:-3'
FUNC_BASIS='wb97m-v/def2-tzvpd'
for state in reactant ts product; do
  pdb2reaction dft -i result_mep/segments/seg_01/${state}.pdb \
      -l "$LIGAND_CHARGE" --func-basis "$FUNC_BASIS" --dft-engine gpu -o dft_${state}
done
```

To take only the segment with the highest local barrier
(`rate_limiting_step.segment` in `summary.json`):

```bash
SUMMARY=result_mep/summary.json
RLS_SEG=$(python - "$SUMMARY" <<'PY'
import json, sys
with open(sys.argv[1], encoding="utf-8") as handle:
    print(int(json.load(handle)["rate_limiting_step"]["segment"]))
PY
)
printf -v RLS_DIR 'seg_%02d' "$RLS_SEG"
TS_FILE="result_mep/segments/${RLS_DIR}/ts.pdb"
test -f "$TS_FILE"
pdb2reaction dft -i "$TS_FILE" -l 'SAM:1,GPP:-3' \
    --func-basis 'wb97m-v/def2-tzvpd' --dft-engine gpu
```

Combine the energies with [`energy-diagram`](utilities.md#energy-diagram).

## Next step

- Install, and aarch64 handling:
  [backends.md](../pdb2reaction-install-backends/backends.md#dft-pyscf-gpu4pyscf).
- `-b dft`, `all --dft`, and GPU memory: [docs](../../docs/dft-backend.md).
- Geometries for the single points: [tsopt.md](tsopt.md), [irc.md](irc.md).
- Flags and defaults: `--help-advanced` and [SKILL.md](SKILL.md#where-flags-and-defaults-live).
