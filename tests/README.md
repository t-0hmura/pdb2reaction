# Tests Layout

`pdb2reaction` uses `tests/` as the single test root.

## Unit and CI Tests

CPU-only tests live directly under `tests/*.py`.

```bash
pytest tests/ -v --tb=short -x
```

These tests cover CLI parsing, helper contracts, logging summaries, geometry
regressions. They are suitable for CI.

## Manual Smoke Tests

GPU smoke fixtures and commands live in `tests/smoke/`.

```bash
cp -a tests/smoke /path/to/writable-scratch/p2r-smoke
cd /path/to/writable-scratch/p2r-smoke
bash run.sh
```

`tests/smoke/run.sh` assumes the caller has already configured the Python
environment and CUDA runtime. It verifies that distribution metadata, the
imported module, and the module CLI report one consistent version, without
pinning a release number. Every CLI case is executed through
`python -m pdb2reaction` from the same interpreter. Scheduler wrappers and
environment activation stay out of the repository.

The default strict lane numerically compares UMA and ORB analytical Hessians
with finite differences. Run the same required checker in the isolated MACE
and AIMNet2 environments:

```bash
bash run_backend_hessian.sh mace
bash run_backend_hessian.sh aimnet2
```
