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

The strict lane numerically compares UMA and ORB analytical Hessians with
finite differences. ORB cannot share the UMA environment (see
`docs/installation.md`), so `run.sh` is run in two environments: the UMA
environment runs it without test53, test53b, and the `orb` argument of
`run_backend_hessian.sh`; the ORB environment runs test53, test53b, and
`bash run_backend_hessian.sh orb`.
Run the same required checker in the isolated MACE and AIMNet2 environments:

```bash
bash run_backend_hessian.sh mace
bash run_backend_hessian.sh aimnet2
```
