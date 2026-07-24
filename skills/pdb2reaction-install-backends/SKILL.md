---
name: pdb2reaction-install-backends
description: Install recipes for pdb2reaction core + MLIP / DFT / xTB backends, with CUDA / PyTorch / e3nn / aarch64 quirks. TRIGGER on install / setup / `pip install` / `conda env` / `ImportError` / CUDA-version mismatch / "GPU not detected" / `huggingface` auth / e3nn version conflict questions. SKIP when pdb2reaction already imports cleanly and the user is invoking subcommands — the CLI skill covers usage.
---

# pdb2reaction — Install and Backends

## Purpose

`pdb2reaction` is a Python package that depends on:

- a recent **PyTorch** wheel matching your CUDA driver,
- one or more **MLIP backends** (UMA / Orb / MACE / AIMNet2),
- optional **PySCF / GPU4PySCF** for DFT single points,
- optional **xtb** for ALPB solvent corrections.

Bundled and installed automatically with the package: `pysisyphus` (a
GPU-tensor fork), `thermoanalysis`. They must not be installed separately.

Files in this skill directory:

| File | When to consult |
|---|---|
| `SKILL.md` (this file) | Orientation, decision tree, install order |
| `env-cuda.md` | After confirming you have GPU + driver, before installing torch |
| `core.md` | Installing `pdb2reaction` itself |
| `uma.md` | Installing UMA (default backend; HuggingFace auth required) |
| `orb.md` | Installing Orb-v3 |
| `mace.md` | Installing MACE — **separate environment** required |
| `aimnet2.md` | Installing AIMNet2 |
| `dft.md` | PySCF + GPU4PySCF (handled separately per `[dft]` extra) |
| `xtb.md` | xTB ALPB solvent correction layer |

## Install order

1. **Check the env** — see `pdb2reaction-env-detect/SKILL.md` to discover
   driver / scheduler / CUDA / conda before doing anything.
2. **CUDA + torch** — `env-cuda.md` decides which torch wheel to pull.
3. **`pdb2reaction` core** — `core.md`.
4. **At least one MLIP backend** — start with UMA (`uma.md`); add others
   only as you need them.
5. **DFT (optional)** — `dft.md`. Skip if you only need MLIP energies.
6. **xtb (optional)** — `xtb.md`. Skip unless you need implicit solvent.

## Decision tree: which backend?

| Goal | Recommendation |
|---|---|
| TS + IRC on known organic + 1st-row metal cluster | start with **UMA-s-1.2** (`uma.md`); if accuracy is borderline, add **MACE-OMOL-0** in a separate env (`mace.md`) |
| Fast screen across many candidates | **Orb-v3** (`orb.md`) |
| Small organics, no metals | **AIMNet2** (`aimnet2.md`) — limited element coverage, light |
| Couple a non-MLIP engine (GFN-xTB / DFTB+ / ORCA / any ASE calc) | **`--calc-file my_calc.py`** — custom backend (see below) |
| DFT//MLIP single-point energies | add **`dft.md`** regardless of MLIP choice |

## Custom backend — any ASE Calculator (`--calc-file`)

To drive energies/gradients with an engine that is not a built-in MLIP
(GFN-xTB, DFTB+, ORCA, Psi4, …), write a Python file exposing a
`get_calculator()` factory that returns an [ASE](https://wiki.fysik.dtu.dk/ase/)
Calculator, and pass it with `--calc-file` (selects the `custom` backend,
overriding `-b`):

```bash
# my_calc.py:
#   from ase.calculators.emt import EMT
#   def get_calculator(charge=0, spin=1, device="auto", **kwargs):
#       return EMT()              # swap for tblite.ase.TBLite(...) etc.
pdb2reaction sp -i model.xyz --calc-file my_calc.py -q 0 -m 1
```

- Works on every listed calculator-backed geometry subcommand (`sp` / `opt` /
  `tsopt` / `freq` / `irc` / `scan` / `scan2d` / `scan3d` / `path-opt` / `path-search`) and
  **the `all` pipeline**
  (forwarded to each child stage). Rename the factory with `--calc-factory NAME`.
- Energy/forces follow the ASE eV / eV·Å⁻¹ contract; Hessians use the
  finite-difference path. Frozen atoms (`--freeze-links` / `--freeze-atoms`) are
  honored. Full guide: `docs/backends.md` (Custom backend section).

## Why two envs for MACE

`fairchem-core` (UMA) and `mace-torch` pin **different `e3nn` versions**
which cannot coexist. Solution: keep UMA in your default env (`<env_a>`)
and put MACE in a second env (`<env_b>`). Other backends (Orb, AIMNet2,
DFT, xTB) can sit in either.

`pdb2reaction` itself is the same code in both envs; only the calculator
plugin set differs.

## Conda env templates

Replace `<...>` with the values you discovered in `env-detect`. The
templates assume `python=3.11`; `pdb2reaction` requires Python ≥ 3.11.

`env_pdb2reaction.yml` (UMA / Orb / AIMNet2 / DFT / xTB):

```yaml
name: <your_env>
channels: [conda-forge, nvidia]
dependencies:
  - python=3.11
  - xtb                                # if you need xtb.md features (binary; called via subprocess)
  - pip
  - pip:
      - --extra-index-url https://download.pytorch.org/whl/<cu_index>
      - torch==2.8.0
      - pdb2reaction[orb,aimnet,dft]   # extras: see core.md / per-backend md
```

`env_pdb2reaction_mace.yml` (MACE only, separate env). `fairchem-core`
is a **core** dependency of `pdb2reaction` (not an extra), so it gets
pulled in by `pip install pdb2reaction` and must be removed *after*
install — see `mace.md`:

```yaml
name: <your_mace_env>
channels: [conda-forge]
dependencies:
  - python=3.11
  - pip
  - pip:
      - --extra-index-url https://download.pytorch.org/whl/<cu_index>
      - torch==2.8.0
      - pdb2reaction
# After conda env create, run inside the env:
#   pip uninstall -y fairchem-core      # remove UMA's e3nn pin
#   pip install mace-torch              # pulls the e3nn version MACE needs
```

`<cu_index>` is one of `cpu`, `cu126`, `cu128`, `cu129` — the indexes in
PyTorch's official 2.8.0 install matrix, which is the version family this
release pins. See `env-cuda.md`; do not infer the correct index only from the
"CUDA Version" banner printed by `nvidia-smi`.

## Verify the install

```bash
pdb2reaction --version
pdb2reaction --help                     # subcommand list
python -c "import pdb2reaction.core.defaults as d; print(sorted(n for n in dir(d) if not n.startswith('_')))"

# backend smoke checks (only those you installed)
python -c "from pdb2reaction.backends import create_calculator; create_calculator(backend='uma',  charge=0, spin=1)"
python -c "from pdb2reaction.backends import create_calculator; create_calculator(backend='orb',  charge=0, spin=1)"
python -c "from pdb2reaction.backends import create_calculator; create_calculator(backend='mace', charge=0, spin=1)"   # only in MACE env
python -c "from pdb2reaction.backends import create_calculator; create_calculator(backend='aimnet2', charge=0, spin=1)"
```

If any of these fails with `ImportError`, the corresponding backend's
`md` page lists missing dependencies. If a CUDA-related error appears,
go back to `env-cuda.md`.

## Common failure → fix

| Symptom | Likely cause | Where to look |
|---|---|---|
| `import torch` fails with `libcudart.so.12 not found` | Incomplete/mixed wheel install or a library-path collision (not enough evidence by itself to blame the driver) | `env-cuda.md` diagnostics |
| `e3nn` version conflict on `pip install` | UMA + MACE in the same env | Use a separate env for MACE (this file, `mace.md`) |
| `gpu4pyscf` import fails on aarch64 | PyPI wheel is x86_64-only | `dft.md` — build `gpu4pyscf` from source (https://github.com/pyscf/gpu4pyscf) or use CPU PySCF |
| `huggingface_hub.errors.GatedRepoError` on UMA load | UMA model is gated, not authenticated | `uma.md` — `hf auth login` (older `huggingface_hub < 0.34`: `huggingface-cli login`) |
| `OSError: libcusolver.so.11 not found` | CUDA wheel dependency missing or shadowed by environment libraries | `env-cuda.md` — compare a clean environment before editing library paths |
| `RuntimeError: CUDA out of memory` during freq | Hessian evaluation too large for VRAM | switch to `--hessian-calc-mode FiniteDifference` (see `pdb2reaction-cli/freq.md`) or reduce the active region |

## See also
`pyproject.toml` lists the canonical extras and version pins. To inspect
without opening the file:

```bash
python -c "import importlib.metadata as m; print(m.metadata('pdb2reaction').get_all('Provides-Extra'))"
python -c "import importlib.metadata as m; print(m.requires('pdb2reaction'))"
```
