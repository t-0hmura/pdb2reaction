---
name: pdb2reaction-install-backends
description: Install recipes for pdb2reaction core + MLIP / DFT / xTB backends, with CUDA / PyTorch / e3nn / aarch64 quirks. TRIGGER on install / setup / `pip install` / `conda env` / `ImportError` / CUDA-version mismatch / "GPU not detected" / `huggingface` auth / e3nn version conflict questions. SKIP when pdb2reaction already imports cleanly and the user is invoking subcommands — the CLI skill covers usage.
---

# pdb2reaction — Install and Backends

## Overview

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

```
Need TS + IRC validation on a known organic + 1st-row metal cluster?
    └── start with UMA-s-1.1 (uma.md)
        └── if accuracy is borderline, also install MACE-OMOL-0 (mace.md, separate env)

Need a fast screen across many candidates?
    └── Orb-v3 (orb.md)

Working on small organic molecules, no metals?
    └── AIMNet2 (aimnet2.md) — limited element coverage, but light

Need DFT//MLIP refinement?
    └── add dft.md regardless of MLIP choice
```

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
      - torch
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
      - torch
      - pdb2reaction
# After conda env create, run inside the env:
#   pip uninstall -y fairchem-core      # remove UMA's e3nn pin
#   pip install mace-torch              # pulls the e3nn version MACE needs
```

`<cu_index>` is one of `cpu`, `cu118`, `cu121`, `cu124`, `cu126`, `cu128`, `cu129` — see
`env-cuda.md` for the driver version → index mapping.

## Verify the install

```bash
pdb2reaction --version
pdb2reaction --help                     # subcommand list
python -c "import pdb2reaction.defaults as d; print(sorted(n for n in dir(d) if not n.startswith('_')))"

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
| `import torch` fails with `libcudart.so.12 not found` | Torch wheel CUDA index mismatches the driver | `env-cuda.md` (driver → cu index table) |
| `e3nn` version conflict on `pip install` | UMA + MACE in the same env | Use a separate env for MACE (this file, `mace.md`) |
| `gpu4pyscf` import fails on aarch64 | `gpu4pyscf-cuda12x` is x86_64 only | `dft.md` — fall back to CPU PySCF |
| `huggingface_hub.errors.GatedRepoError` on UMA load | UMA model is gated, not authenticated | `uma.md` — `hf auth login` (legacy: `huggingface-cli login`) |
| `OSError: libcusolver.so.11 not found` | torch's bundled CUDA libs missing or shadowed | `env-cuda.md` — `LD_LIBRARY_PATH` order |
| `RuntimeError: CUDA out of memory` during freq | Hessian batch too large | reduce `hessian_calc_mode` batch in `pdb2reaction.defaults.UMA_CALC_KW` |

## Live source

`pyproject.toml` lists the canonical extras and version pins. To inspect
without opening the file:

```bash
python -c "import importlib.metadata as m; print(m.metadata('pdb2reaction').get_all('Provides-Extra'))"
python -c "import importlib.metadata as m; print(m.requires('pdb2reaction'))"
```