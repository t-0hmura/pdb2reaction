---
name: pdb2reaction-install-backends
description: "Install and environment setup for pdb2reaction: the core package, MLIP backends (UMA, ORB, MACE, AIMNet2), the optional DFT calculator (PySCF, GPU4PySCF), the xTB implicit-solvent correction, CUDA and PyTorch pairing, aarch64 limits, and probing an unknown compute environment (scheduler, GPU, CUDA, conda, modules). SKILL.md gives the install order, backend choice, verification commands, a custom ASE calculator via --calc-file, and a failure-to-fix table; backends.md holds the per-backend steps and the environment probes. TRIGGER on pip or conda install, ImportError, CUDA or driver mismatch, GPU not detected, Hugging Face authentication, MACE dependency conflicts, or when the compute environment is unknown. SKIP when pdb2reaction already imports cleanly and the question is about running subcommands (pdb2reaction-cli) or writing job scripts (pdb2reaction-hpc)."
---

# Install pdb2reaction

Check the machine → `torch==2.13.0` from the PyTorch index that matches the driver → `pip install pdb2reaction` (UMA's `fairchem-core` is a core dependency) → optional extras, or a separate env for MACE → `hf auth login` for UMA → `pdb2reaction --version`.

pdb2reaction needs a PyTorch wheel that matches the NVIDIA driver and at least
one MLIP backend; DFT and the `xtb` executable are optional. The install is done
when `pdb2reaction --version` prints the version and the smoke check in
[Verify the install](#verify-the-install) returns without an error.

## Install order

1. On a new or unknown host, run the probes in
   [backends.md](backends.md#probe-the-compute-environment) first.
2. Create a conda env with Python 3.12; pdb2reaction needs 3.11 or newer, and
   ORB needs 3.11 or 3.12. If you will use DMF (`--mep-mode dmf`), install
   `cyipopt` now; `pydmf` comes with pdb2reaction.
3. Install PyTorch. `nvidia-smi` shows `CUDA Version` at its top right, the
   newest CUDA the driver supports; choose a wheel at or below it (`cu126`,
   `cu130`, or `cu132`). `cu130` is the recommended choice. Details are in
   [CUDA and PyTorch](backends.md#cuda-and-pytorch).
4. Install pdb2reaction and headless Chrome for Plotly PNG export. The bundled
   `pysisyphus` and `thermoanalysis` come with it; do not install them
   separately.
5. Accept the UMA license on Hugging Face and log in (see [UMA](backends.md#uma)).
6. Add optional backends. MACE needs its own env
   ([MACE](backends.md#mace-separate-environment)); `xtb` comes from conda-forge
   ([xTB](backends.md#xtb-solvent-correction)).

```bash
conda create -n <YOUR_ENV> python=3.12 -y
conda activate <YOUR_ENV>
conda install -c conda-forge cyipopt -y            # only for --mep-mode dmf
pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cu130
pip install pdb2reaction                           # UMA via fairchem-core
plotly_get_chrome -y                               # headless Chrome; needs network
hf auth login                                      # Read token; once per machine and env
pip install --only-binary=dm-tree 'pdb2reaction[orb,aimnet,dft]'   # optional extras
```

| Extra | Pulls in | When you need it |
|---|---|---|
| (none) | UMA via `fairchem-core`, base deps | Default; `-b uma` works |
| `[orb]` | `orb-models` | `-b orb` |
| `[aimnet]` | `aimnet>=0.2.0` | `-b aimnet2` |
| `[dft]` | PySCF, CUDA 13 GPU4PySCF and CuPy on Linux x86_64 | `-b dft`, `--dft`, `pdb2reaction dft` with a `cu130` or `cu132` wheel |
| `[dft-cuda12]` | PySCF, CUDA 12 GPU4PySCF and CuPy | Same, with a `cu126` wheel |
| `[mcp]` | `mcp[cli]` | Running the MCP server |
| `[ci]`, `[dev]`, `[docs]` | test, pytest, and Sphinx dependencies | Contributing and building docs |

There is no `[mace]` extra, because `mace-torch` and `fairchem-core` need
different `e3nn` versions.

## Choose a backend

| `-b` | Default model | Install | Notes |
|---|---|---|---|
| `uma` (default) | `uma-s-1p2`; `uma-m-1p1` via `--backend-model` | Core dependency; gated weights need a Hugging Face login | fp32 by default; the only backend with several workers |
| `orb` | `orb_v3_conservative_omol` | `[orb]` extra | fp64 by default |
| `mace` | `MACE-OMOL-0` | Separate env | fp64 by default |
| `aimnet2` | `aimnet2` | `[aimnet]` extra | fp32 only; default model covers 14 elements, no Mg or first-row transition metals |

Select the backend with `-b` on any calculation command, for example
`pdb2reaction opt -i input.pdb -q 0 -b orb`. Check each candidate model's card
for supported elements, charge, multiplicity, and training domain, then compare
energies, forces, frequencies, runtime, and memory on a representative system.
Add the `[dft]` extra when you need DFT//MLIP single points. For an engine that
is not a built-in MLIP, use a custom calculator.

## Custom backend (--calc-file)

Write a Python file with a `get_calculator()` factory that returns an
[ASE](https://wiki.fysik.dtu.dk/ase/) Calculator (GFN-xTB, DFTB+, ORCA, Psi4,
and so on), and pass it with `--calc-file`. It selects the `custom` backend and
overrides `-b`.

```bash
# my_calc.py:
#   from ase.calculators.emt import EMT
#   def get_calculator(charge=0, spin=1, device="auto", **kwargs):
#       return EMT()              # swap for tblite.ase.TBLite(...) etc.
pdb2reaction sp -i model.xyz --calc-file my_calc.py -q 0 -m 1
```

`--calc-file` works on `sp`, `opt`, `tsopt`, `freq`, `irc`, `scan`, `scan2d`,
`scan3d`, `path-opt`, `path-search`, and `all`, which passes it to every stage.
Rename the factory with `--calc-file-func-name NAME`. Energies and forces follow
the ASE contract (eV, eV/Å). Hessians always use finite differences, and frozen
atoms (`--freeze-links`, `--freeze-atoms`) are honored. The full guide is the
"Custom backend" section of [MLIP Backends](../../docs/backends.md).

## Verify the install

```bash
pdb2reaction --version                 # installed version
pdb2reaction --help                    # lists the 18 subcommands
hf auth whoami                         # Hugging Face account used for UMA downloads
python -c "import torch; print('CUDA:', torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'N/A')"
python -c "from pdb2reaction.backends import create_calculator; create_calculator(backend='uma', charge=0, spin=1)"
python -c "import importlib.metadata as m; print(m.metadata('pdb2reaction').get_all('Provides-Extra'))"
python -c "import importlib.metadata as m; print(m.requires('pdb2reaction'))"
```

Success is a version string, your Hugging Face user name, `CUDA: True` with the
GPU name, and a smoke check that returns silently. Repeat the smoke check with
`backend='orb'`, `'aimnet2'`, or `'mace'` (MACE env only) for each backend you
installed. An `ImportError` points to that backend's section in
[backends.md](backends.md); a CUDA error points to
[CUDA and PyTorch](backends.md#cuda-and-pytorch).

## Common failure → fix

| Symptom | Likely cause | Fix |
|---|---|---|
| `import torch` fails with `libcudart.so.12 not found` | Incomplete or mixed wheel install, or a library-path collision; not by itself evidence against the driver | [CUDA and PyTorch](backends.md#cuda-and-pytorch) diagnostics |
| `torch.cuda.is_available()` is `False` | CPU wheel, hidden GPU, driver and wheel mismatch, or mixed libraries | Run `python -m torch.utils.collect_env` and `python -m pip check` before changing the wheel |
| `OSError: libcusolver.so.11 not found` | CUDA wheel dependency missing or shadowed by `LD_LIBRARY_PATH` | Compare with a clean environment ([CUDA and PyTorch](backends.md#cuda-and-pytorch)) |
| `e3nn` conflict on `pip install` | UMA and MACE in one env | Separate env for MACE |
| `huggingface_hub.errors.GatedRepoError`, `401`, or `403` on UMA load | License not accepted or not logged in | Accept the license on `facebook/UMA`, run `hf auth login`, confirm with `hf auth whoami`; on HPC, compute nodes must be able to write the Hugging Face cache |
| `gpu4pyscf` import fails on aarch64 | The PyPI wheel is x86_64 only | Build GPU4PySCF from source or run with `--dft-engine cpu` |
| `DMF mode (--mep-mode dmf) requires ase, cyipopt, and pydmf>=1.2` or `No module named 'dmf'` | `cyipopt` missing, or `pydmf` lacks its GPU part | `conda install -c conda-forge cyipopt`; if the import still fails, `pip install 'pydmf[torch]'` |
| Plot export fails, or `path-search` / `all` finish without the diagram PNG | No headless Chrome for Plotly | `plotly_get_chrome -y` |
| `RuntimeError: CUDA out of memory` during `freq` | Hessian too large for VRAM | Keep the default `FiniteDifference` Hessian, use a smaller model with `--backend-model`, or a larger GPU ([freq](../pdb2reaction-cli/freq.md)) |

## When the environment is unknown

Other pdb2reaction skills assume placeholders such as `<YOUR_QUEUE>`, `<NCPU>`,
`<NGPU>`, `<CUDA_MODULE>`, and `<YOUR_ENV>` are already known. When they are
not, for example on the first run on a new host or for an agent without prior
context, run the probes in
[Probe the compute environment](backends.md#probe-the-compute-environment)
and fill the placeholders from the output.

## Next step

- Run commands: [pdb2reaction-cli](../pdb2reaction-cli/SKILL.md).
- Write PBS or SLURM job scripts: [pdb2reaction-hpc](../pdb2reaction-hpc/SKILL.md).
- Per-backend steps, CUDA diagnostics, and environment probes: [backends.md](backends.md).
- Docs: [Installation](../../docs/installation.md) and [MLIP Backends](../../docs/backends.md).
