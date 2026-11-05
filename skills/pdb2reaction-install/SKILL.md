---
name: pdb2reaction-install
description: "Install and environment setup for pdb2reaction: the core package, MLIP backends (UMA, ORB, MACE, AIMNet2), the optional DFT calculator (PySCF, GPU4PySCF), the xTB implicit-solvent correction, CUDA and PyTorch pairing, aarch64 limits, and probing an unknown compute environment (scheduler, GPU, CUDA, conda, modules). SKILL.md gives the install order, backend choice, verification commands, a custom ASE calculator via --calc-file, and a failure-to-fix table; backends.md holds the per-backend steps and the environment probes. TRIGGER on pip or conda install, ImportError, CUDA or driver mismatch, GPU not detected, Hugging Face authentication, MACE dependency conflicts, a `Missing system prerequisite:` message (compiler, Chrome libraries), or when the compute environment is unknown. SKIP when pdb2reaction already imports cleanly and the question is about running subcommands (pdb2reaction-cli) or writing job scripts (pdb2reaction-hpc)."
---

# Install pdb2reaction

Check the machine → `torch==2.13.0` from the PyTorch index that matches the driver → `pip install pdb2reaction` (UMA's `fairchem-core` is a core dependency) → optional extras, or a separate env for ORB or MACE → `hf auth login` for UMA → `pdb2reaction --version`.

pdb2reaction needs a PyTorch wheel that matches the NVIDIA driver and at least
one MLIP backend; DFT and the `xtb` executable are optional. The install is done
when `pdb2reaction --version` prints the version and the smoke check in
[Verify the install](#verify-the-install) returns without an error.

## Install order

1. On a new or unknown host, run the probes in
   [backends.md](backends.md#probe-the-compute-environment) first.
2. Create a conda env with Python 3.12 from conda-forge; pdb2reaction needs
   Python 3.12 (recommended) or 3.13. If you will use DMF
   (`--mep-mode dmf`) or `--solvent`, create it with `cyipopt`, `xtb`, and
   PyTorch from conda-forge instead (below) and skip step 3; `pydmf` comes with
   pdb2reaction.
3. Install PyTorch. `nvidia-smi` shows `CUDA Version` at its top right, the
   newest CUDA the driver supports; choose a wheel at or below it (`cu126`,
   `cu130`, or `cu132`). `cu130` is the recommended choice. Details are in
   [CUDA and PyTorch](backends.md#cuda-and-pytorch).
4. Install pdb2reaction and headless Chrome for Plotly image export (PNG, JPEG,
   SVG, PDF). The bundled `pysisyphus` and `thermoanalysis` come with it; do not
   install them separately.
5. Accept the UMA license on Hugging Face and log in (see [UMA](backends.md#uma)).
6. Add optional backends. ORB and MACE each need their own env
   ([ORB](backends.md#orb-separate-environment),
   [MACE](backends.md#mace-separate-environment)); `xtb` comes from conda-forge
   ([xTB](backends.md#xtb-solvent-correction)).

```bash
conda create -n <YOUR_ENV> -c conda-forge python=3.12 -y
conda activate <YOUR_ENV>
pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cu130
pip install pdb2reaction                           # UMA via fairchem-core
plotly_get_chrome -y                               # headless Chrome; needs network
hf auth login                                      # Read token; once per machine (shared by all envs)
pip install 'pdb2reaction[aimnet,dft]'       # optional extras
```

For DMF or `--solvent`, replace the first three lines with the following
(`cuda129` instead of `cuda130` when the driver is below CUDA 13.0):

```bash
conda create -n <YOUR_ENV> -c conda-forge python=3.12 cyipopt xtb "pytorch=2.13.0=cuda130*" -y
conda activate <YOUR_ENV>
```

| Extra | Pulls in | When you need it |
|---|---|---|
| (none) | UMA via `fairchem-core`, base deps | Default; `-b uma` works |
| `[aimnet]` | `aimnet>=0.2.0` | `-b aimnet2` |
| `[dft]` | PySCF, CUDA 13 GPU4PySCF and CuPy on Linux x86_64 | `-b dft`, `--dft`, `pdb2reaction dft` with a `cu130` or `cu132` wheel or the conda-forge `cuda130` build |
| `[dft-cuda12]` | PySCF, CUDA 12 GPU4PySCF and CuPy | Same, with a `cu126` wheel or the `cuda129` build |
| `[ci]`, `[dev]`, `[docs]` | test, pytest, and Sphinx dependencies | Contributing and building docs |

Install ORB and MACE directly in separate environments (see [Choose a backend](#choose-a-backend)).

## System packages

`pip install pdb2reaction` does not install a compiler or the libraries of Chrome, and a bare
container or cluster node often lacks both. On a GPU, UMA compiles code at run
time for speed and needs a C compiler (`gcc`) even without `--compile`; with
`--compile`, UMA and ORB need a C compiler on a GPU and a C++ compiler (`g++`) on a CPU.
Use the existing compiler; check with
`${CC:-gcc} --version` (`${CXX:-g++} --version` for C++).
If gcc or CUDA is missing, load them from the system (on HPC,
`module load gcc <CUDA_MODULE>`). If the system has no C compiler or CUDA, or
you do not want to use them, install both into the env instead:

```bash
conda install -c conda-forge -y cuda-nvcc
conda activate <YOUR_ENV>
```

Then rerun the command.

Image export (PNG, JPEG, SVG, PDF) needs the operating-system libraries of the Chrome that
`plotly_get_chrome -y` downloads. On HPC, ask the administrator for missing Chrome libraries.
The Chrome libraries are, for Ubuntu/Debian (Ubuntu 22.04 and
Debian 12: `libasound2` instead of `libasound2t64`):

```bash
sudo apt-get update
sudo apt-get install -y ca-certificates fonts-liberation libasound2t64 libatk-bridge2.0-0 \
  libatk1.0-0 libatspi2.0-0 libcairo2 libcups2 libcurl4 libdbus-1-3 libexpat1 libgbm1 \
  libglib2.0-0 libgtk-3-0 libnspr4 libnss3 libpango-1.0-0 libudev1 libvulkan1 libx11-6 libxcb1 \
  libxcomposite1 libxdamage1 libxext6 libxfixes3 libxkbcommon0 libxrandr2 wget xdg-utils
```

and for Rocky/RHEL:

```bash
sudo dnf install -y alsa-lib at-spi2-atk at-spi2-core atk ca-certificates cairo cups-libs \
  dbus-libs expat glib2 gtk3 libX11 libXcomposite libXdamage libXext libXfixes libXrandr \
  libgcc liberation-fonts libxcb libxkbcommon mesa-libgbm nspr nss nss-util pango systemd-libs \
  vulkan-loader wget xdg-utils
```

A missing package stops the command with a `Missing system prerequisite:` message
that lists the same commands. Docs: [System prerequisites](../../docs/installation.md#system-prerequisites).

## Choose a backend

| `-b` | Default model | Install | Notes |
|---|---|---|---|
| `uma` (default) | `uma-s-1p2p1`; `uma-s-1p2p1-ft-d4rev`, `uma-s-1p2`, `uma-m-1p1` via `--backend-model` | Core dependency; gated weights need a Hugging Face login | fp32 by default; the only backend with several workers |
| `orb` | `orbmol-v2`; `orb_v3_conservative_omol` via `--backend-model` | `orb-models` in a separate env | fp64 by default |
| `mace` | `MACE-OMOL-0`; `MACE-POLAR-1-S`/`-M`/`-L` (`mace-torch[polar]` >= 0.3.17) via `--backend-model` | Separate env | fp64 by default |
| `aimnet2` | `aimnet2`; `aimnet2-2025`, `aimnet2-nse` (open-shell), `aimnet2-pd` via `--backend-model` | `[aimnet]` extra | fp32 only; default model covers 14 elements, no Mg or first-row transition metals |

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
Use a different factory name with `--calc-file-func-name NAME`. Energies and forces follow
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
GPU name, and a smoke check that returns without an error. Repeat the smoke check with
`backend='orb'`, `'aimnet2'`, or `'mace'` (MACE env only) for each backend you
installed. An `ImportError` points to that backend's section in
[backends.md](backends.md); a CUDA error points to
[CUDA and PyTorch](backends.md#cuda-and-pytorch).

## Common failure → fix

| Symptom | Likely cause | Fix |
|---|---|---|
| `import torch` fails with `libcudart.so.<N> not found` (N is 12 for cu126, 13 for cu130) | Incomplete or mixed wheel install, or a library-path collision; not necessarily a driver problem | [CUDA and PyTorch](backends.md#cuda-and-pytorch) diagnostics |
| `torch.cuda.is_available()` is `False` | CPU wheel, hidden GPU, driver and wheel mismatch, or mixed libraries | Run `python -m torch.utils.collect_env` and `python -m pip check` before changing the wheel |
| `OSError: libcusolver.so.<N> not found` (N is 11 for cu126, 12 for cu130) | CUDA wheel dependency missing or shadowed by `LD_LIBRARY_PATH` | Compare with a clean environment ([CUDA and PyTorch](backends.md#cuda-and-pytorch)) |
| `e3nn` conflict on `pip install` | UMA and MACE in one env | Separate env for MACE |
| `huggingface_hub.errors.GatedRepoError`, `401`, or `403` on UMA load | License not accepted or not logged in | Accept the license on `facebook/UMA`, run `hf auth login`, confirm with `hf auth whoami`; on HPC, compute nodes must be able to write the model cache (`$FAIRCHEM_CACHE_DIR`, default `~/.cache/fairchem`) |
| `gpu4pyscf` import fails on aarch64 | The PyPI wheel is x86_64 only | Build GPU4PySCF from source or run with `--dft-engine cpu` |
| `DMF mode (--mep-mode dmf) requires ase, cyipopt, and pydmf>=1.2` or `No module named 'dmf'` | `cyipopt` missing, or `pydmf` lacks its GPU part | Recreate the env with `cyipopt` and PyTorch from conda-forge ([Install order](#install-order)); if the import still fails, `pip install 'pydmf[torch]'` |
| `Missing system prerequisite: C compiler (gcc).` or `C++ compiler (g++).`; a PyTorch error `Failed to find C compiler` or `No working C++ compiler found` | No compiler on `PATH`; UMA on a GPU, and UMA and ORB with `--compile`, compile code at run time for speed | Load gcc and CUDA (`module load`) or install `cuda-nvcc`: [System packages](#system-packages) |
| Plot export fails (`Missing system prerequisite: Chrome for PNG export.`), or `path-search` / `all` finish without the diagram PNG | No headless Chrome for Plotly, or an operating-system library it needs is missing | `plotly_get_chrome -y`, then the libraries in [System packages](#system-packages) |
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
