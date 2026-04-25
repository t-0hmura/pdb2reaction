# Installing pdb2reaction itself (core.md)

`pdb2reaction` is a pure-Python package; no native C/C++ build is needed.
The bundled `pysisyphus` (GPU-tensor fork) and `thermoanalysis` install
automatically with `pip install pdb2reaction`.

## Prerequisites

- Python ≥ 3.10
- A working PyTorch install matching your CUDA driver — see `env-cuda.md`
- (For DFT) PySCF / GPU4PySCF — see `dft.md`
- (For ALPB) xtb — see `xtb.md`

## Install from PyPI (recommended)

```bash
conda activate <YOUR_ENV>
pip install pdb2reaction                         # core only (UMA + Orb)
pip install 'pdb2reaction[orb,aimnet,dft]'        # extras as needed
```

Available extras (canonical list lives in `pyproject.toml`):

| Extra | Pulls in | When you need it |
|---|---|---|
| (none) | UMA via `fairchem-core`, base deps | Default; UMA backend works out of the box |
| `[orb]` | `orb-models` | Using `-b orb` |
| `[aimnet]` | `aimnet>=0.1.0` | Using `-b aimnet2` |
| `[dft]` | `pyscf`, `gpu4pyscf-cuda12x` (x86_64), `cupy-cuda12x`, `cutensor-cu12`, `basis-set-exchange` | `pdb2reaction dft` subcommand |
| `[ci]` | CPU-only test deps (no GPU libs) | Running unit tests / docs builds |
| `[dev]` | `pytest` family | Contributing |

`[mace]` does **not** exist as an extra because MACE conflicts with
`fairchem-core`'s `e3nn` pin — install `mace-torch` manually in a
**separate environment** (see `mace.md`).

To inspect the live extras list without opening `pyproject.toml`:

```bash
python -c "import importlib.metadata as m; print(m.metadata('pdb2reaction').get_all('Provides-Extra'))"
```

## Install from source (development)

```bash
git clone <repo-url-from-the-toolkit-README> pdb2reaction
cd pdb2reaction
pip install -e '.[orb,aimnet,dft]'
```

`pip install -e .` will pick up edits to the source tree without
re-installing. Useful when you are debugging the toolkit itself.

## Verify the install

```bash
pdb2reaction --version
pdb2reaction --help

# Sanity import
python -c "
import pdb2reaction
import pdb2reaction.defaults as d
print('version :', pdb2reaction.__version__)
print('defaults:', sorted(n for n in dir(d) if not n.startswith('_'))[:10], '...')
"
```

`pdb2reaction --help` should list ~17 subcommands (`all`, `extract`,
`path-search`, `path-opt`, `opt`, `tsopt`, `freq`, `irc`, `dft`, `scan`,
`scan2d`, `scan3d`, `trj2fig`, `energy-diagram`, `add-elem-info`,
`fix-altloc`, `bond-summary`).

## Where things live after install

```bash
python -c "import pdb2reaction, os; print(os.path.dirname(pdb2reaction.__file__))"
```

Inside that directory:

| File / dir | Purpose |
|---|---|
| `cli.py` | Click entry point |
| `defaults.py` | All default kwarg dicts (read with `import pdb2reaction.defaults`) |
| `backends/` | UMA / Orb / MACE / AIMNet2 calculator factories |
| `pysisyphus/` | Bundled GPU-tensor pysisyphus fork |
| `thermoanalysis/` | Bundled QRRHO thermochemistry |
| `extract.py`, `path_search.py`, `tsopt.py`, `irc.py`, `freq.py`, `dft.py`, `all.py` | Subcommand implementations |

## Upgrading

```bash
pip install --upgrade pdb2reaction
pdb2reaction --version                    # confirm new version
```

When upgrading **across minor versions**, also re-check:

- `python -c "import pdb2reaction.defaults as d; print(d.RSIRFO_KW)"` — keys may have moved
- `pdb2reaction <subcommand> --help` — flag set may have changed
- `summary.json` schema — see `pdb2reaction-workflows-output/SKILL.md`

## Uninstall / clean rebuild

```bash
pip uninstall pdb2reaction
# Or scrap the env entirely:
conda env remove -n <YOUR_ENV>
```

## See also

- `env-cuda.md` — torch / CUDA setup that must come first
- `uma.md`, `orb.md`, `mace.md`, `aimnet2.md` — backend-specific extras
- `dft.md` — `[dft]` install details and aarch64 fallback