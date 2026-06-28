# Installing pdb2reaction itself (core.md)

`pdb2reaction` is a pure-Python package; no native C/C++ build is needed.
The bundled `pysisyphus` (GPU-tensor fork) and `thermoanalysis` install
automatically with `pip install pdb2reaction`.

## Prerequisites

- Python ≥ 3.11
- A working PyTorch install matching your CUDA driver — see `env-cuda.md`
- (For DFT) PySCF / GPU4PySCF — see `dft.md`
- (For ALPB) xtb — see `xtb.md`

## Install from PyPI (recommended)

```bash
conda activate <YOUR_ENV>
pip install pdb2reaction                         # core (UMA via fairchem-core); Orb requires the [orb] extra
pip install 'pdb2reaction[orb,aimnet,dft]'        # extras as needed
```

Available extras (canonical list lives in `pyproject.toml`):

| Extra | Pulls in | When you need it |
|---|---|---|
| (none) | UMA via `fairchem-core`, base deps | Default; UMA backend works out of the box |
| `[orb]` | `orb-models` | Using `-b orb` |
| `[aimnet]` | `aimnet>=0.2.0` | Using `-b aimnet2` |
| `[dft]` | `pyscf`, `gpu4pyscf-cuda12x` (x86_64), `cupy-cuda12x`, `basis-set-exchange` | `pdb2reaction dft` subcommand |
| `[mcp]` | `mcp[cli]>=1.0` | Running the MCP server |
| `[ci]` | CPU-only test deps; CI installs the torch CPU wheel separately via the PyTorch CPU index | Running unit tests / docs builds |
| `[dev]` | `pytest` family | Contributing |
| `[docs]` | `sphinx`, `myst-parser`, `furo`, `sphinx-copybutton`, `sphinx-autobuild` | Building the Sphinx docs site |

`[mace]` does **not** exist as an extra because MACE conflicts with
`fairchem-core`'s `e3nn` pin — install `mace-torch` manually in a
**separate environment** (see `mace.md`).

To inspect the live extras list without opening `pyproject.toml`:

```bash
python -c "import importlib.metadata as m; print(m.metadata('pdb2reaction').get_all('Provides-Extra'))"
```

## Install from source (development)

```bash
git clone https://github.com/t-0hmura/pdb2reaction.git pdb2reaction
cd pdb2reaction
pip install -e '.[orb,aimnet,dft]'
```

## Verify the install

```bash
pdb2reaction --version
pdb2reaction --help

# Sanity import
python -c "
import pdb2reaction
import pdb2reaction.core.defaults as d
print('version :', pdb2reaction.__version__)
print('defaults:', sorted(n for n in dir(d) if not n.startswith('_'))[:10], '...')
"
```

`pdb2reaction --help` should list ~18 subcommands (`all`, `extract`,
`path-search`, `path-opt`, `opt`, `tsopt`, `freq`, `irc`, `dft`, `sp`, `scan`,
`scan2d`, `scan3d`, `trj2fig`, `energy-diagram`, `add-elem-info`,
`fix-altloc`, `bond-summary`).

## Where things live after install

```bash
python -c "import pdb2reaction, os; print(os.path.dirname(pdb2reaction.__file__))"
```

Inside that directory:

| File / dir | Purpose |
|---|---|
| `cli/app.py` | Click entry point |
| `core/defaults.py` | All default kwarg dicts (read with `import pdb2reaction.core.defaults`) |
| `backends/` | UMA / Orb / MACE / AIMNet2 calculator factories |
| `workflows/{extract,path_search,tsopt,irc,freq,dft,all}.py` | Subcommand implementations |

`pysisyphus/` and `thermoanalysis/` are **not** inside `pdb2reaction/` — they
install as separate top-level packages (siblings of `pdb2reaction/` in
`site-packages/`).

## Upgrading

```bash
pip install --upgrade pdb2reaction
pdb2reaction --version                    # confirm new version
```

When upgrading **across minor versions**, also re-check:

- `python -c "import pdb2reaction.core.defaults as d; print(d.RSIRFO_KW)"` — keys may have moved
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