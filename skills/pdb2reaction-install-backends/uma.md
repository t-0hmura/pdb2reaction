# UMA backend (uma.md)

UMA (**U**niversal **M**odel for **A**toms, Meta FAIR) is the default
backend for `pdb2reaction`. It covers the broadest element / chemistry
range of the four backends.

## Install

UMA pulls in via `fairchem-core`, which is a **core dependency**, so
you don't need an extras flag:

```bash
pip install pdb2reaction                              # fairchem-core comes along
```

Confirm:

```bash
python -c "import fairchem; print('fairchem :', fairchem.__version__)"
python -c "from pdb2reaction.backends import create_calculator; create_calculator(backend='uma', charge=0, spin=1)"
```

## HuggingFace authentication (required)

UMA model weights are gated on HuggingFace and need an authenticated
download:

```bash
pip install huggingface_hub[cli]
hf auth login                # paste a Read token from huggingface.co/settings/tokens
```

Requires `huggingface_hub >= 0.34`; on older versions use the legacy
`huggingface-cli login` (still works but is being deprecated).

The token is cached in `~/.cache/huggingface/`. Once it's there, future
runs (and PBS jobs) pick it up automatically.

If you hit `huggingface_hub.errors.GatedRepoError` or
`401 Client Error: Unauthorized`, re-run `hf auth login` and make sure
the token has access to the UMA model repos
(`facebook/UMA-S-1.1`, `facebook/UMA-M-1.1`).

## CLI usage

`uma` is the default — `pdb2reaction all -i ...` uses UMA-s-1.1 unless
overridden:

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    -b uma                       # explicit, identical to default
```

Available models (select with the `--backend-model` CLI flag, or via the
`calc.model:` field in `--config` YAML / `pdb2reaction.core.defaults.UMA_CALC_KW`).
Two equivalent notations are common:

| config string (`calc.model`) | paper notation | HuggingFace repo | Notes |
|---|---|---|---|
| `uma-s-1p1` (default) | UMA-s-1.1 | `facebook/UMA-S-1.1` | Smaller / faster, sufficient for most workflows |
| `uma-s-1p2` | UMA-s-1.2 | `facebook/UMA-S-1.2` | Successor checkpoint (still small) |
| `uma-m-1p1` | UMA-m-1.1 | `facebook/UMA-M-1.1` | Larger, slightly more accurate, ~3× slower |

`p` is the dot replacement used by fairchem-core's config parser
(`1p1` ↔ `1.1`). Pass the model name via `--backend-model` (e.g.
`--backend-model uma-s-1p2`) or the `calc.model` YAML key.

Inspect the full default kwarg dict:

```bash
python -c "import pdb2reaction.core.defaults as d; print(d.UMA_CALC_KW)"
```

## Backend-specific flags

UMA accepts these calculator kwargs (canonical list in
`backends/__init__.py:_BACKEND_ACCEPTED_KEYS['uma']`):

| Key | Purpose |
|---|---|
| `charge`, `spin` | Total charge and spin multiplicity |
| `device` | `'cuda'`, `'cpu'`, or `'auto'` |
| `model` | `'uma-s-1p1'`, `'uma-s-1p2'`, or `'uma-m-1p1'` |
| `task_name` | `'omol'` (default — organic molecules + 1st-row metals) |
| `freeze_atoms` | Indices of atoms held fixed (link-atom parents, frozen residues) |
| `hessian_calc_mode` | `'FiniteDifference'` (default) or `'Analytical'` |
| `return_partial_hessian`, `hessian_double` | Memory / numerical-precision toggles |
| `workers`, `workers_per_node` | Multi-GPU inference (uses Ray) |
| `max_neigh`, `radius` | Neighbor-list cutoffs |

These are passed through `--config` YAML or the appropriate flag of the
subcommand; see `pdb2reaction-cli/SKILL.md`.

## Multi-GPU inference (advanced)

Under heavy MEP search load you can shard inference across multiple GPUs.
Configure via a YAML config (`--config`):

```yaml
# multi_worker.yaml
calc:
  workers: 4
  workers_per_node: 4
```

```bash
pdb2reaction all -i ... -b uma --config multi_worker.yaml
```

This spawns a Ray worker pool. Limitations:

- `pdb2reaction` does **not** launch a cross-node Ray cluster by
  itself; for multi-node use, start the Ray cluster externally
  under your scheduler (see `pdb2reaction/docs/hpc-example.md` for
  a PBS template) and `pdb2reaction` will join it via `RAY_ADDRESS`.
  In a single-node configuration the Ray pool is started locally.
- All workers in a single-node pool must see the same GPUs
  (e.g. `CUDA_VISIBLE_DEVICES=0,1,2,3`).
- Adds overhead for small graphs — disable for systems below ~100 atoms.
- **Analytical Hessian unavailable**: when `workers > 1`, requesting
  `hessian_calc_mode=Analytical` raises a `RuntimeError`; use
  `FiniteDifference` (the default) or drop to `workers = 1`. See
  `pdb2reaction/docs/uma-pysis.md` for the full caveat.

## Known gotchas

| Symptom | Cause / fix |
|---|---|
| `e3nn` install conflict | UMA's `fairchem-core` pin clashes with `mace-torch`. Use a separate env for MACE (see `mace.md`). |
| `uma-m-1.1` runs out of VRAM during freq | `'FiniteDifference'` is already the default; if you've overridden it with `--hessian-calc-mode Analytical`, drop the flag or use `uma-s-1.1`. |
| First call is slow (10–30 s) | One-time model download + JIT compile. The cache lives at `~/.cache/huggingface/hub/`. |
| Multi-worker run crashes with `Ray actor died` | Mismatched CUDA versions across processes; fall back to single-worker. |

## See also

- `env-cuda.md` — torch + CUDA setup (UMA needs CUDA-enabled torch).
- `core.md` — `pdb2reaction` itself.
- `mace.md` — alternate backend, requires a **separate** env.
- [`pdb2reaction-cli/tsopt.md`](../pdb2reaction-cli/tsopt.md), [`pdb2reaction-cli/freq.md`](../pdb2reaction-cli/freq.md) — `--hessian-calc-mode` choices.