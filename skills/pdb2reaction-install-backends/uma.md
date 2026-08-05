# UMA backend (uma.md)

UMA (**U**niversal **M**odel for **A**toms, Meta FAIR) is the default
backend integration for `pdb2reaction`. Model-domain coverage and accuracy
must be checked for the chemistry being studied; the default status is not a
claim that it is always the most accurate or broadest choice.

## Install

UMA ships via `fairchem-core`, a **core dependency**, so
you don't need an extras flag:

```bash
pip install pdb2reaction                              # fairchem-core comes along
```

Confirm:

```bash
python -c "from importlib.metadata import version; print('fairchem-core:', version('fairchem-core'))"
python -c "from pdb2reaction.backends import create_calculator; create_calculator(backend='uma', charge=0, spin=1)"
```

## HuggingFace authentication (required)

Every UMA checkpoint lives in one **gated** HuggingFace repo,
`facebook/UMA`. Two steps unlock the download.

**1. Accept the license.** Open <https://huggingface.co/facebook/UMA>
and accept the FAIR Chemistry License v1 (approval is manual). Use the
same HF account that owns the token in step 2.

**2. Authenticate.**

```bash
pip install huggingface_hub[cli]
hf auth login                # paste a Read token from huggingface.co/settings/tokens
```

Requires `huggingface_hub >= 0.34`; on older versions use the legacy
`huggingface-cli login` (still works but is being deprecated).

The token is cached in `~/.cache/huggingface/`. Once it's there, future
runs (and PBS jobs) pick it up automatically.

If you hit `huggingface_hub.errors.GatedRepoError` or
`401 Client Error: Unauthorized`, check that the account behind your
token has been granted access on <https://huggingface.co/facebook/UMA>
(step 1), then re-run `hf auth login`. `facebook/UMA` is the only repo
that needs access: the individual models are checkpoint *files* inside
it, not separate repos.

## CLI usage

`uma` is the default — `pdb2reaction all -i ...` uses UMA-s-1.2 unless
overridden:

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    -b uma                       # explicit, identical to default
```

Release-tested model examples (select with the `--backend-model` CLI flag, or via the
`calc.model:` field in `--config` YAML / `pdb2reaction.core.defaults.UMA_CALC_KW`).
The CLI accepts a model string; the upstream `facebook/UMA` repository is the
authoritative and evolving checkpoint inventory. These tested examples are:

| config string (`calc.model`) | paper notation | checkpoint in `facebook/UMA` | Notes |
|---|---|---|---|
| `uma-s-1p1` | UMA-s-1.1 | `checkpoints/uma-s-1p1.pt` | Small 1.1 checkpoint |
| `uma-s-1p2` (default) | UMA-s-1.2 | `checkpoints/uma-s-1p2.pt` | Small 1.2 checkpoint selected by pdb2reaction |
| `uma-m-1p1` | UMA-m-1.1 | `checkpoints/uma-m-1p1.pt` | Medium 1.1 checkpoint; benchmark cost and result quality for the target system |

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
| `task_name` | `'omol'` (default — the OMol25 molecular/polymer domain, which spans organic/inorganic molecules, transition-metal complexes, and electrolytes; validate applicability to the target system) |
| `freeze_atoms` | Indices of atoms held fixed (cap-atom parents, frozen residues) |
| `hessian_calc_mode` | `'FiniteDifference'` (default) or `'Analytical'` |
| `return_partial_hessian`, `hessian_double` | Memory / numerical-precision toggles |
| `precision` | `'fp32'` (UMA's default — the fairchem baseline) or `'fp64'` for full-precision base inference (needs fairchem-core ≥ 2.0). `--precision fp64` also forces a fp64 Hessian, and can change TSopt / Hessian numerics non-trivially |
| `workers`, `workers_per_node` | Multi-GPU inference (uses Ray) |
| `max_neigh`, `radius` | Neighbor-list cutoffs |

These are passed through `--config` YAML or the appropriate flag of the
subcommand; see `pdb2reaction-cli/SKILL.md`.

## Multi-worker inference (advanced)

Under heavy MEP search load, fairchem can run a parallel predictor worker
pool. GPU/process placement depends on the Ray/scheduler setup and must be
validated on the target cluster.
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
  under your scheduler (see `docs/hpc-example.md` for
  a PBS template) and `pdb2reaction` will join it via `RAY_ADDRESS`.
  In a single-node configuration the Ray pool is started locally.
- All workers in a single-node pool must see the same GPUs
  (e.g. `CUDA_VISIBLE_DEVICES=0,1,2,3`).
- Adds process/communication overhead; benchmark a representative task before
  assuming more workers reduce wall time.
- **Analytical Hessian unavailable**: when `workers > 1`, requesting
  `hessian_calc_mode=Analytical` raises `BackendError` (a `RuntimeError`
  subclass) rather than changing the explicitly requested method; use
  `FiniteDifference` (the default) or drop to `workers = 1`. See
  `docs/uma-pysis.md` for the full caveat.

## Known gotchas

| Symptom | Cause / fix |
|---|---|
| `e3nn` install conflict | UMA's `fairchem-core` pin clashes with `mace-torch`. Use a separate env for MACE (see `mace.md`). |
| Frequency calculation runs out of VRAM | Compare Hessian modes and compatible model sizes on a representative pilot, reduce the Hessian target, or move Hessian assembly to CPU. |
| First call is slower than later calls | Check whether checkpoint download/cache population or model initialization dominates. The Hugging Face cache normally lives at `~/.cache/huggingface/hub/`. |
| Multi-worker run crashes with `Ray actor died` | The message is non-specific. Capture the actor traceback and scheduler/Ray logs, verify identical environments and visible devices on every worker, then reproduce with one worker before changing CUDA packages. |

## See also

- `env-cuda.md` — torch + CUDA setup (UMA needs CUDA-enabled torch).
- `core.md` — `pdb2reaction` itself.
- `mace.md` — alternate backend, requires a **separate** env.
- [`pdb2reaction-cli/tsopt.md`](../pdb2reaction-cli/tsopt.md), [`pdb2reaction-cli/freq.md`](../pdb2reaction-cli/freq.md) — `--hessian-calc-mode` choices.
