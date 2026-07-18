# Reproducibility and determinism

MLIP inference on a GPU may be **non-bit-reproducible by default**: parallel
reductions (atomic adds, scatter operations) accumulate in a
hardware-scheduling-dependent order, so two runs with identical inputs differ
at the floating-point ULP level. Small drift was observed in the project's UMA
smoke benchmark, but its size does not guarantee the same for another backend, model,
GPU, software stack, or long optimizer trajectory. Assess scientific
reproducibility using chemically relevant observables, not file identity alone.

If you need same-stack repeatability (e.g. golden-file regression tests), use
the `--deterministic` flag and keep the backend/model, package versions, GPU,
inputs, and command fixed.

## `--deterministic`

`--deterministic` is accepted by every compute subcommand
(`opt`, `tsopt`, `freq`, `irc`, `scan`, `scan2d`, `scan3d`, `path-opt`,
`path-search`, `all`, `sp`). It turns on `torch.use_deterministic_algorithms`
plus an `index_reduce_` shim. This requests strict determinism from PyTorch; it
cannot control an arbitrary custom ASE calculator or every third-party custom
kernel.

```bash
pdb2reaction opt -i input.pdb -q 0 --deterministic
pdb2reaction all -i r.pdb p.pdb -q -1 --tsopt --deterministic
```

- It is **process-global**: setting it on `all` applies to every in-process
  stage; you do not pass it per stage.
- It is **slower**: deterministic scatter/reduce kernels have lower throughput
  than the default ones. Use it only when you need strict same-stack
  repeatability.
- It **fails loudly for operations PyTorch identifies**: if the current build cannot provide a
  deterministic kernel for an operation in your run, the command raises rather
  than silently producing non-reproducible output.
- The environment variable `PDB2REACTION_STRICT_DETERMINISTIC=1` is the equivalent
  entry point for CI or the direct Python API (`create_calculator`).

### Backend scope

| Backend | `--deterministic` |
|---|---|
| `uma` | Same-stack end-to-end smoke gate compares output coordinates from two runs |
| `orb` / `mace` | PyTorch strict mode is enabled, but repeatability must be smoke-tested for the installed third-party backend version |
| `aimnet2` | **not supported — rejected** (see below) |
| `custom` | The user-supplied ASE calculator owns determinism; the flag cannot guarantee it |

## Precision and reproducibility

Changing `--precision` changes numerical error and may change an optimizer
trajectory, but fp64 by itself does **not** make a GPU run bit-identical:
reduction-order non-determinism is independent of precision. Strict mode
targets same-stack repeatability, not cross-version or cross-hardware identity.

Model precision and Hessian dtype are separate knobs. Hessian evaluation
defaults to fp64, while `calc.hessian_double: false` can explicitly request the
model's native dtype. Passing `--precision fp64` forces the Hessian back to fp64
so optimizer linear algebra cannot silently run below the model precision.

### Choosing precision by backend and purpose

`--precision` selects the floating-point precision of MLIP inference
(`fp32` | `fp64`, case-insensitive). It is backend-agnostic — the CLI routes the
value into each backend's native key (UMA `precision`, ORB `precision`, MACE
`default_dtype`; for `aimnet2`, `fp32` is a no-op and `fp64` is rejected because
its model inputs are cast to float32 upstream). Left unset, UMA runs fp32 (the
screening-speed baseline) while ORB and MACE run fp64 — see
[Backends → Precision](backends.md#precision). Which value to pick depends on the
purpose. GPU class affects performance, not which result qualifies as a
validated TS:

| Purpose | Recommended | Why |
| --- | --- | --- |
| Routine run | Leave unset (`auto`) | Keeps tested defaults: UMA/AIMNet2 fp32, ORB/MACE fp64. |
| Speed screening | Explicit `--precision fp32` only when needed | This downgrades ORB/MACE; do not trust those finite-difference Hessians as final. |
| Final TS/Hessian | Keep ORB/MACE fp64; consider UMA fp64 if noise matters | Validate with an independent frequency calculation and IRC regardless of precision. |

fp64 has a non-trivial effect on TS optimization and Hessians for the
OMol-trained UMA backend. Measure the hardware cost and document the production
setting. `--deterministic` is a separate same-stack repeatability control and
does not make a lower-precision PES more accurate.

## AIMNet2 limitations

AIMNet2 does not support these features:

- **`--precision fp64`** — AIMNet2's model inputs are cast to float32 upstream,
  so an "fp64" run would not actually be fp64.
- **`--deterministic`** — AIMNet2 computes forces via a custom CUDA kernel
  that lies outside `torch.use_deterministic_algorithms` control, so its forces
  are not bit-reproducible (energy is). PyTorch's deterministic mode neither
  detects nor controls the custom op, so the limitation is reported explicitly.

For strict same-stack repeatability, use a supported backend with
`--deterministic` and verify a two-run smoke test on the installed stack. UMA
has the project's end-to-end gate; ORB and MACE still require a gate for the
installed third-party version.
