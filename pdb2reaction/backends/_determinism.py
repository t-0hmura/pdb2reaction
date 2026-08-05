"""Opt-in strict-determinism setup for the MLIP backends.

Activated by the ``--deterministic`` CLI flag (present on every compute
subcommand) or the ``PDB2REACTION_STRICT_DETERMINISTIC=1`` environment variable.
Process-global and idempotent: the first activation patches the
non-deterministic ``Tensor.index_reduce_(reduce="mean")`` op (which has no
deterministic CUDA kernel in torch <= 2.8 and crashes under
``use_deterministic_algorithms``) to a ``scatter_reduce`` detour, then enables
``torch.use_deterministic_algorithms(True)`` + cuDNN determinism + fixed seeds.

Design notes:
- **Patch first.** Native strict mode is *known* to crash on
  ``index_reduce_``; probing it first would only waste a guaranteed failure.
  The shim is applied before strict mode is turned on.
- **Fail loud.** If the patch target is gone under a torch upgrade, or strict
  mode rejects an op with no deterministic kernel, raise — never silently
  degrade to ``warn_only``.
- **Scope.** Requests same-software/hardware-stack determinism for operations
  controlled by PyTorch. It cannot guarantee cross-version/cross-hardware
  identity or control arbitrary third-party/custom kernels and calculators.
- **Cost.** Deterministic scatter/reduce kernels may be slower, and the runtime
  monkey-patch is version-sensitive. Default OFF.
"""
from __future__ import annotations

import os
import threading

import numpy as np

_DONE = False
_ORIG_INDEX_REDUCE = None
_SETUP_LOCK = threading.RLock()


def is_deterministic_requested() -> bool:
    """True when ``PDB2REACTION_STRICT_DETERMINISTIC=1`` (the env-var entry point used by
    CI and the direct Python API; the CLI uses ``--deterministic``)."""
    return os.environ.get("PDB2REACTION_STRICT_DETERMINISTIC") == "1"


def is_deterministic_active() -> bool:
    """True once strict-deterministic mode has actually been applied this process
    (set by ``setup_deterministic``; covers both the --deterministic flag and the
    env-var path)."""
    return _DONE


def _index_reduce_mean_deterministic(self, dim, index, source, reduce, include_self=True):
    """Drop-in for ``Tensor.index_reduce_`` routing mean reductions through
    ``scatter_reduce`` (which has a deterministic CUDA kernel)."""
    if reduce != "mean":
        return _ORIG_INDEX_REDUCE(self, dim, index, source, reduce=reduce, include_self=include_self)
    if index.dim() == 1 and source.dim() > 1:
        normalized_dim = int(dim) % source.dim()
        index_shape = [1] * source.dim()
        index_shape[normalized_dim] = index.numel()
        idx_exp = index.reshape(index_shape).expand_as(source)
    else:
        idx_exp = index
    result = self.scatter_reduce(dim, idx_exp, source, reduce="mean", include_self=include_self)
    self.copy_(result)
    return self


def setup_deterministic() -> None:
    """Enable strict-deterministic mode once (idempotent). Patch first, fail loud.

    Raises ``RuntimeError`` if determinism cannot be honoured on this torch
    build (patch target missing, shim self-check mismatch, or strict mode
    rejected) rather than degrading silently.
    """
    global _DONE, _ORIG_INDEX_REDUCE
    with _SETUP_LOCK:
        if _DONE:
            return
        import torch

        orig = getattr(torch.Tensor, "index_reduce_", None)
        if orig is None:
            raise RuntimeError(
                "--deterministic: torch.Tensor.index_reduce_ is missing on this "
                f"torch ({torch.__version__}); the determinism shim needs updating."
            )

        # Verify the replacement before changing any process-global state.
        checks = (
            (torch.zeros(3), 0, torch.tensor([0, 0, 1]), torch.tensor([1.0, 3.0, 5.0])),
            (torch.zeros(3, 2), 0, torch.tensor([0, 0, 1]), torch.arange(6.0).reshape(3, 2)),
            (torch.zeros(2, 3), 1, torch.tensor([0, 0, 1]), torch.arange(6.0).reshape(2, 3)),
        )
        for target, dim, idx, src in checks:
            for include_self in (False, True):
                ref = orig(
                    target.clone(),
                    dim,
                    idx,
                    src,
                    reduce="mean",
                    include_self=include_self,
                )
                got = _index_reduce_mean_deterministic(
                    target.clone(), dim, idx, src, "mean", include_self
                )
                if not torch.allclose(ref, got):
                    raise RuntimeError(
                        "--deterministic: scatter_reduce detour diverges from native "
                        "index_reduce_(mean); the shim is unsafe on this torch build."
                    )

        prior_orig = _ORIG_INDEX_REDUCE
        prior_tensor_method = torch.Tensor.index_reduce_
        env_present = "CUBLAS_WORKSPACE_CONFIG" in os.environ
        prior_env = os.environ.get("CUBLAS_WORKSPACE_CONFIG")
        prior_deterministic = torch.are_deterministic_algorithms_enabled()
        prior_warn_only = (
            torch.is_deterministic_algorithms_warn_only_enabled()
            if hasattr(torch, "is_deterministic_algorithms_warn_only_enabled")
            else False
        )
        prior_cudnn_deterministic = torch.backends.cudnn.deterministic
        prior_cudnn_benchmark = torch.backends.cudnn.benchmark
        prior_cpu_rng = torch.get_rng_state()
        prior_numpy_rng = np.random.get_state()
        cuda_available = torch.cuda.is_available()
        prior_cuda_rng = torch.cuda.get_rng_state_all() if cuda_available else None

        try:
            os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
            torch.use_deterministic_algorithms(True, warn_only=False)
            torch.backends.cudnn.deterministic = True
            torch.backends.cudnn.benchmark = False
            torch.manual_seed(0)
            np.random.seed(0)
            if cuda_available:
                torch.cuda.manual_seed_all(0)
            _ORIG_INDEX_REDUCE = orig
            torch.Tensor.index_reduce_ = _index_reduce_mean_deterministic
            _DONE = True
        except Exception as exc:
            torch.Tensor.index_reduce_ = prior_tensor_method
            _ORIG_INDEX_REDUCE = prior_orig
            _DONE = False
            if env_present:
                assert prior_env is not None
                os.environ["CUBLAS_WORKSPACE_CONFIG"] = prior_env
            else:
                os.environ.pop("CUBLAS_WORKSPACE_CONFIG", None)
            torch.use_deterministic_algorithms(
                prior_deterministic, warn_only=prior_warn_only
            )
            torch.backends.cudnn.deterministic = prior_cudnn_deterministic
            torch.backends.cudnn.benchmark = prior_cudnn_benchmark
            torch.set_rng_state(prior_cpu_rng)
            np.random.set_state(prior_numpy_rng)
            if prior_cuda_rng is not None:
                torch.cuda.set_rng_state_all(prior_cuda_rng)
            raise RuntimeError(
                "--deterministic: strict deterministic setup failed; "
                "process state was restored."
            ) from exc
