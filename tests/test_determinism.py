from __future__ import annotations

import os

import pytest
import torch


def test_failed_shim_self_check_does_not_commit_state(monkeypatch) -> None:
    from pdb2reaction.backends import _determinism

    original = torch.Tensor.index_reduce_
    monkeypatch.setattr(_determinism, "_DONE", False)
    monkeypatch.setattr(_determinism, "_ORIG_INDEX_REDUCE", None)
    monkeypatch.setattr(torch, "allclose", lambda *_args, **_kwargs: False)

    with pytest.raises(RuntimeError, match="shim is unsafe"):
        _determinism.setup_deterministic()

    assert _determinism._ORIG_INDEX_REDUCE is None
    assert _determinism.is_deterministic_active() is False
    assert torch.Tensor.index_reduce_ is original


def test_failed_strict_activation_rolls_back_and_can_retry(monkeypatch) -> None:
    from pdb2reaction.backends import _determinism

    prior_done = _determinism._DONE
    prior_orig = _determinism._ORIG_INDEX_REDUCE
    prior_method = torch.Tensor.index_reduce_
    env_present = "CUBLAS_WORKSPACE_CONFIG" in os.environ
    prior_env = os.environ.get("CUBLAS_WORKSPACE_CONFIG")
    prior_deterministic = torch.are_deterministic_algorithms_enabled()
    prior_warn_only = torch.is_deterministic_algorithms_warn_only_enabled()
    prior_cudnn_deterministic = torch.backends.cudnn.deterministic
    prior_cudnn_benchmark = torch.backends.cudnn.benchmark
    prior_rng = torch.get_rng_state().clone()
    prior_cuda_rng = (
        torch.cuda.get_rng_state_all() if torch.cuda.is_available() else None
    )
    original_manual_seed = torch.manual_seed
    failed = False

    def fail_once(seed):
        nonlocal failed
        if not failed:
            failed = True
            raise RuntimeError("injected seed failure")
        return original_manual_seed(seed)

    monkeypatch.setattr(torch, "manual_seed", fail_once)
    _determinism._DONE = False
    _determinism._ORIG_INDEX_REDUCE = None
    try:
        with pytest.raises(RuntimeError, match="process state was restored"):
            _determinism.setup_deterministic()

        assert _determinism._DONE is False
        assert _determinism._ORIG_INDEX_REDUCE is None
        assert torch.Tensor.index_reduce_ is prior_method
        assert torch.are_deterministic_algorithms_enabled() is prior_deterministic
        assert torch.is_deterministic_algorithms_warn_only_enabled() is prior_warn_only
        assert torch.backends.cudnn.deterministic is prior_cudnn_deterministic
        assert torch.backends.cudnn.benchmark is prior_cudnn_benchmark
        assert torch.equal(torch.get_rng_state(), prior_rng)
        assert ("CUBLAS_WORKSPACE_CONFIG" in os.environ) is env_present
        assert os.environ.get("CUBLAS_WORKSPACE_CONFIG") == prior_env

        _determinism.setup_deterministic()
        assert _determinism.is_deterministic_active() is True
    finally:
        torch.Tensor.index_reduce_ = prior_method
        _determinism._DONE = prior_done
        _determinism._ORIG_INDEX_REDUCE = prior_orig
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
        torch.set_rng_state(prior_rng)
        if prior_cuda_rng is not None:
            torch.cuda.set_rng_state_all(prior_cuda_rng)
