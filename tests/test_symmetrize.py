"""Regression coverage for ``pdb2reaction.core.utils.symmetrize_inplace``.

The chunked symmetrizer averages each on-diagonal block via a transposed view.
When an edge length leaves a 1x1 trailing block (``N % chunk == 1``), an
in-place ``add_`` on that strided self-overlapping view is rejected by Torch
("some elements of the input tensor and the written-to tensor refer to a single
memory location"), which used to abort freq / TSopt / IRC Hessian assembly. The
fix writes the diagonal block out-of-place, matching the off-diagonal path. The
canonical reproducer is ``N=513`` (171 active atoms x 3) against the default
``chunk=512``.
"""

import pytest
import torch

from pdb2reaction.core.utils import symmetrize_inplace

# Edge lengths that bracket the chunk=512 boundary: 512/513 straddle it,
# 170/171 reproduce the active-atom count behind the original crash, 1024
# spans two full chunks, and 1 covers the scalar degenerate case.
EDGE_LENGTHS = [1, 170, 171, 512, 513, 600, 1024]
PRECISIONS = [torch.float32, torch.float64]


def _tol(dtype):
    return 1e-6 if dtype == torch.float32 else 1e-12


@pytest.mark.parametrize("dtype", PRECISIONS)
@pytest.mark.parametrize("n", EDGE_LENGTHS)
def test_matches_naive_symmetrization(n, dtype):
    """Output equals ``0.5 * (A + A.T)`` and is itself symmetric."""
    torch.manual_seed(7)
    a = torch.randn(n, n, dtype=dtype)
    expected = (a + a.t()).mul(0.5)
    result = symmetrize_inplace(a.clone())
    assert result.shape == (n, n)
    assert result.dtype == dtype
    assert torch.allclose(result, expected, rtol=0.0, atol=_tol(dtype))
    assert torch.allclose(result, result.t(), rtol=0.0, atol=_tol(dtype))


@pytest.mark.parametrize("dtype", PRECISIONS)
def test_n513_does_not_raise(dtype):
    """The 1x1 trailing-block case (N=513) must symmetrize without raising."""
    torch.manual_seed(7)
    a = torch.randn(513, 513, dtype=dtype)
    result = symmetrize_inplace(a.clone())  # previously raised RuntimeError
    expected = (a + a.t()).mul(0.5)
    assert torch.allclose(result, expected, rtol=0.0, atol=_tol(dtype))
