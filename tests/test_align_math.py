# tests/test_align_math.py
"""Tests for pure math utilities in pdb2reaction.align_freeze_atoms."""

import numpy as np
import pytest

from pysisyphus.Geometry import Geometry
import pdb2reaction.workflows.align_freeze as align_freeze
from pdb2reaction.workflows.align_freeze import (
    align_second_to_first_kabsch_inplace,
    kabsch_R_t,
    _rodrigues,
    _rotation_align_vectors,
    _orth_proj_perp,
)


@pytest.mark.parametrize(
    ("freeze_atoms", "expected"),
    [
        ([], "used 3 atoms"),
        ([0, 1, 2], "used 3 freeze atoms"),
    ],
)
def test_kabsch_log_names_the_alignment_selection(
    monkeypatch, freeze_atoms, expected
):
    coords = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    )
    ref = Geometry(["H", "H", "H"], coords.reshape(-1), coord_type="cart")
    mob = Geometry(
        ["H", "H", "H"],
        (coords + np.array([0.2, -0.1, 0.3])).reshape(-1),
        coord_type="cart",
    )
    ref.freeze_atoms = np.asarray(freeze_atoms, dtype=int)
    mob.freeze_atoms = np.asarray(freeze_atoms, dtype=int)
    messages = []
    monkeypatch.setattr(
        align_freeze,
        "emit",
        lambda message, **kwargs: messages.append(message),
    )

    align_second_to_first_kabsch_inplace(ref, mob, verbose=True)

    assert expected in messages[-1]


class TestKabsch:
    def test_identity_case(self):
        """Same point sets should give identity rotation and zero translation."""
        P = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        R, t = kabsch_R_t(P, P)
        np.testing.assert_allclose(R, np.eye(3), atol=1e-10)
        np.testing.assert_allclose(t, np.zeros(3), atol=1e-10)

    def test_pure_translation(self):
        """Translated point set: R = I, t = shift."""
        P = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        shift = np.array([3.0, -2.0, 1.0])
        Q = P - shift  # Q @ I + shift = P => t = shift
        R, t = kabsch_R_t(P, Q)
        np.testing.assert_allclose(R, np.eye(3), atol=1e-10)
        np.testing.assert_allclose(t, shift, atol=1e-10)

    def test_known_rotation(self):
        """90-degree rotation about z-axis."""
        P = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]])
        # Rotate P by -90 degrees about z to get Q
        theta = -np.pi / 2
        Rz = np.array([[np.cos(theta), -np.sin(theta), 0],
                        [np.sin(theta),  np.cos(theta), 0],
                        [0, 0, 1]])
        Q = P @ Rz.T  # Row-vector rotation
        R, t = kabsch_R_t(P, Q)
        # Q @ R + t ≈ P
        reconstructed = Q @ R + t
        np.testing.assert_allclose(reconstructed, P, atol=1e-10)

    def test_reflection_handling(self):
        """Kabsch should handle reflection case (det(R) < 0 before correction)."""
        P = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0],
                       [1.0, 1.0, 1.0]])
        Q = P.copy()
        Q[:, 2] *= -1  # Reflect z
        R, t = kabsch_R_t(P, Q)
        assert np.linalg.det(R) > 0  # Should be proper rotation

    def test_shape_validation(self):
        with pytest.raises(ValueError, match="shape"):
            kabsch_R_t(np.zeros((3, 2)), np.zeros((3, 2)))


class TestRodrigues:
    def test_zero_angle(self):
        R = _rodrigues(np.array([0.0, 0.0, 1.0]), 0.0)
        np.testing.assert_allclose(R, np.eye(3), atol=1e-10)

    def test_90_about_z(self):
        R = _rodrigues(np.array([0.0, 0.0, 1.0]), np.pi / 2)
        expected = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)
        np.testing.assert_allclose(R, expected, atol=1e-10)

    def test_rotation_is_orthogonal(self):
        R = _rodrigues(np.array([1.0, 1.0, 1.0]), 1.23)
        np.testing.assert_allclose(R @ R.T, np.eye(3), atol=1e-10)
        np.testing.assert_allclose(np.linalg.det(R), 1.0, atol=1e-10)

    def test_zero_axis(self):
        R = _rodrigues(np.zeros(3), 1.0)
        np.testing.assert_allclose(R, np.eye(3), atol=1e-10)


class TestRotationAlignVectors:
    def test_same_direction(self):
        a = np.array([1.0, 0.0, 0.0])
        R = _rotation_align_vectors(a, a)
        np.testing.assert_allclose(R, np.eye(3), atol=1e-10)

    def test_opposite_direction(self):
        a = np.array([1.0, 0.0, 0.0])
        b = np.array([-1.0, 0.0, 0.0])
        R = _rotation_align_vectors(a, b)
        result = R @ a
        np.testing.assert_allclose(result, b, atol=1e-10)

    def test_orthogonal_vectors(self):
        a = np.array([1.0, 0.0, 0.0])
        b = np.array([0.0, 1.0, 0.0])
        R = _rotation_align_vectors(a, b)
        result = R @ a
        np.testing.assert_allclose(result, b, atol=1e-10)

    def test_arbitrary_vectors(self):
        a = np.array([1.0, 2.0, 3.0])
        b = np.array([-1.0, 0.5, 2.0])
        R = _rotation_align_vectors(a, b)
        result = R @ a
        # Should point in same direction as b
        cos_angle = np.dot(result, b) / (np.linalg.norm(result) * np.linalg.norm(b))
        np.testing.assert_allclose(cos_angle, 1.0, atol=1e-10)


class TestOrthProjPerp:
    def test_z_axis(self):
        P = _orth_proj_perp(np.array([0.0, 0.0, 1.0]))
        # Projecting z onto xy-plane should kill z component
        result = P @ np.array([1.0, 2.0, 3.0])
        np.testing.assert_allclose(result, [1.0, 2.0, 0.0], atol=1e-10)

    def test_is_idempotent(self):
        P = _orth_proj_perp(np.array([1.0, 1.0, 0.0]))
        v = np.array([1.0, 2.0, 3.0])
        np.testing.assert_allclose(P @ P @ v, P @ v, atol=1e-10)

    def test_zero_vector(self):
        P = _orth_proj_perp(np.zeros(3))
        np.testing.assert_allclose(P, np.eye(3), atol=1e-10)
