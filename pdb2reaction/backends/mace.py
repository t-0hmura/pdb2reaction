# pdb2reaction/backends/mace.py

"""
MACE (mace-torch) backend for pdb2reaction.

Requires: ``pip install mace-torch`` (v0.3.8+ coexists with fairchem-core;
older versions need a separate env due to e3nn pinning).
"""

from __future__ import annotations

import os
import tempfile
import urllib.request
import warnings
from typing import Optional, Sequence

import numpy as np

from .base import MLIPCalculator, BackendError

MACE_MP_ALIASES_FALLBACK = [
    "small", "medium", "large",
    "small-0b", "medium-0b",
    "small-0b2", "medium-0b2", "large-0b2",
    "medium-0b3", "medium-mpa-0",
    "small-omat-0", "medium-omat-0",
    "mace-matpes-pbe-0", "mace-matpes-r2scan-0",
    "mh-0", "mh-1",
]


class MACECalculator(MLIPCalculator):
    """MACE backend via mace.calculators."""

    def __init__(
        self,
        *,
        model: str = "MACE-OMOL-0",
        default_dtype: str = "float64",
        # Base class parameters
        charge: int = 0,
        spin: int = 1,
        device: str = "auto",
        freeze_atoms: Optional[Sequence[int]] = None,
        hessian_calc_mode: str = "FiniteDifference",
        return_partial_hessian: bool = False,
        hessian_double: bool = True,
        print_timing: bool = True,
        out_hess_torch: bool = False,
        **kwargs,
    ):
        try:
            import torch
        except Exception as exc:
            raise BackendError(
                "MACE backend requires torch and mace-torch. "
                "Install with: pip install mace-torch "
                "(v0.3.8+ coexists with fairchem-core; older versions need separate env)"
            ) from exc

        # Warn about potential UMA/MACE conflict (resolved in MACE v0.3.8+)
        try:
            import fairchem.core  # noqa: F401
            try:
                from importlib.metadata import version as _pkg_version
                _mace_ver = _pkg_version("mace-torch")
                from packaging.version import Version
                _old_mace = Version(_mace_ver) < Version("0.3.8")
            except Exception:
                _old_mace = True  # assume old if we cannot determine version
            if _old_mace:
                warnings.warn(
                    "fairchem-core is installed alongside mace-torch < v0.3.8. "
                    "These packages conflict due to e3nn version pinning. "
                    "Upgrade to mace-torch >= 0.3.8 (PR #589) or use "
                    "separate environments for UMA and MACE backends.",
                    stacklevel=2,
                )
        except ImportError:
            pass

        super().__init__(
            charge=charge,
            spin=spin,
            device=device,
            freeze_atoms=freeze_atoms,
            hessian_calc_mode=hessian_calc_mode,
            return_partial_hessian=return_partial_hessian,
            hessian_double=hessian_double,
            print_timing=print_timing,
            out_hess_torch=out_hess_torch,
            **kwargs,
        )

        self._torch = torch
        if str(device).lower() == "auto":
            device = "cuda" if torch.cuda.is_available() else "cpu"
        self.device_str = str(device)

        self.model_spec = str(model)
        self.default_dtype = str(default_dtype)

        self._calc = self._build_calc(self.model_spec)

    @staticmethod
    def _download_to_tmp(url: str) -> str:
        target = os.path.join(
            tempfile.gettempdir(),
            os.path.basename(str(url)).split("?")[0] or "mace.model",
        )
        if not os.path.exists(target):
            urllib.request.urlretrieve(str(url), target)
        return target

    def _build_calc(self, model_spec: str):
        try:
            from mace.calculators import mace_anicc, mace_mp, mace_off, mace_omol
        except Exception as exc:
            raise BackendError(
                "Could not import mace.calculators. "
                "Install with: pip install mace-torch (separate env required)"
            ) from exc

        mp_aliases = []
        try:
            from mace.calculators.foundations_models import mace_mp_urls
            mp_aliases = list(mace_mp_urls.keys())
        except Exception:
            mp_aliases = list(MACE_MP_ALIASES_FALLBACK)

        spec = str(model_spec).strip()
        spec_l = spec.lower()
        mp_alias_lookup = {str(x).lower(): x for x in mp_aliases}

        def _mk_from_path(path_or_url):
            from mace.calculators.mace import MACECalculator as _MC
            model_path = str(path_or_url)
            if model_path.startswith("http://") or model_path.startswith("https://"):
                model_path = self._download_to_tmp(model_path)
            return _MC(model_paths=model_path, device=self.device_str, default_dtype=self.default_dtype)

        def _safe_anicc(kwargs):
            try:
                return mace_anicc(**kwargs)
            except Exception as exc:
                msg = str(exc).lower()
                if ("no nvidia driver" in msg) or ("cuda" in msg and "backend" in msg):
                    raise BackendError(
                        "MACE ANICC could not be loaded. Try a CUDA-enabled node or another MACE model."
                    ) from exc
                raise

        # Prefix forms
        if spec_l.startswith("mp:"):
            alias = spec.split(":", 1)[1].strip() or None
            if alias is not None:
                alias = mp_alias_lookup.get(str(alias).lower(), alias)
            return mace_mp(model=alias, device=self.device_str, default_dtype=self.default_dtype)

        if spec_l.startswith("off:"):
            alias = spec.split(":", 1)[1].strip() or None
            if alias is not None:
                alias = str(alias).lower()
            return mace_off(model=alias, device=self.device_str, default_dtype=self.default_dtype)

        if spec_l.startswith("omol:"):
            alias = spec.split(":", 1)[1].strip() or None
            if alias == "":
                alias = None
            if alias is not None:
                alias_l = str(alias).lower()
                if alias_l in ("mace-omol-0", "mace_omol_0", "maceomol0"):
                    alias = "extra_large"
                else:
                    alias = alias_l
            return mace_omol(model=alias, device=self.device_str, default_dtype=self.default_dtype)

        if spec_l.startswith("anicc"):
            path = None
            if ":" in spec:
                path = spec.split(":", 1)[1].strip() or None
            kwargs = {"device": self.device_str}
            if path:
                kwargs["model_path"] = path
            return _safe_anicc(kwargs)

        # Alias forms
        if spec_l in mp_alias_lookup:
            return mace_mp(model=mp_alias_lookup[spec_l], device=self.device_str, default_dtype=self.default_dtype)

        if spec_l in ("off-small", "off-medium", "off-large"):
            alias = spec_l.split("-", 1)[1]
            return mace_off(model=alias, device=self.device_str, default_dtype=self.default_dtype)

        if spec_l in ("omol-extra_large", "extra_large", "mace-omol-0", "mace_omol_0", "maceomol0"):
            return mace_omol(model="extra_large", device=self.device_str, default_dtype=self.default_dtype)

        if spec_l in ("anicc", "ani", "ani500k"):
            return _safe_anicc({"device": self.device_str})

        # Local file / URL
        if os.path.exists(spec) or spec.startswith("http://") or spec.startswith("https://"):
            return _mk_from_path(spec)

        raise BackendError(f"Unknown MACE model spec '{spec}'.")

    def _compute_energy_forces_ev(self, elem, coord_ang):
        from ase import Atoms

        atoms = Atoms(symbols=list(elem), positions=np.asarray(coord_ang, dtype=np.float64))
        atoms.info["charge"] = int(self.charge)
        atoms.info["spin"] = int(self.mult)
        atoms.calc = self._calc

        energy = float(atoms.get_potential_energy())
        forces = np.asarray(atoms.get_forces(), dtype=np.float64)
        return energy, forces

    # ------------------------------------------------------------------
    # Analytical Hessian (delegated to MACECalculator.get_hessian)
    # ------------------------------------------------------------------

    def _supports_analytical_hessian(self) -> bool:
        return True

    def _compute_analytical_hessian_ev(self, elem, coord_ang):
        """Return MACE's analytical Hessian as a torch.Tensor on GPU.

        The public ``MACECalculator.get_hessian(atoms=...)`` API calls
        ``hessian.detach().cpu().numpy()`` internally, forcing a CPU
        round-trip.  To preserve the GPU tensor (so that the base
        dispatcher's torch path can use it and ``out_hess_torch=True``
        works), we reproduce the single-model code path directly:

            batch = calc._atoms_to_batch(atoms)
            out   = calc.models[0](
                calc._clone_batch(batch).to_dict(),
                compute_hessian=True, compute_stress=False,
                training=calc.use_compile,
            )
            return out["hessian"]  # torch.Tensor on calc.device

        This uses the same private helpers that ``get_hessian`` uses and
        is therefore tightly coupled to mace-torch >=0.3.8.  A fallback
        to the public numpy API is kept for older or non-standard MACE
        calculators.
        """
        from ase import Atoms

        if not hasattr(self._calc, "get_hessian"):
            raise BackendError(
                "Installed MACE calculator does not expose get_hessian(). "
                "Upgrade mace-torch (>=0.3.8) or use "
                "--hessian-calc-mode FiniteDifference."
            )

        atoms = Atoms(symbols=list(elem),
                      positions=np.asarray(coord_ang, dtype=np.float64))
        atoms.info["charge"] = int(self.charge)
        atoms.info["spin"] = int(self.mult)
        atoms.calc = self._calc

        calc = self._calc
        torch = self._torch

        try:
            have_internal_path = (
                hasattr(calc, "_atoms_to_batch")
                and hasattr(calc, "_clone_batch")
                and getattr(calc, "models", None) is not None
                and len(calc.models) == 1
            )
            if have_internal_path:
                batch = calc._atoms_to_batch(atoms)
                out = calc.models[0](
                    calc._clone_batch(batch).to_dict(),
                    compute_hessian=True,
                    compute_stress=False,
                    training=getattr(calc, "use_compile", False),
                )
                H = out["hessian"]
            else:
                # Fallback: public API (CPU round-trip).
                H = calc.get_hessian(atoms=atoms)
        except BackendError:
            raise
        except Exception as exc:
            msg = str(exc).lower()
            if "out of memory" in msg and "cuda" in msg:
                raise BackendError(
                    "MACE analytical Hessian failed due to CUDA out-of-memory. "
                    "Retry with --hessian-calc-mode FiniteDifference."
                ) from exc
            raise BackendError(
                f"MACE analytical Hessian failed: {exc}"
            ) from exc

        dof = len(elem) * 3
        if isinstance(H, torch.Tensor):
            # Common shapes: (1, N, 3, N, 3), (N, 3, N, 3), (3N, 3N)
            if H.dim() == 5 and H.size(0) > 0:
                H = H[0]
            return H.reshape(dof, dof)

        H_np = np.asarray(H, dtype=np.float64)
        if H_np.ndim == 5 and H_np.shape[0] > 0:
            H_np = H_np[0]
        return H_np.reshape(dof, dof)


class MACEASECalculator:
    """Factory that returns a MACE ASE calculator for DMF."""

    def __new__(
        cls,
        *,
        model: str = "MACE-OMOL-0",
        device: str = "auto",
        default_dtype: str = "float64",
    ):
        calc = MACECalculator(
            model=model, device=device, default_dtype=default_dtype,
            charge=0, spin=1,
        )
        return calc._calc
