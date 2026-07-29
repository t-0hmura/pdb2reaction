# pdb2reaction/backends/mace.py

"""
MACE (mace-torch) backend for pdb2reaction.

Requires: ``pip install mace-torch`` in a dedicated environment.  Current
``mace-torch`` and ``fairchem-core`` releases require incompatible ``e3nn``
versions, so the MACE and UMA runtimes must not share an environment.
"""

from __future__ import annotations

import hashlib
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
        return_partial_hessian: bool = True,
        hessian_double: bool = True,
        print_timing: bool = True,
        out_hess_torch: bool = True,
        **kwargs,
    ):
        try:
            import torch
        except Exception as exc:
            raise BackendError(
                "MACE backend requires torch and mace-torch. "
                "Install mace-torch in the dedicated MACE environment; "
                "it cannot share an environment with fairchem-core (UMA) "
                "because their e3nn requirements conflict."
            ) from exc

        # Detect the installed distribution without importing fairchem.  Importing
        # it here can itself fail after mace-torch has replaced e3nn, hiding the
        # more useful explanation below.  Current releases conflict regardless of
        # the mace-torch version: mace-torch pins e3nn==0.4.4 while fairchem-core
        # requires e3nn>=0.5.
        try:
            from importlib.metadata import PackageNotFoundError, version

            version("fairchem-core")
        except (PackageNotFoundError, ImportError):
            pass
        else:
            warnings.warn(
                "fairchem-core and mace-torch are installed in the same "
                "environment, but their e3nn requirements conflict. This "
                "environment is unsupported and may fail at import or runtime. "
                "Use separate environments for the UMA and MACE backends; see "
                "the MACE installation guide.",
                RuntimeWarning,
                stacklevel=2,
            )

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
        url_text = str(url)
        basename = os.path.basename(url_text).split("?")[0] or "mace.model"
        digest = hashlib.sha256(url_text.encode("utf-8")).hexdigest()[:16]
        target = os.path.join(
            tempfile.gettempdir(),
            f"pdb2reaction-mace-{digest}-{basename}",
        )
        if not os.path.exists(target):
            fd, staged = tempfile.mkstemp(
                prefix=f".{os.path.basename(target)}.",
                dir=tempfile.gettempdir(),
            )
            os.close(fd)
            try:
                urllib.request.urlretrieve(url_text, staged)
                os.replace(staged, target)
            finally:
                if os.path.exists(staged):
                    os.unlink(staged)
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
            kwargs = {
                "device": self.device_str,
                "default_dtype": self.default_dtype,
            }
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
            return _safe_anicc(
                {
                    "device": self.device_str,
                    "default_dtype": self.default_dtype,
                }
            )

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

    # Analytical Hessian (delegated to MACECalculator.get_hessian)

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
                # detach so model activations / saved tensors free; we only
                # need the Hessian numerically downstream, no gradients.
                H = out["hessian"].detach()
                del out
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
