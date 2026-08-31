"""Stateful PySCF/GPU4PySCF calculator backend."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

import numpy as np

from pysisyphus.constants import ANG2BOHR, AU2EV, BOHR2ANG

from pdb2reaction.backends.base import BackendError, MLIPCalculator
from pdb2reaction.core.dft_settings import DFTSettings, PCM_DIELECTRIC, resolve_dft_settings


def _to_numpy(value):
    if value is None:
        return None
    getter = getattr(value, "get", None)
    if callable(getter):
        value = getter()
    return np.asarray(value)


def _safe_config_value(value: Any) -> bool:
    if value is None or isinstance(value, (str, int, float, bool)):
        return True
    if isinstance(value, (list, tuple)):
        return all(_safe_config_value(item) for item in value)
    if isinstance(value, Mapping):
        return all(isinstance(key, str) and _safe_config_value(item) for key, item in value.items())
    return False


def _apply_attributes(obj: Any, values: Mapping[str, Any], path: str) -> None:
    for name, value in values.items():
        if not _safe_config_value(value):
            raise BackendError(f"{path}.{name} is not a YAML-safe scalar/list/mapping value.")
        if not hasattr(obj, name):
            raise BackendError(f"Unknown PySCF attribute {path}.{name}.")
        current = getattr(obj, name)
        if callable(current):
            raise BackendError(f"PySCF callable attribute {path}.{name} cannot be replaced.")
        try:
            setattr(obj, name, deepcopy_value(value))
        except (TypeError, ValueError) as exc:
            raise BackendError(f"Invalid value for PySCF attribute {path}.{name}: {exc}") from exc


def deepcopy_value(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: deepcopy_value(item) for key, item in value.items()}
    if isinstance(value, list):
        return [deepcopy_value(item) for item in value]
    if isinstance(value, tuple):
        return tuple(deepcopy_value(item) for item in value)
    return value


def _settings_from_mapping(values: Mapping[str, Any]) -> DFTSettings:
    raw = dict(values)
    return resolve_dft_settings(
        {
            "backend": "dft",
            "charge": raw.get("charge", 0),
            "spin": raw.get("multiplicity", 1),
            "dft": raw,
        }
    )


class PySCFDFTSession:
    """One live SCF method, one last-good density, and one exact-coordinate cache."""

    def __init__(self, settings: DFTSettings):
        self.settings = settings
        self._scanner = None
        self._using_rks_lowmem = False
        self._last_key = None
        self._last_symbols: Optional[Tuple[str, ...]] = None
        self._last_coords: Optional[np.ndarray] = None
        self._last_mm_coords: Optional[np.ndarray] = None
        self._last_mm_charges: Optional[np.ndarray] = None
        self._cache: Dict[str, Any] = {}
        self._last_good: Optional[Dict[str, Any]] = None
        self._pending_checkpoint = False
        self._checkpoint_attempted = False
        self.checkpoint_status: Optional[Dict[str, Any]] = None
        self.metrics: list[Dict[str, Any]] = []

    @staticmethod
    def _key(
        symbols: Sequence[str],
        coords_ang: np.ndarray,
        mm_coords_ang: Optional[np.ndarray],
        mm_charges: Optional[np.ndarray],
    ) -> str:
        digest = hashlib.sha256()
        digest.update("\0".join(map(str, symbols)).encode("utf-8"))
        digest.update(np.ascontiguousarray(coords_ang, dtype=np.float64).tobytes())
        if mm_coords_ang is not None:
            digest.update(np.ascontiguousarray(mm_coords_ang, dtype=np.float64).tobytes())
        if mm_charges is not None:
            digest.update(np.ascontiguousarray(mm_charges, dtype=np.float64).tobytes())
        return digest.hexdigest()

    def _make_mol(self, symbols: Sequence[str], coords_ang: np.ndarray):
        try:
            from pyscf import gto, lib
        except ImportError as exc:
            raise BackendError(
                "PySCF is required for --backend dft. Install the project dft extra."
            ) from exc

        lib.num_threads(int(self.settings.nprocs))
        mol = gto.Mole()
        mol.atom = [
            (str(symbol), tuple(map(float, coord)))
            for symbol, coord in zip(symbols, coords_ang)
        ]
        mol.unit = "Angstrom"
        mol.basis = self.settings.basis
        mol.charge = int(self.settings.charge)
        mol.spin = int(self.settings.multiplicity) - 1
        mol.verbose = int(self.settings.verbose)
        if self.settings.memory_mb is not None:
            mol.max_memory = int(self.settings.memory_mb)
        _apply_attributes(mol, self.settings.pyscf.get("mol", {}), "pyscf.mol")
        if (
            not mol.ecp
            and self.settings.basis.casefold().startswith("def2")
            and any(int(gto.charge(str(symbol))) >= 37 for symbol in symbols)
        ):
            mol.ecp = self.settings.basis
        mol.build()
        return mol

    def _build_method(
        self,
        mol,
        mm_coords_ang: Optional[np.ndarray],
        mm_charges: Optional[np.ndarray],
    ) -> None:
        from pyscf import dft, qmmm, scf

        unrestricted = int(self.settings.multiplicity) != 1
        self._using_rks_lowmem = self.settings.use_rks_lowmem
        if self._using_rks_lowmem:
            try:
                from gpu4pyscf.dft import rks_lowmem
            except ImportError as exc:
                raise BackendError(
                    "gpu4pyscf.dft.rks_lowmem is required by the default "
                    "closed-shell GPU low-memory path. Install a compatible "
                    "GPU4PySCF release or explicitly use --no-lowmem."
                ) from exc
            xc = "HF" if self.settings.is_hf else self.settings.functional
            mf = rks_lowmem.RKS(mol, xc=xc)
        elif self.settings.is_hf:
            mf = scf.UHF(mol) if unrestricted else scf.RHF(mol)
        else:
            mf = dft.UKS(mol) if unrestricted else dft.RKS(mol)
            mf.xc = self.settings.functional

        mf.conv_tol = float(self.settings.conv_tol)
        mf.max_cycle = int(self.settings.max_cycle)
        mf.chkfile = None
        _apply_attributes(mf, self.settings.pyscf.get("mf", {}), "pyscf.mf")
        if hasattr(mf, "grids"):
            mf.grids.level = int(self.settings.grid_level)
            _apply_attributes(mf.grids, self.settings.pyscf.get("grids", {}), "pyscf.grids")

        if self.settings.density_fit:
            kwargs = {}
            if self.settings.auxbasis:
                kwargs["auxbasis"] = self.settings.auxbasis
            density_fit_cfg = self.settings.pyscf.get("density_fit", {})
            kwargs.update(
                {key: value for key, value in density_fit_cfg.items() if key != "enabled"}
            )
            mf = mf.density_fit(**kwargs)
            _apply_attributes(mf.with_df, self.settings.pyscf.get("with_df", {}), "pyscf.with_df")

        if self.settings.solvent_model == "pcm":
            mf = mf.PCM()
            solvent_cfg = self.settings.pyscf.get("with_solvent", {})
            if "eps" not in solvent_cfg:
                try:
                    mf.with_solvent.eps = PCM_DIELECTRIC[self.settings.solvent.casefold()]
                except KeyError as exc:
                    raise BackendError(
                        "PCM solvent requires pyscf.with_solvent.eps for an unknown solvent name "
                        f"({self.settings.solvent!r})."
                    ) from exc
            _apply_attributes(mf.with_solvent, solvent_cfg, "pyscf.with_solvent")
        elif self.settings.solvent_model == "smd":
            mf = mf.SMD()
            mf.with_solvent.solvent = self.settings.solvent
            _apply_attributes(
                mf.with_solvent,
                self.settings.pyscf.get("with_solvent", {}),
                "pyscf.with_solvent",
            )

        if self.settings.engine == "gpu" and not self._using_rks_lowmem:
            try:
                mf = mf.to_gpu()
            except (ImportError, ModuleNotFoundError) as exc:
                raise BackendError(
                    "GPU4PySCF is required for --engine gpu. Install a compatible gpu4pyscf package."
                ) from exc
            except Exception as exc:
                raise BackendError(f"Could not move the PySCF method to GPU4PySCF: {exc}") from exc

        if self.settings.embedcharge:
            if mm_coords_ang is None or mm_charges is None:
                raise BackendError("DFT electrostatic embedding requires MM coordinates and charges.")
            if self.settings.engine == "gpu":
                from gpu4pyscf import qmmm as gpu_qmmm

                mf = gpu_qmmm.mm_charge(mf, mm_coords_ang, mm_charges, unit="Angstrom")
            else:
                mf = qmmm.mm_charge(mf, mm_coords_ang, mm_charges, unit="Angstrom")

        mf.chkfile = None
        self._scanner = mf if self._using_rks_lowmem else mf.as_scanner()

    def _update_mm_mol(
        self,
        mm_coords_ang: Optional[np.ndarray],
        mm_charges: Optional[np.ndarray],
    ) -> None:
        if not self.settings.embedcharge or self._scanner is None:
            return
        if mm_coords_ang is None or mm_charges is None:
            raise BackendError("DFT electrostatic embedding requires MM coordinates and charges.")
        if self.settings.engine == "gpu":
            from gpu4pyscf.qmmm import mm_mole
        else:
            from pyscf.qmmm import mm_mole
        self._scanner.mm_mol = mm_mole.create_mm_mol(
            mm_coords_ang, mm_charges, unit="Angstrom"
        )

    def _run_scf(
        self,
        mol,
        mm_coords_ang: Optional[np.ndarray],
        mm_charges: Optional[np.ndarray],
    ) -> Tuple[float, str]:
        had_scanner = self._scanner is not None
        guess_source = (
            "checkpoint"
            if self._pending_checkpoint
            else "previous_density"
            if had_scanner
            else "fresh"
        )
        if self._scanner is None:
            self._build_method(mol, mm_coords_ang, mm_charges)
            if self._last_good is not None:
                self._restore_last_good_into_scanner(mol)
            self._update_mm_mol(mm_coords_ang, mm_charges)
            if self._using_rks_lowmem:
                dm0 = (
                    self._scanner.make_rdm1()
                    if self._last_good is not None
                    else None
                )
                energy = float(self._scanner.kernel(dm0=dm0))
            else:
                energy = float(self._scanner(mol))
        elif self._using_rks_lowmem:
            # rks_lowmem intentionally has no scanner. Rebuild its geometry-
            # bound method while retaining the previous GPU density as dm0.
            dm0 = self._scanner.make_rdm1()
            self._scanner = None
            self._build_method(mol, mm_coords_ang, mm_charges)
            self._update_mm_mol(mm_coords_ang, mm_charges)
            energy = float(self._scanner.kernel(dm0=dm0))
        else:
            self._update_mm_mol(mm_coords_ang, mm_charges)
            energy = float(self._scanner(mol))
        if bool(getattr(self._scanner, "converged", False)):
            self._pending_checkpoint = False
            return energy, guess_source

        # A failed reused guess must not poison the next geometry. Rebuild once
        # and retry from PySCF's native fresh guess.
        self._scanner = None
        saved_last_good = self._last_good
        self._last_good = None
        try:
            self._build_method(mol, mm_coords_ang, mm_charges)
        finally:
            self._last_good = saved_last_good
        self._update_mm_mol(mm_coords_ang, mm_charges)
        energy = float(
            self._scanner.kernel()
            if self._using_rks_lowmem
            else self._scanner(mol)
        )
        if not bool(getattr(self._scanner, "converged", False)):
            self._scanner = None
            advice = (
                " If sufficient GPU and host memory are available, retry with "
                "--no-lowmem; standard density-fitted SCF may converge more robustly."
                if self.settings.lowmem
                else ""
            )
            raise BackendError(
                "PySCF SCF did not converge with either the reused density or a fresh guess."
                + advice
            )
        self._pending_checkpoint = False
        return energy, "fresh_retry"

    def _restore_last_good_into_scanner(self, mol) -> None:
        if self._scanner is None or self._last_good is None:
            return
        def engine_array(value):
            if self.settings.engine == "gpu":
                import cupy

                return cupy.asarray(value)
            return np.asarray(value).copy()

        self._scanner.mo_coeff = engine_array(self._last_good["mo_coeff"])
        self._scanner.mo_occ = engine_array(self._last_good["mo_occ"])
        self._scanner.mo_energy = engine_array(self._last_good["mo_energy"])
        self._scanner._last_mol_fp = mol.ao_loc.copy()

    def _commit_last_good(self) -> None:
        self._last_good = {
            "energy": float(self._scanner.e_tot),
            "mol": self._scanner.mol,
            "mo_coeff": _to_numpy(self._scanner.mo_coeff).copy(),
            "mo_occ": _to_numpy(self._scanner.mo_occ).copy(),
            "mo_energy": _to_numpy(self._scanner.mo_energy).copy(),
        }

    def evaluate(
        self,
        symbols: Sequence[str],
        coords_ang: np.ndarray,
        *,
        need_forces: bool = True,
        mm_coords_ang: Optional[np.ndarray] = None,
        mm_charges: Optional[np.ndarray] = None,
    ) -> Dict[str, Any]:
        coords = np.asarray(coords_ang, dtype=np.float64).reshape(-1, 3)
        mm_coords = (
            None if mm_coords_ang is None else np.asarray(mm_coords_ang, dtype=np.float64).reshape(-1, 3)
        )
        charges = None if mm_charges is None else np.asarray(mm_charges, dtype=np.float64).reshape(-1)
        if not self._checkpoint_attempted:
            self._checkpoint_attempted = True
            if self.settings.checkpoint_path:
                self.load_scf_checkpoint(
                    self.settings.checkpoint_path,
                    (symbols, coords),
                )
        key = self._key(symbols, coords, mm_coords, charges)
        if key == self._last_key and (not need_forces or self._cache.get("forces_au") is not None):
            cached = dict(self._cache)
            cached["cache_hit"] = True
            return cached

        mol = self._make_mol(symbols, coords)
        if key != self._last_key:
            energy, guess_source = self._run_scf(mol, mm_coords, charges)
            self._last_key = key
            self._last_symbols = tuple(map(str, symbols))
            self._last_coords = coords.copy()
            self._last_mm_coords = None if mm_coords is None else mm_coords.copy()
            self._last_mm_charges = None if charges is None else charges.copy()
            self._cache = {
                "energy_au": energy,
                "forces_au": None,
                "mm_forces_au": None,
                "guess_source": guess_source,
                "cache_hit": False,
            }
            if self.settings.save_scf_checkpoint:
                self._commit_last_good()
            else:
                # The live SCF method already owns the density used by its
                # next geometry. Avoid a second quadratic host copy by default.
                self._last_good = None
            self.metrics.append(
                {
                    "guess_source": guess_source,
                    "cycles": int(getattr(self._scanner, "cycles", -1)),
                    "converged": True,
                }
            )

        if need_forces and self._cache.get("forces_au") is None:
            grad = self._scanner.nuc_grad_method()
            if hasattr(grad, "auxbasis_response") and self.settings.density_fit:
                grad.auxbasis_response = True
            gradient = _to_numpy(grad.kernel())
            self._cache["forces_au"] = -gradient.reshape(-1, 3)
            if self.settings.embedcharge:
                dm = self._scanner.make_rdm1()
                mm_gradient = _to_numpy(grad.grad_hcore_mm(dm)) + _to_numpy(grad.grad_nuc_mm())
                self._cache["mm_forces_au"] = -mm_gradient.reshape(-1, 3)

        return dict(self._cache)

    def hessian(self, symbols: Sequence[str], coords_ang: np.ndarray) -> np.ndarray:
        if self.settings.embedcharge:
            raise BackendError(
                "Analytical PySCF Hessian with movable MM point charges is unavailable; "
                "use the complete full-force finite-difference Hessian."
            )
        self.evaluate(symbols, coords_ang, need_forces=True)
        try:
            hess = self._scanner.Hessian()
            if hasattr(hess, "auxbasis_response") and self.settings.density_fit:
                hess.auxbasis_response = 2
            raw = _to_numpy(hess.kernel())
        except Exception as exc:
            raise BackendError(f"PySCF analytical Hessian failed: {exc}") from exc
        n_atoms = len(symbols)
        return raw.transpose(0, 2, 1, 3).reshape(3 * n_atoms, 3 * n_atoms)

    @staticmethod
    def _coerce_atoms(atoms) -> Tuple[Tuple[str, ...], np.ndarray]:
        if hasattr(atoms, "get_chemical_symbols") and hasattr(atoms, "get_positions"):
            return tuple(atoms.get_chemical_symbols()), np.asarray(atoms.get_positions(), dtype=float)
        if hasattr(atoms, "atoms") and hasattr(atoms, "cart_coords"):
            return tuple(map(str, atoms.atoms)), np.asarray(atoms.cart_coords, dtype=float).reshape(-1, 3) * BOHR2ANG
        if isinstance(atoms, tuple) and len(atoms) == 2:
            return tuple(map(str, atoms[0])), np.asarray(atoms[1], dtype=float).reshape(-1, 3)
        raise TypeError("atoms must be ASE Atoms, pysisyphus Geometry, or (symbols, coordinates).")

    def _metadata(self, symbols: Sequence[str], coords_ang: np.ndarray) -> Dict[str, Any]:
        identity = json.dumps(self.settings.scientific_identity(), sort_keys=True, separators=(",", ":"))
        return {
            "schema": 1,
            "symbols": list(map(str, symbols)),
            "coordinates_angstrom": np.asarray(coords_ang, dtype=float).tolist(),
            "scientific_identity_sha256": hashlib.sha256(identity.encode("utf-8")).hexdigest(),
        }

    def save_scf_checkpoint(self, path, atoms) -> Path:
        if not self.settings.save_scf_checkpoint:
            raise BackendError("SCF checkpoint saving is disabled; enable --save-scf-checkpoint.")
        symbols, coords = self._coerce_atoms(atoms)
        if self._key(symbols, coords, None, None) != self._last_key:
            self.evaluate(symbols, coords, need_forces=False)
        from pyscf.scf import chkfile

        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        tmp = destination.with_name(destination.name + ".tmp")
        metadata_path = destination.with_suffix(destination.suffix + ".json")
        metadata_tmp = metadata_path.with_name(metadata_path.name + ".tmp")
        record = self._last_good
        if record is None:
            self._commit_last_good()
            record = self._last_good
        chkfile.dump_scf(
            record["mol"],
            str(tmp),
            float(record["energy"]),
            record["mo_energy"],
            record["mo_coeff"],
            record["mo_occ"],
        )
        metadata_tmp.write_text(
            json.dumps(self._metadata(symbols, coords), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(tmp, destination)
        os.replace(metadata_tmp, metadata_path)
        return destination

    def save_last_checkpoint(self, path=None) -> Optional[Path]:
        """Persist the last converged electronic state without re-evaluation."""

        destination = path or self.settings.checkpoint_path
        if destination is None or self._last_symbols is None or self._last_coords is None:
            return None
        return self.save_scf_checkpoint(
            destination,
            (self._last_symbols, self._last_coords),
        )

    def load_scf_checkpoint(self, path, atoms) -> bool:
        destination = Path(path)
        metadata_path = destination.with_suffix(destination.suffix + ".json")
        if not destination.is_file() or not metadata_path.is_file():
            self.checkpoint_status = {
                "loaded": False,
                "reason": "checkpoint_or_metadata_missing",
                "path": str(destination),
            }
            return False
        symbols, coords = self._coerce_atoms(atoms)
        try:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            expected = self._metadata(symbols, coords)
            if metadata.get("symbols") != expected["symbols"]:
                self.checkpoint_status = {
                    "loaded": False, "reason": "symbols_mismatch", "path": str(destination)
                }
                return False
            if metadata.get("scientific_identity_sha256") != expected["scientific_identity_sha256"]:
                self.checkpoint_status = {
                    "loaded": False, "reason": "settings_mismatch", "path": str(destination)
                }
                return False
            saved_coords = np.asarray(metadata.get("coordinates_angstrom"), dtype=float)
            if saved_coords.shape != coords.shape or not np.allclose(saved_coords, coords, atol=1.0e-8, rtol=0.0):
                self.checkpoint_status = {
                    "loaded": False, "reason": "coordinates_mismatch", "path": str(destination)
                }
                return False
            from pyscf.scf import chkfile

            loaded_mol, record = chkfile.load_scf(str(destination))
            self._last_good = {
                "energy": float(record["e_tot"]),
                "mol": loaded_mol,
                "mo_energy": _to_numpy(record["mo_energy"]).copy(),
                "mo_coeff": _to_numpy(record["mo_coeff"]).copy(),
                "mo_occ": _to_numpy(record["mo_occ"]).copy(),
            }
            self._scanner = None
            self._pending_checkpoint = True
            self._checkpoint_attempted = True
            self._last_key = None
            self._cache = {}
            self.checkpoint_status = {
                "loaded": True, "reason": "matched", "path": str(destination)
            }
            return True
        except (OSError, ValueError, KeyError, TypeError):
            self.checkpoint_status = {
                "loaded": False, "reason": "invalid_checkpoint", "path": str(destination)
            }
            return False


class DFTCalculator(MLIPCalculator):
    """PySisyphus calculator adapter using a stateful PySCF session."""

    def __init__(self, *, dft_settings: Mapping[str, Any], **kwargs):
        settings = _settings_from_mapping(dft_settings)
        kwargs.pop("charge", None)
        kwargs.pop("spin", None)
        super().__init__(
            charge=settings.charge,
            spin=settings.multiplicity,
            **kwargs,
        )
        self.settings = settings
        self.session = PySCFDFTSession(settings)
        self.device_str = settings.engine
        self._closed = False

    def _compute_energy_forces_ev(self, elem, coord_ang):
        result = self.session.evaluate(elem, coord_ang, need_forces=True)
        energy_ev = float(result["energy_au"]) * AU2EV
        forces_ev_ang = np.asarray(result["forces_au"]) * AU2EV * ANG2BOHR
        return energy_ev, forces_ev_ang

    def _supports_analytical_hessian(self) -> bool:
        return not self.settings.embedcharge

    def _compute_analytical_hessian_ev(self, elem, coord_ang):
        hessian_au = self.session.hessian(elem, coord_ang)
        return hessian_au * AU2EV * ANG2BOHR * ANG2BOHR

    def _build_fd_hessian_cpu(self, elem, coord_ang, *, eps_ang=1.0e-3):
        result = super()._build_fd_hessian_cpu(
            elem, coord_ang, eps_ang=eps_ang
        )
        self.session.evaluate(elem, coord_ang, need_forces=True)
        return result

    def save_scf_checkpoint(self, path, atoms) -> Path:
        return self.session.save_scf_checkpoint(path, atoms)

    def load_scf_checkpoint(self, path, atoms) -> bool:
        return self.session.load_scf_checkpoint(path, atoms)

    def close(self) -> None:
        if self._closed:
            return
        if self.settings.save_scf_checkpoint:
            self.session.save_last_checkpoint()
        self._closed = True

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass


try:
    from ase.calculators.calculator import Calculator as ASECalculator, all_changes

    class DFTASECalculator(ASECalculator):
        implemented_properties = ["energy", "forces"]

        def __init__(self, *, dft_settings: Mapping[str, Any], **kwargs):
            super().__init__(**kwargs)
            self.settings = _settings_from_mapping(dft_settings)
            self.session = PySCFDFTSession(self.settings)
            self._closed = False

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            result = self.session.evaluate(
                atoms.get_chemical_symbols(),
                atoms.get_positions(),
                need_forces="forces" in properties,
            )
            self.results = {"energy": float(result["energy_au"]) * AU2EV}
            if "forces" in properties:
                self.results["forces"] = np.asarray(result["forces_au"]) * AU2EV * ANG2BOHR

        def save_scf_checkpoint(self, path=None, atoms=None):
            if atoms is not None:
                self.session.evaluate(
                    atoms.get_chemical_symbols(),
                    atoms.get_positions(),
                    need_forces=False,
                )
            return self.session.save_last_checkpoint(path)

        def load_scf_checkpoint(self, path, atoms) -> bool:
            return self.session.load_scf_checkpoint(path, atoms)

        def close(self) -> None:
            if self._closed:
                return
            if self.settings.save_scf_checkpoint:
                self.session.save_last_checkpoint()
            self._closed = True

        def __del__(self):
            try:
                self.close()
            except Exception:
                pass

except ImportError:  # pragma: no cover - ASE is a project dependency
    DFTASECalculator = None


__all__ = ["DFTASECalculator", "DFTCalculator", "PySCFDFTSession"]
