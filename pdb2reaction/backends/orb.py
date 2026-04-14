# pdb2reaction/backends/orb.py

"""
ORB (orb-models) backend for pdb2reaction.

Requires: ``pip install "pdb2reaction[orb]"`` (orb-models).
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

from .base import (
    MLIPCalculator,
    BackendError,
    _prepare_model_for_autograd_hessian,
    _restore_model_after_autograd_hessian,
)

# Deprecated model aliases
ORB_DEPRECATED_MODEL_ALIASES = {
    "orb-v1": "orb-v2",
    "orb-d3-v1": "orb-d3-v2",
    "orb-d3-sm-v1": "orb-d3-sm-v2",
    "orb-d3-xs-v1": "orb-d3-xs-v2",
    "orb-v1-mptraj-only": "orb-mptraj-only-v2",
    "orb-mptraj-only-v1": "orb-mptraj-only-v2",
}


def _is_conservative_orb_model(model_name: str) -> bool:
    norm = str(model_name).replace("_", "-").lower()
    return ("conservative" in norm) and ("direct" not in norm)


def _unique_ordered(items):
    seen = set()
    out = []
    for item in items:
        if item not in seen:
            seen.add(item)
            out.append(item)
    return out


class OrbCalculator(MLIPCalculator):
    """ORB backend via orb-models."""

    def __init__(
        self,
        *,
        model: str = "orb_v3_conservative_omol",
        # orb_models expects a precision string such as "float32-high",
        # "float32-highest", or "float64".  The historical default
        # "float32" is not a valid value and silently falls through to
        # the slow "highest" matmul precision, which also enables the
        # donated-buffer aot_autograd optimisation that blocks
        # double-backward Hessians.  Default to "float32-high" so the
        # fast path is active.
        precision: str = "float32-high",
        compile_model: bool = False,
        # Base class parameters
        charge: int = 0,
        spin: int = 1,
        device: str = "auto",
        freeze_atoms: Optional[Sequence[int]] = None,
        hessian_calc_mode: str = "FiniteDifference",
        return_partial_hessian: bool = False,
        hessian_double: bool = True,
        print_timing: bool = True,
        **kwargs,
    ):
        try:
            import torch
            from orb_models.forcefield import pretrained as orb_pretrained
        except Exception as exc:
            raise BackendError(
                "ORB backend requires orb-models and torch. "
                "Install with: pip install \"pdb2reaction[orb]\""
            ) from exc

        super().__init__(
            charge=charge,
            spin=spin,
            device=device,
            freeze_atoms=freeze_atoms,
            hessian_calc_mode=hessian_calc_mode,
            return_partial_hessian=return_partial_hessian,
            hessian_double=hessian_double,
            print_timing=print_timing,
            **kwargs,
        )

        self._torch = torch
        self._pretrained = orb_pretrained

        if str(device).lower() == "auto":
            device = "cuda" if torch.cuda.is_available() else "cpu"
        self.device_str = device

        self.model_name = str(model)
        self.precision = str(precision)
        self.compile_model = bool(compile_model)

        if not _is_conservative_orb_model(self.model_name):
            raise BackendError(
                f"Only conservative Orb models are supported. Requested '{self.model_name}'."
            )

        self._loader = self._resolve_loader(self.model_name)
        self._model_obj, self._adapter = self._load_model()
        self._ase_calc = self._build_ase_calculator()

    def _resolve_loader(self, model_name: str):
        norm_dash = str(model_name).replace("_", "-").lower()
        if norm_dash in ORB_DEPRECATED_MODEL_ALIASES:
            model_name = ORB_DEPRECATED_MODEL_ALIASES[norm_dash]

        if not _is_conservative_orb_model(model_name):
            raise BackendError(f"Only conservative Orb models are supported. Requested '{model_name}'.")

        if hasattr(self._pretrained, "ORB_PRETRAINED_MODELS"):
            model_map = getattr(self._pretrained, "ORB_PRETRAINED_MODELS")
            cands = [model_name, model_name.replace("_", "-"), model_name.replace("-", "_")]
            cands.extend([x.lower() for x in cands])
            for cand in cands:
                if cand in model_map:
                    return model_map[cand]
            lower_map = {str(k).lower(): v for k, v in model_map.items()}
            for cand in cands:
                if str(cand).lower() in lower_map:
                    return lower_map[str(cand).lower()]

        for cand in (model_name, model_name.replace("-", "_"), model_name.lower().replace("-", "_")):
            if hasattr(self._pretrained, cand):
                return getattr(self._pretrained, cand)

        raise BackendError(f"Unknown Orb model '{model_name}'.")

    def _load_model(self):
        bases = [
            {"device": self.device_str, "precision": self.precision, "compile": self.compile_model},
            {"device": self.device_str, "precision": self.precision},
            {"device": self.device_str},
            {},
        ]
        last_exc = None
        for kwargs in bases:
            try:
                out = self._loader(**kwargs)
                if isinstance(out, tuple) and len(out) >= 2:
                    return out[0], out[1]
                return out, None
            except Exception as exc:
                last_exc = exc
                continue
        raise BackendError(f"Failed to load Orb model '{self.model_name}': {last_exc}")

    def _build_ase_calculator(self):
        def _try_construct(cls):
            arg_combos = []
            if self._adapter is not None:
                arg_combos.append(((self._model_obj, self._adapter), {"device": self.device_str}))
                arg_combos.append(((self._model_obj, self._adapter), {}))
            arg_combos.append(((self._model_obj,), {"device": self.device_str}))
            arg_combos.append(((self._model_obj,), {}))
            for args, kw in arg_combos:
                try:
                    return cls(*args, **kw)
                except TypeError:
                    continue
            return None

        try:
            from orb_models.forcefield.inference.calculator import ORBCalculator
            calc = _try_construct(ORBCalculator)
            if calc is not None:
                return calc
        except ImportError:
            pass

        try:
            from orb_models.forcefield.calculator import ORBCalculator
            calc = _try_construct(ORBCalculator)
            if calc is not None:
                return calc
        except ImportError:
            pass

        raise BackendError("Failed to build ORBCalculator.")

    def _compute_energy_forces_ev(self, elem, coord_ang):
        from ase import Atoms

        atoms = Atoms(symbols=list(elem), positions=np.asarray(coord_ang, dtype=np.float64))
        atoms.info["charge"] = int(self.charge)
        atoms.info["spin"] = int(self.mult)
        atoms.calc = self._ase_calc

        energy = float(atoms.get_potential_energy())
        forces = np.asarray(atoms.get_forces(), dtype=np.float64)
        return energy, forces

    # ------------------------------------------------------------------
    # Analytical Hessian — Orb-specific get_hessian() override
    # ------------------------------------------------------------------
    #
    # Why an override instead of providing only
    # ``_compute_analytical_hessian_ev``?  The base-class dispatcher
    # calls ``_compute_energy_forces_ev`` (which routes through the ASE
    # ORBCalculator) *immediately* before ``_compute_analytical_hessian_ev``.
    # That first eval-mode predict compiles the Orb model's backward
    # graph with aot_autograd / donated-buffer optimisations, which
    # then blocks the double-backward that Hessian autograd needs.
    # To sidestep that, we override ``get_hessian`` so that, when the
    # user requests analytical mode, the energy, forces, and Hessian
    # are all computed in a *single* train-mode predict call on a
    # fresh graph.  FD mode delegates to ``super().get_hessian()``
    # unchanged, preserving the validated FD path.

    def _supports_analytical_hessian(self) -> bool:
        return True

    @staticmethod
    def _energy_key(out_dict):
        """Pick the energy field from an Orb predict() output dict."""
        for key in ("energy", "free_energy", "total_energy", "E"):
            if key in out_dict:
                return key
        for key in out_dict:
            if str(key).lower().startswith("energy"):
                return key
        raise BackendError(
            f"Orb predict() output has no energy key; got {sorted(out_dict)!r}."
        )

    # ------------------------------------------------------------------
    # get_hessian override: analytical path entirely in a single
    # train-mode predict call on a fresh graph, avoiding the base
    # dispatcher's prior eval-mode ASE call that taints the model's
    # backward graph with aot_autograd donated buffers.
    # ------------------------------------------------------------------

    def get_hessian(self, elem, coords):
        import time
        from pysisyphus.constants import BOHR2ANG
        from .base import EV2AU, F_EVAA_2_AU

        mode = (self.hessian_calc_mode or "FiniteDifference").strip().lower()
        use_analytical = mode in ("analytical", "analytic")

        # FD mode (and any unsupported mode): delegate to the base
        # dispatcher unchanged.  The FD path has been validated and
        # we do not want to touch it.
        if not (use_analytical and self._supports_analytical_hessian()):
            return super().get_hessian(elem, coords)

        # ------------------- analytical path -------------------
        coord_ang = np.asarray(coords, dtype=np.float64).reshape(-1, 3) * BOHR2ANG
        n_atoms = coord_ang.shape[0]

        t0 = time.perf_counter()
        result = self._analytical_energy_forces_hessian_ev(elem, coord_ang)
        mode_elapsed = time.perf_counter() - t0

        e_ev = result["energy_ev"]
        F_ev = np.asarray(result["forces_ev_ang"], dtype=np.float64).reshape(-1, 3)
        H_ev = np.asarray(result["hessian_ev_ang2"], dtype=np.float64)

        if H_ev.ndim == 4:
            H_ev = H_ev.reshape(n_atoms * 3, n_atoms * 3)

        F_ev = self._zero_frozen_forces_ev(F_ev)
        H_ev = self._apply_active_trim_np(H_ev, n_atoms)
        H_au = self._hessian_ev_to_au(H_ev)

        out = {
            "energy": e_ev * EV2AU,
            "forces": (F_ev.reshape(-1) * F_EVAA_2_AU),
            "hessian": H_au,
        }

        if self.print_timing:
            import click
            click.echo(
                f"[HessianTiming] mode: Analytical | elapsed: {mode_elapsed:.2f} s"
            )

        partial_meta = self._build_partial_hessian_meta(n_atoms)
        if partial_meta is not None:
            out["within_partial_hessian"] = partial_meta

        return out

    def _analytical_energy_forces_hessian_ev(self, elem, coord_ang):
        """Compute E (eV), F (eV/Å), H (eV/Å²) in a single train-mode predict.

        This runs ``model.train(True)`` + autograd Hessian in a fresh
        graph so no prior eval-mode call has compiled donated buffers
        on the model's backward pass.
        """
        import warnings
        from ase import Atoms

        if self.compile_model:
            warnings.warn(
                "OrbCalculator(compile_model=True) may interact poorly with "
                "analytical Hessian autograd. If this call fails, retry "
                "with compile_model=False or --hessian-calc-mode FiniteDifference.",
                stacklevel=2,
            )

        torch = self._torch
        n_atoms = len(elem)

        atoms = Atoms(
            symbols=list(elem),
            positions=np.asarray(coord_ang, dtype=np.float64),
        )
        atoms.info["charge"] = int(self.charge)
        atoms.info["spin"] = int(self.mult)

        def _raise(exc):
            msg = str(exc).lower()
            if "out of memory" in msg and "cuda" in msg:
                raise BackendError(
                    "Orb analytical Hessian failed due to CUDA out-of-memory. "
                    "Retry with --hessian-calc-mode FiniteDifference."
                ) from exc
            raise BackendError(
                f"Orb analytical Hessian failed: {exc}"
            ) from exc

        # Disable donated-buffer optimisation for this Hessian call so
        # the second backward can walk the graph.
        donated_saved = None
        _functorch_config = None
        try:
            from torch._functorch import config as _functorch_config  # type: ignore
            donated_saved = getattr(_functorch_config, "donated_buffer", None)
            _functorch_config.donated_buffer = False
        except Exception:
            _functorch_config = None

        # Reset torch._dynamo compile cache so no previously-compiled
        # eval-mode graph (from an upstream call) blocks our double
        # backward.
        try:
            import torch._dynamo as _dynamo  # type: ignore
            _dynamo.reset()
        except Exception:
            pass

        state = _prepare_model_for_autograd_hessian(self._model_obj, torch)
        try:
            try:
                # ---------- New API path (adapter) ----------
                if self._adapter is not None and hasattr(self._model_obj, "predict"):
                    graph = self._adapter.from_ase_atoms(
                        atoms=atoms, device=self.device_str
                    )
                    if not hasattr(graph, "node_features"):
                        raise BackendError(
                            "Unexpected Orb graph format: missing node_features."
                        )
                    if "positions" not in graph.node_features:
                        raise BackendError(
                            "Unexpected Orb graph format: "
                            "node_features missing 'positions'."
                        )

                    pos = (
                        graph.node_features["positions"]
                        .detach()
                        .clone()
                        .to(self.device_str)
                        .requires_grad_(True)
                    )
                    graph.node_features["positions"] = pos

                    out_dict = self._model_obj.predict(graph)

                # ---------- Legacy API fallback ----------
                else:
                    try:
                        from orb_models.forcefield import atomic_system  # type: ignore
                    except ImportError as exc:
                        raise BackendError(
                            "Orb analytical Hessian: neither the new adapter "
                            "API nor orb_models.forcefield.atomic_system is "
                            f"available.  Import error: {exc}"
                        ) from exc

                    graph = atomic_system.ase_atoms_to_atom_graphs(
                        atoms,
                        getattr(self._model_obj, "system_config", None),
                        device=self.device_str,
                    )
                    pos_attr = None
                    for cand in ("positions", "pos", "coords", "xyz"):
                        if hasattr(graph, cand):
                            pos_attr = cand
                            break
                    if pos_attr is None:
                        raise BackendError(
                            "Orb analytical Hessian: could not locate position "
                            "tensor on the legacy graph object."
                        )

                    pos = (
                        torch.as_tensor(getattr(graph, pos_attr), device=self.device_str)
                        .detach()
                        .clone()
                        .requires_grad_(True)
                    )
                    setattr(graph, pos_attr, pos)

                    out_dict = self._model_obj.predict(graph, split=False)

                E = out_dict[self._energy_key(out_dict)].squeeze()
                if E.ndim > 0:
                    E = E.reshape(-1)[0]

                # Build Hessian via row-by-row autograd.grad
                H = self._hessian_from_energy_grad(E, pos, n_atoms)

                # Forces: grab them from the predict output instead of
                # recomputing — predict already ran
                # compute_gradient_forces_and_stress internally.
                forces_ev_ang = None
                for fkey in ("forces", "grad_forces", "force"):
                    if fkey in out_dict:
                        forces_ev_ang = (
                            out_dict[fkey].detach().cpu().numpy().astype(np.float64)
                        )
                        break
                if forces_ev_ang is None:
                    raise BackendError(
                        "Orb predict() output contained no forces key; "
                        f"got {sorted(out_dict)!r}."
                    )

                return {
                    "energy_ev": float(E.detach().cpu().item()),
                    "forces_ev_ang": forces_ev_ang,
                    "hessian_ev_ang2": H,
                }

            except BackendError:
                raise
            except Exception as exc:
                _raise(exc)
        finally:
            _restore_model_after_autograd_hessian(self._model_obj, state)
            if _functorch_config is not None and donated_saved is not None:
                try:
                    _functorch_config.donated_buffer = donated_saved
                except Exception:
                    pass
            if str(self.device_str).startswith("cuda"):
                try:
                    torch.cuda.empty_cache()
                except Exception:
                    pass

    def _hessian_from_energy_grad(self, E, pos, n_atoms):
        """Build H = d^2 E / dx^2 via two chained torch.autograd.grad calls.

        We deliberately avoid ``torch.autograd.functional.hessian`` because
        PyTorch >= 2.5 routes it through aot_autograd, which does not
        support the double-backward pass required for a Hessian when the
        model's forward sub-ops are wrapped in ``torch.compile`` /
        aot_autograd.  Instead we compute ``grad(E, pos, create_graph=True)``
        once and then walk each element of the gradient through a second
        plain ``torch.autograd.grad`` to assemble the Hessian row by row.
        ``pos`` must be the exact 3-D ``(N,3)`` positions tensor that is
        embedded in the Orb graph (otherwise the backward chain is broken).
        """
        torch = self._torch
        dof = n_atoms * 3

        g = torch.autograd.grad(E, pos, create_graph=True)[0]
        gflat = g.reshape(-1)
        rows = []
        for i in range(dof):
            row = torch.autograd.grad(
                gflat[i],
                pos,
                retain_graph=(i < dof - 1),
                create_graph=False,
                allow_unused=False,
            )[0]
            rows.append(row.reshape(-1))
        H = torch.stack(rows, dim=0)
        H = H.detach().cpu().numpy().astype(np.float64)
        H = 0.5 * (H + H.T)  # enforce symmetry (numerical noise)
        return H


class OrbASECalculator:
    """Factory that returns an ORB ASE calculator for DMF."""

    def __new__(
        cls,
        *,
        model: str = "orb_v3_conservative_omol",
        device: str = "auto",
        precision: str = "float32",
        compile_model: bool = False,
    ):
        # Build and return the ORB ASE calculator directly
        calc = OrbCalculator(
            model=model, device=device, precision=precision,
            compile_model=compile_model, charge=0, spin=1,
        )
        return calc._ase_calc
