"""Central configuration defaults for pdb2reaction workflows."""

from __future__ import annotations

from copy import deepcopy
from typing import Any, Dict, Mapping, Optional
from pysisyphus.tr_projection import DEFAULT_TR_PROJECTION

# Convergence preset choices (single source of truth so opt/tsopt/
# path_opt/path_search/scan/all all reject unknown --thresh values at
# Click parse time).

THRESH_CHOICES = ("gau_loose", "gau", "gau_tight", "gau_vtight", "baker", "never")


OUT_DIR_OPT = "./result_opt/"
OUT_DIR_SCAN = "./result_scan/"
OUT_DIR_SCAN2D = "./result_scan2d/"
OUT_DIR_SCAN3D = "./result_scan3d/"
OUT_DIR_FREQ = "./result_freq/"
OUT_DIR_IRC = "./result_irc/"
OUT_DIR_TSOPT = "./result_tsopt/"
OUT_DIR_PATH_OPT = "./result_path_opt/"
OUT_DIR_PATH_SEARCH = "./result_path_search/"
OUT_DIR_ALL = "./result_all/"
OUT_DIR_DFT = "./result_dft/"
OUT_DIR_SP = "./result_sp/"

# `all`-pipeline layout dirnames — single source for the deliverable/scratch split.
SEGMENTS_DIRNAME = "segments"  # per-segment deliverables: segments/seg_NN/<subcmd>/
WORK_DIRNAME = "_work"  # pipeline-wide scratch (rm -rf safe)


GEOM_KW_DEFAULT: Dict[str, Any] = {
    "coord_type": "cart",
    "freeze_atoms": [],
    "tr_projection": DEFAULT_TR_PROJECTION,
}

# Calculator defaults (UMA)

# Single source of truth for the default UMA model name. Every code default
# (backend factory, calculators, workflows) references this instead of a
# hardcoded literal so the release model can be bumped in one place.
DEFAULT_UMA_MODEL = "uma-s-1p2"

CALC_KW_DEFAULT: Dict[str, Any] = {
    "backend": "uma",
    "charge": 0,
    "spin": 1,
    "model": DEFAULT_UMA_MODEL,
    "task_name": "omol",
    "device": "auto",
    "max_neigh": None,
    "radius": None,
    "r_edges": False,
    "workers": 1,
    "workers_per_node": 1,
    # "auto" (resolve per backend: uma/aimnet2 fp32, orb/mace fp64 — see
    # backends._BACKEND_DEFAULT_PRECISION) | "fp32" | "fp64". The sentinel keeps
    # "user named no precision" distinguishable from an explicit "--precision fp32";
    # a literal "fp32" here dispatched reduced precision to ORB and MACE.
    "precision": "auto",
    "hessian_calc_mode": "FiniteDifference",
    "out_hess_torch": True,
    "hessian_double": True,
    "print_timing": True,
    "print_vram": True,
    "return_partial_hessian": True,
    # Solvent correction
    "solvent": "none",
    "solvent_model": "alpb",
    "xtb_cmd": "xtb",
    "xtb_acc": 0.2,
}

# Extended UMA calculator defaults with freeze_atoms support
UMA_CALC_KW: Dict[str, Any] = {
    **CALC_KW_DEFAULT,
    "freeze_atoms": None,
}


ORB_BACKEND_DEFAULTS: Dict[str, Any] = {
    "model": "orb_v3_conservative_omol",
    # float64, not orb_models' "float32-high": that mode is TF32 matmul, whose
    # force noise inflates finite-difference Hessians into spurious imaginary
    # modes. Mirrors backends._BACKEND_DEFAULT_PRECISION["orb"] = "fp64".
    "precision": "float64",
    "compile_model": False,
}

MACE_BACKEND_DEFAULTS: Dict[str, Any] = {
    "model": "MACE-OMOL-0",
    "default_dtype": "float64",
}

AIMNET2_BACKEND_DEFAULTS: Dict[str, Any] = {
    "model": "aimnet2",
}

_BACKEND_DEFAULTS_MAP: Dict[str, Dict[str, Any]] = {
    "orb": ORB_BACKEND_DEFAULTS,
    "mace": MACE_BACKEND_DEFAULTS,
    "aimnet2": AIMNET2_BACKEND_DEFAULTS,
}


def apply_backend_defaults(cfg: Dict[str, Any]) -> None:
    """Apply backend-specific defaults to *cfg* in-place.

    Only overwrite keys whose current value equals the UMA default
    (i.e., the caller has not explicitly set them via CLI or YAML).
    """
    backend = str(cfg.get("backend", "uma"))
    if backend == "auto":
        from pdb2reaction.backends import resolve_backend

        backend = resolve_backend(backend)
        cfg["backend"] = backend
    defaults = _BACKEND_DEFAULTS_MAP.get(backend)
    if defaults is None:
        return
    for key, val in defaults.items():
        if key not in cfg or cfg[key] == CALC_KW_DEFAULT.get(key):
            cfg[key] = val


# Optimizer base (common to LBFGS & RFO)

OPT_BASE_KW: Dict[str, Any] = {
    "thresh": "gau",
    "max_cycles": 10000,
    "print_every": 100,
    "min_step_norm": 1e-8,
    "assert_min_step": True,
    "rms_force": None,
    "rms_force_only": False,
    "max_force_only": False,
    "force_only": False,
    "converge_to_geom_rms_thresh": 0.05,
    "overachieve_factor": 0.0,
    "check_eigval_structure": False,
    "energy_plateau": True,
    "energy_plateau_thresh": 1e-4,
    "energy_plateau_window": 50,
    "line_search": True,
    "dump": False,
    "dump_restart": False,
    "prefix": "",
    "out_dir": OUT_DIR_OPT,
}


LBFGS_KW: Dict[str, Any] = {
    **OPT_BASE_KW,
    "keep_last": 7,
    "beta": 1.0,
    "gamma_mult": False,
    "max_step": 0.30,
    "control_step": True,
    "double_damp": True,
    "mu_reg": None,
    "max_mu_reg_adaptions": 10,
    "reject_uphill": False,
    "uphill_tolerance": 1e-4,
    "rejection_step_floor": 1e-7,
    "max_rejections_at_floor": 3,
}


RFO_KW: Dict[str, Any] = {
    **OPT_BASE_KW,
    "trust_radius": 0.10,
    "trust_update": True,
    "trust_min": 1e-4,
    "trust_max": 0.10,
    "max_energy_incr": None,
    "reject_uphill": False,
    "uphill_tolerance": 1e-4,
    "rejection_trust_floor": 1e-7,
    "max_rejections_at_floor": 3,
    "hessian_update": "bfgs",
    "hessian_init": "calc",
    "hessian_recalc": 500,
    "hessian_recalc_adapt": None,
    "small_eigval_thresh": 1e-8,
    "alpha0": 1.0,
    "max_micro_cycles": 50,
    "rfo_overlaps": False,
    "gediis": False,
    "gdiis": True,
    "gdiis_thresh": 2.5e-3,
    "gediis_thresh": 1.0e-2,
    "gdiis_test_direction": True,
    "adapt_step_func": True,
}

# Bias (harmonic well) defaults

BIAS_KW: Dict[str, Any] = {
    "k": 300,
}


BOND_KW: Dict[str, Any] = {
    "device": "auto",
    "bond_factor": 1.20,
    "margin_fraction": 0.05,
    "delta_fraction": 0.05,
}


OPT_MODE_ALIASES = (
    (("grad", "light", "lbfgs"), "lbfgs"),
    (("hess", "heavy", "rfo"), "rfo"),
)

# DMF (Direct Max Flux + (C)FB-ENM) defaults

DMF_KW: Dict[str, Any] = {
    "backend": "gpu",  # DMF backend: "gpu" (dmf.torch / CUDA, default) or "cpu" (dmf / NumPy)
    "max_cycles": 300,
    "correlated": True,
    "sequential": True,
    "fbenm_only_endpoints": False,
    "fbenm_options": {
        "delta_scale": 0.2,
        "bond_scale": 1.25,
        "fix_planes": True,
    },
    "cfbenm_options": {
        "bond_scale": 1.25,
        "corr0_scale": 1.10,
        "corr1_scale": 1.50,
        "corr2_scale": 1.60,
        "eps": 0.05,
        "pivotal": True,
        "single": True,
        "remove_fourmembered": True,
    },
    "dmf_options": {
        "remove_rotation_and_translation": False,
        "mass_weighted": False,
        "parallel": False,
        "eps_vel": 0.01,
        "eps_rot": 0.01,
        "beta": 10.0,
        "update_teval": False,
    },
    "k_fix": 300.0,
}


def fresh_dmf_config(
    overrides: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    """Return an invocation-local DMF configuration.

    ``DMF_KW`` contains nested mappings.  A shallow ``dict(DMF_KW)`` copy
    therefore shares those nested mappings with the canonical defaults and
    can leak one request's YAML values into later in-process requests.  This
    factory owns the recursive copy and merge so every request starts from an
    independent tree.
    """

    config = deepcopy(DMF_KW)

    def _merge(dst: Dict[str, Any], src: Mapping[str, Any]) -> None:
        for key, value in src.items():
            if isinstance(value, Mapping) and isinstance(dst.get(key), dict):
                _merge(dst[key], value)
            else:
                dst[key] = deepcopy(value)

    if overrides:
        _merge(config, overrides)
    return config

# GrowingString (path representation) defaults

GS_KW: Dict[str, Any] = {
    "fix_first": True,
    "fix_last": True,
    "max_nodes": 20,
    "perp_thresh": 5e-3,
    "reparam_check": "rms",
    "reparam_every": 1,
    "reparam_every_full": 1,
    "param": "equi",
    "max_micro_cycles": 10,
    "reset_dlc": True,
    "climb": True,
    "climb_rms": 5e-4,
    "climb_lanczos": True,
    "climb_lanczos_rms": 5e-4,
    "climb_fixed": False,
    "scheduler": None,
}

# StringOptimizer (optimization control) defaults

STOPT_KW: Dict[str, Any] = {
    "type": "string",
    "thresh": "gau_loose",
    "stop_in_when_full": 300,
    "align": False,
    "scale_step": "global",
    "max_cycles": 300,
    "dump": False,
    "dump_restart": False,
    "reparam_thresh": 0.0,
    "coord_diff_thresh": 0.0,
    "out_dir": OUT_DIR_PATH_OPT,
    "print_every": 10,
}


SEARCH_KW: Dict[str, Any] = {
    "max_depth": 10,
    "stitch_rmsd_thresh": 1.0e-4,
    "bridge_rmsd_thresh": 1.0e-4,
    "max_nodes_segment": 20,
    "max_nodes_bridge": 5,
    "kink_max_nodes": 3,
    "max_seq_kink": 2,
    "refine_mode": None,
}


IRC_KW: Dict[str, Any] = {
    "step_length": 0.10,
    "max_cycles": 125,
    "downhill": False,
    "forward": True,
    "backward": True,
    "root": 0,
    "hessian_init": "calc",
    "displ": "energy",
    "displ_energy": 1.0e-3,
    "displ_length": 0.10,
    "rms_grad_thresh": 1.0e-3,
    "hard_rms_grad_thresh": None,
    "energy_thresh": 1.0e-6,
    "energy_increase_thresh": 0.0,
    "imag_below": 0.0,
    "force_inflection": True,
    "check_bonds": False,
    "out_dir": OUT_DIR_IRC,
    "prefix": "",
    "hessian_update": "bofill",
    "hessian_recalc": None,
    # Low-level periodic checkpointing is disabled by default. The checkpoint
    # contains coordinates/energies/gradients only, never a dense Hessian.
    "dump_fn": "irc_data.h5",
    "dump_every": None,
    "max_pred_steps": 500,
    "loose_cycles": 3,
    "corr_func": "mbs",
    # Opt-in: bypass gradient and energy endpoint criteria and trace to
    # max_cycles. Numerical/integration failure and interruption still stop.
    "never_stop": False,
}


FREQ_KW: Dict[str, Any] = {
    "amplitude_ang": 0.8,
    "n_frames": 20,
    "max_write": 10,
    "sort": "value",
    "out_dir": OUT_DIR_FREQ,
}


THERMO_KW: Dict[str, Any] = {
    "temperature": 298.15,
    "pressure_atm": 1.0,
    "symmetry_number": 1,
    "dump": False,
}


TSOPT_MODE_ALIASES = (
    (("grad", "light", "dimer"), "dimer"),
    (("hess", "heavy", "rsprfo"), "rsprfo"),
    (("rsirfo",), "rsirfo"),
    (("trim",), "trim"),
)

# Saddle certification counts imaginary modes; it does not assess their
# character. Warn, without changing status, when the leading imaginary mode is
# very soft so the mode shape and IRC connectivity receive explicit review.
TS_IMAG_SOFT_WARN_CM = 50.0


DIMER_KW: Dict[str, Any] = {
    "length": 0.0189,
    "rotation_max_cycles": 15,
    "rotation_method": "fourier",
    "rotation_thresh": 1e-4,
    "rotation_tol": 1,
    "rotation_max_element": 0.001,
    "rotation_interpolate": True,
    "rotation_disable": False,
    "rotation_disable_pos_curv": True,
    "rotation_remove_trans": True,
    "trans_force_f_perp": True,
    "bonds": None,
    "N_hessian": None,
    "bias_rotation": False,
    "bias_translation": False,
    "bias_gaussian_dot": 0.1,
    "seed": None,
    "write_orientations": True,
    "forward_hessian": True,
}


HESSIAN_DIMER_KW: Dict[str, Any] = {
    "thresh_loose": "gau_loose",
    "thresh": "baker",
    "update_interval_hessian": 500,
    "neg_freq_thresh_cm": 5.0,
    "flatten_amp_ang": 0.10,
    "flatten_max_iter": 50,
    "flatten_sep_cutoff": 0.0,
    "flatten_k": 10,
    "flatten_loop_bofill": False,
    "mem": 100000,
    "device": "auto",
    "root": 0,
}

# LBFGS for TS inner loop (baker threshold)

LBFGS_TS_KW: Dict[str, Any] = {
    **LBFGS_KW,
    "thresh": "baker",
    # A dimer follows an effective saddle-search force; physical energy need
    # not decrease monotonically along that trajectory.
    "reject_uphill": False,
}

# Hessian Guided Dimer CLI-level defaults
#   (extends HESSIAN_DIMER_KW with nested dimer/lbfgs configs)

HESSIAN_DIMER_CLI_KW: Dict[str, Any] = {
    **HESSIAN_DIMER_KW,
    "dimer": {**DIMER_KW},
    "lbfgs": {**LBFGS_TS_KW},
}

# RFO-family shared defaults for TS optimization (hess/heavy → RS-P-RFO default; also RS-I-RFO / TRIM)

# Keys to drop from RFO_KW for TS optimizers: the gdiis/gediis family the base
# Optimizer rejects (it has no **kwargs), plus the reject_uphill group (accepted
# by TSHessianOptimizer but deliberately not forced on TS opts).
_RFO_ONLY_KEYS = {
    "gediis",
    "gdiis",
    "gdiis_thresh",
    "gediis_thresh",
    "gdiis_test_direction",
    "adapt_step_func",
    "rfo_overlaps",
    "reject_uphill",
    "uphill_tolerance",
    "rejection_trust_floor",
    "max_rejections_at_floor",
}

RSIRFO_KW: Dict[str, Any] = {
    **{k: v for k, v in RFO_KW.items() if k not in _RFO_ONLY_KEYS},
    "thresh": "baker",
    "check_eigval_structure": False,
    "trust_max": 0.10,
    "roots": [0],
    "hessian_ref": None,
    "rx_modes": None,
    "prim_coord": None,
    "rx_coords": None,
    "hessian_update": "bofill",
    "hessian_recalc_reset": True,
    "max_micro_cycles": 50,
    "augment_bonds": False,
    "assert_neg_eigval": False,
    "track_mode_by_overlap": False,
    # Keep the default on the standard restricted-step/root-following path.
    # Exact PHVA below remains the terminal saddle-order authority.
    "reject_mode_loss": False,
    "mode_loss_trust_floor": 1e-5,
    "max_mode_loss_rejections": 5,
    "verify_saddle": True,
    # Keep exact optimizer-side PHVA verification aligned with the final
    # frequency analysis (HESSIAN_DIMER_KW["neg_freq_thresh_cm"]).
    "saddle_imaginary_threshold_cm": 5.0,
    "saddle_recovery_step": 0.01,
    "saddle_recovery_check_interval": 50,
    "saddle_recovery_max_cycles": 0,
}

# Freq calc defaults (alias of UMA_CALC_KW)

FREQ_CALC_KW: Dict[str, Any] = dict(UMA_CALC_KW)
