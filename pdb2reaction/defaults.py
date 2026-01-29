# pdb2reaction/defaults.py

"""
Common default configurations for pdb2reaction workflows.

This module centralizes ALL shared settings used across opt, scan, freq, tsopt,
path_opt, path_search, and other command-line tools.

All default dictionaries are defined here to avoid redundant definitions across modules.
"""

from typing import Any, Dict

# -----------------------------------------------
# Geometry defaults
# -----------------------------------------------

GEOM_KW_DEFAULT: Dict[str, Any] = {
    "coord_type": "cart",       # coordinate representation (cart | dlc | redund)
    "freeze_atoms": [],         # list[int], 0-based atom indices to freeze
}

# -----------------------------------------------
# Calculator defaults (UMA)
# -----------------------------------------------

CALC_KW_DEFAULT: Dict[str, Any] = {
    # Charge and multiplicity
    "charge": 0,              # int, total charge
    "spin": 1,                # int, multiplicity (2S+1)

    # Model selection
    "model": "uma-s-1p1",     # str, UMA pretrained model ID
    "task_name": "omol",      # str, dataset/task tag

    # Device & graph construction
    "device": "auto",         # str, "cuda" | "cpu" | "auto"
    "max_neigh": None,        # Optional[int], max neighbors per atom
    "radius": None,           # Optional[float], cutoff radius in Å
    "r_edges": False,         # bool, use radial edges

    # Parallelism
    "workers": 1,             # int, UMA predictor workers (>1 spawns parallel predictor)
    "workers_per_node": 1,    # int, workers per node when workers>1

    # Hessian computation
    "hessian_calc_mode": "FiniteDifference",  # "Analytical" | "FiniteDifference"
    "out_hess_torch": False,  # bool, return torch.Tensor Hessian (CUDA) or numpy (CPU)
    "hessian_double": False,  # bool, use float64 for Hessian
    "return_partial_hessian": False,  # bool, return only active DOF
}

# -----------------------------------------------
# Optimizer base (common to LBFGS & RFO)
# -----------------------------------------------

OPT_BASE_KW: Dict[str, Any] = {
    # Convergence threshold preset
    "thresh": "gau",            # "gau_loose" | "gau" | "gau_tight" | "gau_vtight" | "baker" | "never"

    # Convergence criteria (forces in Hartree/bohr, steps in bohr)
    # +------------+------------------------------------------------------------+---------+--------+-----------+-------------+
    # |  Preset    | Purpose                                                    | max|F|  | RMS(F) | max|step| | RMS(step)   |
    # +------------+------------------------------------------------------------+---------+--------+-----------+-------------+
    # | gau_loose  | Loose/quick preoptimization; rough path searches           | 2.5e-3  | 1.7e-3 | 1.0e-2    | 6.7e-3      |
    # | gau        | Standard "Gaussian-like" tightness for routine work        | 4.5e-4  | 3.0e-4 | 1.8e-3    | 1.2e-3      |
    # | gau_tight  | Tighter; better structures / freq / TS refinement          | 1.5e-5  | 1.0e-5 | 6.0e-5    | 4.0e-5      |
    # | gau_vtight | Very tight; benchmarking/high-precision final structures   | 2.0e-6  | 1.0e-6 | 6.0e-6    | 4.0e-6      |
    # | baker*     | Baker-style rule (special; see below)                      | 3.0e-4* |   —    | 3.0e-4*   |     —       |
    # | never      | Disable built-in convergence (debug/external stopping)     |   —     |   —    |    —      |    —        |
    # +------------+------------------------------------------------------------+---------+--------+-----------+-------------+
    # * Baker rule: converged if (max|F| < 3.0e-4) AND (|ΔE| < 1.0e-6 OR max|step| < 3.0e-4).

    "max_cycles": 10000,         # hard cap on optimization cycles
    "print_every": 100,          # progress print frequency in cycles

    # Step-size safeguarding
    "min_step_norm": 1e-8,       # minimum ||step|| before raising ZeroStepLength
    "assert_min_step": True,     # enforce the min_step_norm check

    # Convergence criteria toggles
    "rms_force": None,           # Optional[float], if set, derive thresholds from this RMS(F)
    "rms_force_only": False,     # only check RMS(force)
    "max_force_only": False,     # only check max(|force|)
    "force_only": False,         # check RMS(force) and max(|force|) only

    # Extra convergence mechanisms
    "converge_to_geom_rms_thresh": 0.05,  # RMSD to reference geometry (Growing-NT)
    "overachieve_factor": 0.0,            # consider converged if forces < thresh/this_factor
    "check_eigval_structure": False,      # TS search: require expected negative modes

    # Line search
    "line_search": True,         # enable polynomial line search

    # Dumping / restart / bookkeeping
    "dump": False,               # write optimization trajectory
    "dump_restart": False,       # False | int, write restart YAML every N cycles (False disables)
    "prefix": "",                # file name prefix
    "out_dir": "./result_opt/",  # output directory
}

# -----------------------------------------------
# LBFGS-specific
# -----------------------------------------------

LBFGS_KW: Dict[str, Any] = {
    **OPT_BASE_KW,

    # History / memory
    "keep_last": 7,              # number of (s, y) pairs to retain

    # Preconditioner / initial scaling
    "beta": 1.0,                 # β in -(H + βI)^{-1} g
    "gamma_mult": False,         # estimate β from previous cycle (Nocedal Eq. 7.20)

    # Step-size control
    "max_step": 0.30,            # maximum allowed component-wise step
    "control_step": True,        # scale step to satisfy |max component| <= max_step

    # Safeguards
    "double_damp": True,         # double-damping to enforce s·y > 0

    # Regularized L-BFGS (μ_reg)
    "mu_reg": None,              # initial regularization; enables regularized L-BFGS if set
    "max_mu_reg_adaptions": 10,  # maximum trial steps for μ adaptation
}

# -----------------------------------------------
# RFO-specific
# -----------------------------------------------

RFO_KW: Dict[str, Any] = {
    **OPT_BASE_KW,

    # Trust-region (step-size) control
    "trust_radius": 0.10,        # initial trust radius (in working coordinates)
    "trust_update": True,        # adapt the trust radius based on step quality
    "trust_min": 0.00,           # lower bound for trust radius
    "trust_max": 0.10,           # upper bound for trust radius
    "max_energy_incr": None,     # abort if ΔE exceeds this after a bad step

    # Hessian model / refresh
    "hessian_update": "bfgs",    # "bfgs" (faster convergence) | "bofill" (more robust)
    "hessian_init": "calc",      # initial Hessian calculation
    "hessian_recalc": 200,       # recompute exact Hessian every N cycles
    "hessian_recalc_adapt": None,# heuristic: trigger exact Hessian recompute based on force norm

    # Numerical hygiene & mode filtering
    "small_eigval_thresh": 1e-8, # treat |λ| < threshold as zero / remove corresponding modes

    # RFO/RS micro-iterations
    "alpha0": 1.0,               # initial α for restricted-step RFO
    "max_micro_cycles": 50,      # max inner iterations to hit the trust radius
    "rfo_overlaps": False,       # mode following via eigenvector overlap across cycles

    # Inter/Extrapolation helpers
    "gediis": False,             # enable GEDIIS (energy-based DIIS)
    "gdiis": True,               # enable GDIIS (gradient-based DIIS)

    # Thresholds for enabling DIIS (semantics matter)
    "gdiis_thresh": 2.5e-3,      # compared to RMS(step)  → enable GDIIS when small
    "gediis_thresh": 1.0e-2,     # compared to RMS(force) → enable GEDIIS when small

    "gdiis_test_direction": True,# compare DIIS step direction to the RFO step

    # Choice of step model
    "adapt_step_func": True,     # switch to shifted-Newton on trust when PD & gradient is small
}

# -----------------------------------------------
# Bias (harmonic well) defaults
# -----------------------------------------------

BIAS_KW: Dict[str, Any] = {
    "k": 100,  # float, harmonic bias strength in eV/Å^2
}

# -----------------------------------------------
# Bond-change detection
# -----------------------------------------------

BOND_KW: Dict[str, Any] = {
    "device": "cuda",            # str, device for UMA graph analysis during bond detection
    "bond_factor": 1.20,         # float, scaling of covalent radii for bond cutoff
    "margin_fraction": 0.05,     # float, fractional margin to tolerate small deviations
    "delta_fraction": 0.05,      # float, change threshold to flag bond formation/breaking
}

# -----------------------------------------------
# Optimizer mode aliases
# -----------------------------------------------

OPT_MODE_ALIASES = (
    (("light",), "lbfgs"),
    (("heavy",), "rfo"),
)

# -----------------------------------------------
# UMA Calculator extended defaults (uma_pysis)
# -----------------------------------------------

# Extended UMA calculator defaults with additional keys for Hessian control
UMA_CALC_KW: Dict[str, Any] = {
    **CALC_KW_DEFAULT,
    "out_hess_torch": True,   # bool, return Hessian as torch.Tensor (CUDA) or numpy (CPU)
    "freeze_atoms": None,     # Optional[Sequence[int]], list of 0-based freeze atom indices
}

# -----------------------------------------------
# Scan-specific defaults (scan, scan2d, scan3d)
# -----------------------------------------------

# Scan optimizer defaults (used by scan.py, scan2d.py, scan3d.py)
SCAN_OPT_BASE_KW: Dict[str, Any] = {
    **OPT_BASE_KW,
    "dump": False,
    "max_cycles": 10000,
}

SCAN_LBFGS_KW: Dict[str, Any] = {
    **LBFGS_KW,
}

SCAN_RFO_KW: Dict[str, Any] = {
    **RFO_KW,
}
