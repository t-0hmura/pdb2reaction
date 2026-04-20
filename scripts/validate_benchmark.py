#!/usr/bin/env python3
"""
Validation summary script for p2r cluster model calculations.

For each completed run, checks:
1. Did MEP find a reaction path? (bond changes detected)
2. Did TS optimization succeed? (imaginary frequency present)
3. Does IRC reproduce the correct reaction? (bond changes match input R→P)
4. Are optimized R/P structures consistent with input? (no spurious bond changes)
5. Energy barriers vs literature reference

Output: CSV table + human-readable summary

Usage:
    python validate_benchmark.py [--backend uma_s1p1] [--csv results.csv]
"""

import argparse
import csv
import os
import re
import sys
from pathlib import Path

import json

import numpy as np
import yaml

# Bond change detection API (direct import, no subprocess)
from pdb2reaction.bond_changes import compare_structures, summarize_changes
from pysisyphus.helpers import geom_loader

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
# Default path points at the shipped benchmark inputs under
# <repo>/examples/benchmark. Override with PDB2REACTION_BENCHMARK_DIR if
# you keep results in a separate directory. Note: the original validation
# layout expected ``<base>/p2r_cluster_<backend>/<pdb>/<model>``; if your
# own runs sit directly under ``<base>/<pdb>/<model>`` (as the shipped
# example inputs do) there is nothing for this validator to analyse until
# you populate result_* directories by executing each task's run.sh.
VALIDATION_DIR = Path(
    os.environ.get(
        "PDB2REACTION_BENCHMARK_DIR",
        str(Path(__file__).resolve().parent.parent / "examples" / "benchmark"),
    )
)
BACKENDS = ["uma_s1p1", "uma_s1p2", "uma_m1p1", "mace", "orb"]
PDBS = ["1AH7", "1O9S", "1PWZ", "1RTQ", "2E7Z", "4OTA"]

# Condition label -> result directory name
CONDITIONS = {
    "default": "result",
    "nrp": "result_nrp",
    "tsonly": "result_tsonly",
    "ana": "result_ana",
}

# Literature reference barriers (B3LYP single-point, kcal/mol)
# Derived from dataset/Cluster/<pdb>/energy_profile.csv with the LEFT input
# structure (as fed to pdb2reaction via -i) as the baseline (0 kcal/mol).
# barrier = max(middle_states) - E(left_input)
# rxn_e   = E(right_input) - E(left_input)
LITERATURE = {
    # 1AH7: R -> TS -> P (single step per model)
    ("1AH7", "model_A"):  {"barrier": 22.8, "rxn_e": -26.1, "ref": "Liao 2010 JPC"},
    ("1AH7", "model_B"):  {"barrier": 20.9, "rxn_e": 1.5,   "ref": "Liao 2010 JPC"},
    ("1AH7", "model_Ba"): {"barrier": 3.9,  "rxn_e": -7.1,  "ref": "Liao 2010 JPC"},
    # 1O9S: single step per model (R -> TS -> P)
    ("1O9S", "model_A"):  {"barrier": 18.8, "rxn_e": -2.9,  "ref": "Georgieva 2010 JCC"},
    ("1O9S", "model_B"):  {"barrier": 21.7, "rxn_e": 0.5,   "ref": "Georgieva 2010 JCC"},
    ("1O9S", "model_C"):  {"barrier": 15.4, "rxn_e": -16.7, "ref": "Georgieva 2010 JCC"},
    ("1O9S", "model_D"):  {"barrier": 18.9, "rxn_e": -9.2,  "ref": "Georgieva 2010 JCC"},
    # 1PWZ: single step per model (R -> TS -> P)
    ("1PWZ", "model_A"):  {"barrier": 23.0, "rxn_e": 17.5,  "ref": "Hopmann 2008 JCTC"},
    ("1PWZ", "model_B"):  {"barrier": 17.9, "rxn_e": 14.1,  "ref": "Hopmann 2008 JCTC"},
    ("1PWZ", "model_C"):  {"barrier": 18.2, "rxn_e": 5.5,   "ref": "Hopmann 2008 JCTC"},
    # 1RTQ: R -> TS1 -> IM1_1 -> IM1_2 -> TS2 -> IM2 -> TS3 -> IM3 -> TS4 -> P
    # Each step is scored under the left-input baseline convention documented in
    # paper SI-G: barrier = max(middle stationary points) - left input,
    # rxn_e = right input - left input.
    ("1RTQ", "model/step1"): {"barrier": 12.1, "rxn_e": 11.3,  "ref": "Chen 2008 JPC"},
    ("1RTQ", "model/step2"): {"barrier": 4.2,  "rxn_e": -1.2,  "ref": "Chen 2008 JPC"},
    ("1RTQ", "model/step3"): {"barrier": -0.8, "rxn_e": -0.1,  "ref": "Chen 2008 JPC"},
    ("1RTQ", "model/step4"): {"barrier": -0.6, "rxn_e": -11.7, "ref": "Chen 2008 JPC"},
    # 2E7Z/model_A: R -> IM1 -> TS1 -> IM2 -> TS2 -> IM3 -> TS3 -> IM4 -> TS4 -> P
    # step1 IM1(-3.9)/IM2(13.9), step2 IM2/IM3(-13.1), step3 IM3/IM4(-20.6), step4 IM4/P(-17)
    ("2E7Z", "model_A/step1"): {"barrier": 19.8, "rxn_e": 17.8,  "ref": "Liao 2010 PNAS"},
    ("2E7Z", "model_A/step2"): {"barrier": 5.5,  "rxn_e": -27.0, "ref": "Liao 2010 PNAS"},
    ("2E7Z", "model_A/step3"): {"barrier": 9.5,  "rxn_e": -7.5,  "ref": "Liao 2010 PNAS"},
    ("2E7Z", "model_A/step4"): {"barrier": 14.7, "rxn_e": 3.6,   "ref": "Liao 2010 PNAS"},
    # 2E7Z/model_B: neutral-Asp13 Model B single-TS run; left input is the free
    # reactant R, so under the left-input baseline convention (SI-G) the
    # literature barrier is E_TS1 - E_R = 13.7 kcal/mol and the reaction energy
    # R -> IM2 is 0.2 kcal/mol. The alternative IM1-referenced barrier (33.4
    # kcal/mol) from the Kromann 2016 repository is not used here because IM1
    # is not the left input of the pdb2reaction run for this system.
    ("2E7Z", "model_B"):  {"barrier": 13.7, "rxn_e": 0.2,  "ref": "Liao 2010 PNAS"},
    # 4OTA: R -> TS1 -> IM -> TS2 -> P (2 steps per model)
    ("4OTA", "model_A/step1"): {"barrier": 12.8, "rxn_e": 9.8,   "ref": "Sevastik 2007"},
    ("4OTA", "model_A/step2"): {"barrier": -2.8, "rxn_e": -13.5, "ref": "Sevastik 2007"},
    ("4OTA", "model_B/step1"): {"barrier": 6.9,  "rxn_e": -5.7,  "ref": "Sevastik 2007"},
    ("4OTA", "model_B/step2"): {"barrier": 4.1,  "rxn_e": 2.8,   "ref": "Sevastik 2007"},
}

# Experimental barriers where available (kcal/mol)
EXPERIMENTAL = {
    "1O9S": {"barrier": 20.9, "ref": "Georgieva 2010 (cited kcat)"},
}


# ---------------------------------------------------------------------------
# Timing extraction
# ---------------------------------------------------------------------------
_RE_STAGE_TIME = re.compile(
    r"\[time\] Elapsed Time for (.+?):\s+(\d+):(\d+):(\d+)\.(\d+)"
)
_RE_TOTAL_TIME = re.compile(
    r"\[all\] Elapsed for Whole Pipeline:\s+(\d+):(\d+):(\d+)\.(\d+)"
)


def _hms_to_sec(h, m, s, ms):
    return int(h) * 3600 + int(m) * 60 + int(s) + int(ms) / 1000.0


def parse_timing(all_out_path):
    """Extract timing data from all.out log file.

    Returns dict with keys: total_sec, path_opt_sec, ts_opt_sec, irc_sec,
    freq_total_sec, dft_total_sec. DFT time is summed separately from freq
    (each `[time] Elapsed Time for DFT` line is accumulated).
    """
    result = {
        "total_sec": None,
        "path_opt_sec": None,
        "ts_opt_sec": None,
        "irc_sec": None,
        "freq_total_sec": None,
        "dft_total_sec": None,
    }
    try:
        text = Path(all_out_path).read_text()
    except Exception:
        return result

    freq_times = []
    dft_times = []
    for line in text.splitlines():
        m = _RE_STAGE_TIME.search(line)
        if m:
            stage = m.group(1)
            sec = _hms_to_sec(m.group(2), m.group(3), m.group(4), m.group(5))
            if stage == "Path Opt":
                result["path_opt_sec"] = sec
            elif stage == "TS Opt":
                result["ts_opt_sec"] = sec
            elif stage == "IRC":
                result["irc_sec"] = sec
            elif stage == "Freq":
                freq_times.append(sec)
            elif stage == "DFT":
                dft_times.append(sec)
            continue
        m = _RE_TOTAL_TIME.search(line)
        if m:
            result["total_sec"] = _hms_to_sec(m.group(1), m.group(2), m.group(3), m.group(4))

    if freq_times:
        result["freq_total_sec"] = sum(freq_times)
    if dft_times:
        result["dft_total_sec"] = sum(dft_times)

    return result


def fmt_time(seconds):
    """Format seconds as MM:SS or HH:MM:SS string, or '—'."""
    if seconds is None:
        return "—"
    seconds = float(seconds)
    h = int(seconds // 3600)
    m = int((seconds % 3600) // 60)
    s = int(seconds % 60)
    if h > 0:
        return f"{h}:{m:02d}:{s:02d}"
    return f"{m}:{s:02d}"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def n_bond_changes(bc):
    """Total number of bond changes (formed + broken)."""
    if bc is None:
        return 999  # large number to avoid selection
    return len(bc["formed"]) + len(bc["broken"])


# Metals and elements whose coordination changes are considered "soft"
METALS = {"Zn", "W", "Cu", "Fe", "Mg", "Mn", "Co", "Ni", "Ca"}


def is_metal_coord_change(bond_str):
    """Check if a bond change involves a metal coordination (soft change)."""
    parts = bond_str.split("-")
    if len(parts) != 2:
        return False
    # Extract element from e.g. "Zn47" -> "Zn"
    elem1 = "".join(c for c in parts[0] if c.isalpha())
    elem2 = "".join(c for c in parts[1] if c.isalpha())
    return elem1 in METALS or elem2 in METALS


def classify_opt_changes(bc):
    """Classify bond changes from geometry optimization.
    Returns: 'OK', 'SOFT' (only metal coord changes), or 'BROKEN' (covalent bond changed).
    """
    if bc is None:
        return "N/A"
    if bond_changes_empty(bc):
        return "OK"
    all_changes = bc["formed"] + bc["broken"]
    n_hard = sum(1 for b in all_changes if not is_metal_coord_change(b))
    if n_hard == 0:
        return "SOFT"  # only metal coordination changes
    return "BROKEN"


def filter_nonmetal(bc):
    """Filter bond changes to non-metal bonds only."""
    if bc is None:
        return None
    return {
        "formed": [b for b in bc["formed"] if not is_metal_coord_change(b)],
        "broken": [b for b in bc["broken"] if not is_metal_coord_change(b)],
    }


def judge_reaction(bc_input_nm, bc_irc_nm):
    """Judge IRC reproduction of expected reaction (non-metal bonds).
    ✓ = all expected non-metal bond changes found in IRC output
    ~ = some but not all found
    ✗ = none found
    Returns (judge, detail).
    """
    if bc_input_nm is None or bc_irc_nm is None:
        return "—", "PARSE_ERR"
    ess_f = set(bc_input_nm["formed"])
    ess_b = set(bc_input_nm["broken"])
    irc_f = set(bc_irc_nm["formed"])
    irc_b = set(bc_irc_nm["broken"])
    f_ok = ess_f.issubset(irc_f)
    b_ok = ess_b.issubset(irc_b)
    if f_ok and b_ok:
        return "✓", "ALL_MATCH"
    elif (ess_f & irc_f) or (ess_b & irc_b):
        return "~", "PARTIAL"
    else:
        return "✗", "NONE"


def run_bond_summary(f1, f2):
    """Detect bond changes between two XYZ files using the API directly."""
    if not os.path.isfile(f1) or not os.path.isfile(f2):
        return None
    try:
        geom1 = geom_loader(f1)
        geom2 = geom_loader(f2)
        result = compare_structures(geom1, geom2, device="cpu")
        elems = [a.capitalize() for a in geom1.atoms]
        formed = []
        for i, j in sorted(result.formed_covalent):
            formed.append(f"{elems[i]}{i+1}-{elems[j]}{j+1}")
        broken = []
        for i, j in sorted(result.broken_covalent):
            broken.append(f"{elems[i]}{i+1}-{elems[j]}{j+1}")
        return {"formed": sorted(formed), "broken": sorted(broken)}
    except Exception:
        return None


def bond_changes_match(bc1, bc2):
    """Check if two bond-change dicts describe the same reaction."""
    if bc1 is None or bc2 is None:
        return None  # can't determine
    return set(bc1["formed"]) == set(bc2["formed"]) and set(bc1["broken"]) == set(bc2["broken"])


def bond_changes_empty(bc):
    """True if no bond changes detected."""
    if bc is None:
        return None
    return len(bc["formed"]) == 0 and len(bc["broken"]) == 0


def fmt_bc(bc):
    """Format bond changes as short string."""
    if bc is None:
        return "N/A"
    parts = []
    for b in bc["formed"]:
        parts.append(f"+{b}")
    for b in bc["broken"]:
        parts.append(f"-{b}")
    return " ".join(parts) if parts else "(none)"


def parse_run_sh_endpoints(run_dir, result_dir_name="result"):
    """Parse run.sh to extract the -i arguments for the given condition's run."""
    run_sh = Path(run_dir) / "run.sh"
    if not run_sh.exists():
        return None, None
    text = run_sh.read_text()
    for line in text.splitlines():
        stripped = line.lstrip("#").strip()
        if "pdb2reaction" not in stripped:
            continue
        # Must target this exact --out-dir (anchor to end-of-token to avoid prefix match)
        if not re.search(rf"--out-dir\s+{re.escape(result_dir_name)}(\s|$)", stripped):
            continue
        # Extract -i arguments (1 or 2 files — tsonly uses 1 file)
        m = re.search(r"-i\s+(\S+\.(?:xyz|pdb))(?:\s+(\S+\.(?:xyz|pdb)))?", stripped)
        if m:
            r_file = str(Path(run_dir) / m.group(1))
            p_file = str(Path(run_dir) / m.group(2)) if m.group(2) else None
            return r_file, p_file
    return None, None


def find_input_endpoints(run_dir, result_dir_name="result"):
    """Find input R and P files. Prefers run.sh parsing for correct multi-step endpoints."""
    # Primary: parse run.sh
    r, p = parse_run_sh_endpoints(run_dir, result_dir_name)
    if r and p and os.path.isfile(r) and os.path.isfile(p):
        return r, p
    if r and os.path.isfile(r) and not p:
        # tsonly case: single TS input; no paired P. Return (r, None).
        return r, None
    # Fallback: first/last sorted xyz
    run_dir = Path(run_dir)
    xyz_files = sorted(run_dir.glob("*.xyz"))
    if not xyz_files:
        pdb_files = sorted(run_dir.glob("*.pdb"))
        if pdb_files:
            return str(pdb_files[0]), str(pdb_files[-1])
        return None, None
    return str(xyz_files[0]), str(xyz_files[-1])


def _post_seg_parent(run_dir, result_dir_name):
    if result_dir_name == "result_tsonly":
        return None
    if result_dir_name in ("result_nrp", "result_ana"):
        return Path(run_dir) / result_dir_name / "path_opt"
    return Path(run_dir) / result_dir_name / "path_search"


def list_post_seg_dirs(run_dir, result_dir_name):
    """Return all post_seg_NN dirs sorted by NN. Empty list if none/tsonly."""
    parent = _post_seg_parent(run_dir, result_dir_name)
    if parent is None or not parent.exists():
        return []
    return sorted(parent.glob("post_seg_*"), key=lambda p: p.name)


def get_seg_dir(run_dir, result_dir_name):
    """Return the post-segment directory for this condition.

    For multi-segment paths, returns the LAST post_seg_NN. For tsonly,
    returns tsopt_single. Falls back to post_seg_01 even if missing.
    """
    if result_dir_name == "result_tsonly":
        return Path(run_dir) / result_dir_name / "tsopt_single"
    segs = list_post_seg_dirs(run_dir, result_dir_name)
    if segs:
        return segs[-1]
    return _post_seg_parent(run_dir, result_dir_name) / "post_seg_01"


def _struct_or_irc_endpoints(seg_dir):
    """Pick reactant/product from a single post_seg dir (structures/ or irc/)."""
    struct_dir = seg_dir / "structures"
    r_opt = struct_dir / "reactant.xyz"
    p_opt = struct_dir / "product.xyz"
    if r_opt.exists() and p_opt.exists():
        return str(r_opt), str(p_opt)
    irc_dir = seg_dir / "irc"
    first = irc_dir / "finished_first.xyz"
    last = irc_dir / "finished_last.xyz"
    if first.exists() and last.exists():
        return str(first), str(last)
    return None, None


def find_irc_endpoints(run_dir, result_dir_name="result"):
    """Find IRC-optimized R/P. For multi-seg paths uses first seg's R and last seg's P."""
    if result_dir_name == "result_tsonly":
        seg_dir = Path(run_dir) / result_dir_name / "tsopt_single"
        r, p = _struct_or_irc_endpoints(seg_dir)
        return r, p, (r is not None and p is not None)

    segs = list_post_seg_dirs(run_dir, result_dir_name)
    if not segs:
        return None, None, False
    first_r, _ = _struct_or_irc_endpoints(segs[0])
    _, last_p = _struct_or_irc_endpoints(segs[-1])
    if first_r and last_p:
        return first_r, last_p, True
    # Fallback: try the last seg alone
    r, p = _struct_or_irc_endpoints(segs[-1])
    if r and p:
        return r, p, True
    return None, None, False


def assign_irc_to_input(inp_R, inp_P, irc_first, irc_last):
    """
    Assign IRC endpoints to input R/P based on bond changes.
    The IRC endpoint with fewer bond changes vs input R is assigned as IRC_R.
    Returns (irc_R, irc_P, swapped).
    """
    bc_R_first = run_bond_summary(inp_R, irc_first)
    bc_R_last = run_bond_summary(inp_R, irc_last)

    # Cost = number of bond changes (fewer = more similar)
    cost_normal = n_bond_changes(bc_R_first) + n_bond_changes(run_bond_summary(inp_P, irc_last))
    cost_swap = n_bond_changes(bc_R_last) + n_bond_changes(run_bond_summary(inp_P, irc_first))

    if cost_normal <= cost_swap:
        return irc_first, irc_last, False
    else:
        return irc_last, irc_first, True


def get_barrier_from_summary(summary_path, diagram_name="energy_diagram_G_UMA"):
    """Extract barrier from summary.json energy diagram."""
    if not os.path.isfile(summary_path):
        return None, None
    with open(summary_path) as f:
        data = json.load(f)
    if not data:
        return None, None
    for diag in data.get("energy_diagrams", []):
        if diag.get("name") == diagram_name:
            energies = diag.get("energies_kcal", [])
            if len(energies) >= 2:
                barrier = max(energies[1:]) - energies[0]
                rxn_e = energies[-1] - energies[0]
                return barrier, rxn_e
    # Fallback to energy_diagram_UMA
    for diag in data.get("energy_diagrams", []):
        if diag.get("name") == "energy_diagram_UMA":
            energies = diag.get("energies_kcal", [])
            if len(energies) >= 2:
                barrier = max(energies[1:]) - energies[0]
                rxn_e = energies[-1] - energies[0]
                return barrier, rxn_e
    # Fallback to MEP barrier
    segs = data.get("segments", [])
    if segs:
        barrier = segs[0].get("barrier_kcal")
        delta = segs[0].get("delta_kcal")
        return barrier, delta
    return None, None


def check_ts_imag_freq(run_dir, result_dir_name="result"):
    """Count TS imaginary modes by inspecting <seg_dir>/ts/vib/imag_*cm-1*.xyz.

    Each pysisyphus-emitted file corresponds to one imaginary mode; pysisyphus
    already applies a 5 cm^-1 cutoff upstream, so every emitted file is
    counted here without further filtering.

    Returns the count, or None if the vib directory does not exist.
    """
    seg_dir = get_seg_dir(run_dir, result_dir_name)
    vib_dir = seg_dir / "ts" / "vib"
    if not vib_dir.exists():
        return None
    return sum(1 for _ in vib_dir.glob("imag_*cm-1*.xyz"))


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------
def find_rp_inputs_for_tsonly(run_dir):
    """For tsonly, the run.sh only specifies a single TS file.
    Use the R/P input files that appear in the 'default' pdb2reaction all command
    (--out-dir result) to compare against IRC-relaxed endpoints.
    """
    r, p = parse_run_sh_endpoints(run_dir, "result")
    if r and p and os.path.isfile(r) and os.path.isfile(p):
        return r, p
    # Fallback: first/last sorted xyz from run dir
    run_dir = Path(run_dir)
    xyz_files = sorted(run_dir.glob("*.xyz"))
    if xyz_files:
        return str(xyz_files[0]), str(xyz_files[-1])
    return None, None


def analyze_run(backend, pdb, model_key, run_dir, condition="default"):
    """Analyze a single run directory under a given result condition.

    condition: one of CONDITIONS keys. Selects which result_* subdir to validate.
    """
    result_dir_name = CONDITIONS.get(condition, "result")
    out_file_name = {
        "default": "all.out",
        "nrp": "all_nrp.out",
        "tsonly": "all_tsonly.out",
        "ana": "all_ana.out",
    }.get(condition, "all.out")

    result = {
        "backend": backend,
        "pdb": pdb,
        "model": model_key,
        "condition": condition,
        "status": "UNKNOWN",
        "judge": "N/A",
        "opt_r": "N/A",
        "opt_p": "N/A",
        "opt_r_bc": "N/A",
        "opt_p_bc": "N/A",
        "mep_reaction": "N/A",
        "irc_reaction_match": "N/A",
        "irc_swapped": "N/A",
        "r_shifted": "N/A",
        "p_shifted": "N/A",
        "ts_n_imag": "N/A",
        "barrier_kcal": "N/A",
        "rxn_e_kcal": "N/A",
        "lit_barrier": "N/A",
        "barrier_error": "N/A",
        "input_bc": "N/A",
        "irc_bc": "N/A",
        "r_shift_bc": "N/A",
        "p_shift_bc": "N/A",
        "total_time_sec": "N/A",
        "path_opt_sec": "N/A",
        "ts_opt_sec": "N/A",
        "irc_sec": "N/A",
        "freq_total_sec": "N/A",
        "dft_total_sec": "N/A",
    }

    run_dir = Path(run_dir)
    if not run_dir.exists():
        result["status"] = "DIR_MISSING"
        return result

    # --- Timing extraction (always attempt, even for failed runs) ---
    all_out = run_dir / out_file_name
    if all_out.exists():
        timing = parse_timing(str(all_out))
        if timing["total_sec"] is not None:
            result["total_time_sec"] = f"{timing['total_sec']:.1f}"
        for key in ("path_opt_sec", "ts_opt_sec", "irc_sec", "freq_total_sec", "dft_total_sec"):
            if timing[key] is not None:
                result[key] = f"{timing[key]:.1f}"

    # --- Locate summary.json for this condition ---
    result_dir = run_dir / result_dir_name
    if condition == "tsonly":
        summary_json = result_dir / "tsopt_single" / "summary.json"
        if not summary_json.exists():
            summary_json = result_dir / "summary.json"
    else:
        summary_json = result_dir / "summary.json"
        if not summary_json.exists():
            summary_json = result_dir / "path_opt" / "summary.json"

    if not summary_json.exists():
        # Check for error in .out file
        if all_out.exists():
            text = all_out.read_text()[-500:]
            if "Error" in text or "Traceback" in text:
                result["status"] = "ERROR"
            elif "SVD did not converge" in text:
                result["status"] = "SVD_FAIL"
            else:
                result["status"] = "RUNNING_OR_INCOMPLETE"
        else:
            result["status"] = "NOT_STARTED"
        return result

    # --- MEP bond changes (from summary.json segments) ---
    with open(summary_json) as f:
        data = json.load(f)
    segs = data.get("segments", []) if data else []
    if segs:
        bc_list = segs[0].get("bond_changes", [])
        has_bc = False
        for item in bc_list:
            if isinstance(item, dict):
                for k, v in item.items():
                    if isinstance(v, list) and v:
                        has_bc = True
        result["mep_reaction"] = "YES" if has_bc else "NO_BC"
    else:
        result["mep_reaction"] = "NO_SEG"

    # --- Cluster systems are all single-step reactions ---
    # If the pipeline split the path into 2+ segments, mark as failure.
    # (tsonly never splits — skip check.)
    if condition != "tsonly" and len(segs) > 1:
        result["status"] = "MULTI_SEGMENT"
        result["judge"] = "✗"
        result["irc_reaction_match"] = f"SPLIT_INTO_{len(segs)}_SEGMENTS"
        return result

    # --- Barrier ---
    barrier, rxn_e = get_barrier_from_summary(str(summary_json))
    if barrier is not None:
        result["barrier_kcal"] = f"{barrier:.1f}"
        if abs(barrier) > 200 or (rxn_e is not None and abs(rxn_e) > 500):
            result["status"] = "ANOMALOUS"
            result["rxn_e_kcal"] = f"{rxn_e:.1f}" if rxn_e else "N/A"
            return result
    if rxn_e is not None:
        result["rxn_e_kcal"] = f"{rxn_e:.1f}"

    # Literature comparison
    lit_key = (pdb, model_key)
    if lit_key in LITERATURE:
        lit = LITERATURE[lit_key]
        result["lit_barrier"] = f"{lit['barrier']:.1f}"
        if barrier is not None and lit["barrier"] is not None:
            result["barrier_error"] = f"{barrier - lit['barrier']:+.1f}"

    # --- Input endpoints ---
    # For tsonly the run.sh line has a single TS file; compare IRC to R/P
    # from the default run's -i args (living in the same run dir).
    if condition == "tsonly":
        inp_R, inp_P = find_rp_inputs_for_tsonly(run_dir)
    else:
        inp_R, inp_P = find_input_endpoints(run_dir, result_dir_name)
    if not inp_R or not inp_P:
        result["status"] = "NO_INPUT"
        return result

    # Input R→P bond changes
    bc_input = run_bond_summary(inp_R, inp_P)
    result["input_bc"] = fmt_bc(bc_input)

    # --- Geometry optimization check (init00 = R, init01 = P) ---
    # tsonly has no initXX_lbfgs_opt (no preopt of endpoints)
    if condition != "tsonly":
        path_opt_dir = run_dir / result_dir_name / "path_opt"
        opt_R = path_opt_dir / "init00_lbfgs_opt" / "final_geometry.xyz"
        opt_P = path_opt_dir / "init01_lbfgs_opt" / "final_geometry.xyz"
        if opt_R.exists():
            bc_opt_r = run_bond_summary(inp_R, str(opt_R))
            result["opt_r_bc"] = fmt_bc(bc_opt_r)
            result["opt_r"] = classify_opt_changes(bc_opt_r)
        if opt_P.exists():
            bc_opt_p = run_bond_summary(inp_P, str(opt_P))
            result["opt_p_bc"] = fmt_bc(bc_opt_p)
            result["opt_p"] = classify_opt_changes(bc_opt_p)

    # --- IRC endpoints ---
    irc_first, irc_last, has_irc = find_irc_endpoints(run_dir, result_dir_name)
    if not has_irc:
        result["irc_reaction_match"] = "NO_IRC"
        result["status"] = "OK_NO_IRC"
        return result

    # Assign IRC endpoints to input R/P by RMSD (2 orderings)
    try:
        irc_R, irc_P, swapped = assign_irc_to_input(inp_R, inp_P, irc_first, irc_last)
        result["irc_swapped"] = "YES" if swapped else "NO"
    except Exception as e:
        result["irc_reaction_match"] = f"RMSD_ERR: {e}"
        result["status"] = "IRC_ASSIGN_ERR"
        return result

    # IRC R→P bond changes (with correct assignment)
    bc_irc = run_bond_summary(irc_R, irc_P)
    result["irc_bc"] = fmt_bc(bc_irc)

    # Match check (non-metal subset matching)
    bc_input_nm = filter_nonmetal(bc_input)
    bc_irc_nm = filter_nonmetal(bc_irc)
    judge, detail = judge_reaction(bc_input_nm, bc_irc_nm)
    result["judge"] = judge
    result["irc_reaction_match"] = detail

    # R shift check (Input R → IRC R, non-metal only)
    bc_r_shift = run_bond_summary(inp_R, irc_R)
    bc_r_shift_nm = filter_nonmetal(bc_r_shift)
    result["r_shift_bc"] = fmt_bc(bc_r_shift_nm)
    if bc_r_shift_nm:
        result["r_shifted"] = "OK" if bond_changes_empty(bc_r_shift_nm) else "SHIFTED"

    # P shift check (Input P → IRC P, non-metal only)
    bc_p_shift = run_bond_summary(inp_P, irc_P)
    bc_p_shift_nm = filter_nonmetal(bc_p_shift)
    result["p_shift_bc"] = fmt_bc(bc_p_shift_nm)
    if bc_p_shift_nm:
        result["p_shifted"] = "OK" if bond_changes_empty(bc_p_shift_nm) else "SHIFTED"

    # --- TS imaginary frequency (informational; does not affect judge) ---
    n_imag = check_ts_imag_freq(run_dir, result_dir_name)
    if n_imag is not None:
        result["ts_n_imag"] = str(n_imag)

    # --- Final status (based on non-metal judgment) ---
    opt_broken = result["opt_r"] == "BROKEN" or result["opt_p"] == "BROKEN"

    if opt_broken:
        result["status"] = "OPT_BROKEN"
    elif judge == "✓" and result["r_shifted"] == "OK" and result["p_shifted"] == "OK":
        result["status"] = "PERFECT"
    elif judge == "✓":
        result["status"] = "GOOD"
    elif judge == "~":
        result["status"] = "PARTIAL_REACT"
    elif judge == "✗":
        result["status"] = "WRONG_REACT"
    elif judge == "—":
        result["status"] = "UNCERTAIN"

    return result


def main():
    parser = argparse.ArgumentParser(description="Validate p2r cluster model results")
    parser.add_argument("--backend", nargs="*", default=BACKENDS, help="Backends to check")
    parser.add_argument("--pdb", nargs="*", default=PDBS, help="PDB systems to check")
    parser.add_argument("--condition", nargs="*", default=list(CONDITIONS.keys()),
                        help="Result conditions: default nrp tsonly ana")
    parser.add_argument("--csv", default=None, help="Output CSV path")
    parser.add_argument("--summary", action="store_true", help="Print summary table only")
    args = parser.parse_args()

    all_results = []

    for backend in args.backend:
        base = VALIDATION_DIR / f"p2r_cluster_{backend}"
        if not base.exists():
            continue
        for pdb in args.pdb:
            pdb_dir = base / pdb
            if not pdb_dir.exists():
                continue
            for model_name in sorted(os.listdir(pdb_dir)):
                model_dir = pdb_dir / model_name
                if not model_dir.is_dir():
                    continue

                step_dirs = sorted([d for d in model_dir.iterdir() if d.is_dir() and d.name.startswith("step")])

                if step_dirs:
                    for step_dir in step_dirs:
                        model_key = f"{model_name}/{step_dir.name}"
                        for cond in args.condition:
                            r = analyze_run(backend, pdb, model_key, step_dir, condition=cond)
                            all_results.append(r)
                            if not args.summary:
                                status_icon = {"PERFECT": "✓✓", "GOOD": "✓", "PARTIAL_REACT": "△",
                                               "WRONG_REACT": "✗", "IRC_FLAT": "~", "ANOMALOUS": "!",
                                               "SVD_FAIL": "F", "ERROR": "E"}.get(r["status"], "?")
                                print(f"{status_icon} {backend:14s} {pdb}/{model_key:20s} [{cond:7s}] "
                                      f"status={r['status']:16s} barrier={r['barrier_kcal']:>8s} "
                                      f"lit={r['lit_barrier']:>6s} err={r['barrier_error']:>7s} "
                                      f"optR={r['opt_r']:7s} optP={r['opt_p']:7s} "
                                      f"irc={r['irc_reaction_match']:10s} R={r['r_shifted']:7s} P={r['p_shifted']:7s}")
                else:
                    model_key = model_name
                    for cond in args.condition:
                        r = analyze_run(backend, pdb, model_key, model_dir, condition=cond)
                        all_results.append(r)
                        if not args.summary:
                            status_icon = {"PERFECT": "✓", "GOOD": "✓", "PARTIAL_REACT": "~",
                                           "WRONG_REACT": "✗", "ANOMALOUS": "!",
                                           "SVD_FAIL": "✗", "ERROR": "✗", "OPT_BROKEN": "✗"}.get(r["status"], "—")
                            print(f"{status_icon} {backend:14s} {pdb}/{model_key:20s} [{cond:7s}] "
                                  f"status={r['status']:16s} barrier={r['barrier_kcal']:>8s} "
                                  f"lit={r['lit_barrier']:>6s} err={r['barrier_error']:>7s} "
                                  f"irc={r['irc_reaction_match']:10s} R={r['r_shifted']:7s} P={r['p_shifted']:7s}")

    # --- Summary ---
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    for backend in args.backend:
        br = [r for r in all_results if r["backend"] == backend]
        if not br:
            continue
        n_total = len(br)
        n_perfect = sum(1 for r in br if r["status"] == "PERFECT")
        n_good = sum(1 for r in br if r["status"] in ("PERFECT", "GOOD"))
        n_partial = sum(1 for r in br if r["status"] == "PARTIAL_REACT")
        n_wrong = sum(1 for r in br if r["status"] == "WRONG_REACT")
        n_opt_broken = sum(1 for r in br if r["status"] == "OPT_BROKEN")
        n_fail = sum(1 for r in br if r["status"] in ("SVD_FAIL", "ERROR", "ANOMALOUS"))
        n_incomplete = sum(1 for r in br if r["status"] in ("RUNNING_OR_INCOMPLETE", "NOT_STARTED", "DIR_MISSING"))

        # MAE for runs with literature comparison
        errors = []
        for r in br:
            if r["barrier_error"] != "N/A":
                try:
                    errors.append(abs(float(r["barrier_error"])))
                except ValueError:
                    pass
        mae = f"{np.mean(errors):.1f}" if errors else "N/A"

        # Timing statistics
        times = []
        for r in br:
            if r["total_time_sec"] != "N/A":
                times.append(float(r["total_time_sec"]))
        time_str = ""
        if times:
            time_str = (f"  Time: avg={fmt_time(np.mean(times))} "
                        f"med={fmt_time(np.median(times))} "
                        f"min={fmt_time(min(times))} max={fmt_time(max(times))}")

        print(f"\n{backend}:")
        print(f"  Total={n_total}  Perfect={n_perfect}  Good={n_good}  "
              f"Partial={n_partial}  OptBroken={n_opt_broken}  Wrong={n_wrong}  "
              f"Fail={n_fail}  Incomplete={n_incomplete}  MAE={mae} kcal/mol")
        if time_str:
            print(time_str)

    # --- CSV output ---
    if args.csv:
        fields = list(all_results[0].keys()) if all_results else []
        with open(args.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(all_results)
        print(f"\nCSV written to {args.csv}")


if __name__ == "__main__":
    main()
