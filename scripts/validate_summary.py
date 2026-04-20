#!/usr/bin/env python3
"""
Compact validation summary: per-system O/△/X/- table.

For each (backend, pdb, model):
  1. Input R→P bond changes (from run.sh -i args) = expected reaction
  2. IRC-opt R→P bond changes (from structures/reactant.xyz → structures/product.xyz)
  3. Input R → structures/reactant.xyz = did R change during optimization?
  4. Input P → structures/product.xyz = did P change during optimization?

Judgment (non-metal bonds only):
  Pass (✓) = IRC-opt R→P contains ALL non-metal bond changes from input R→P
  Partial  = IRC-opt R→P contains SOME non-metal bond changes
  Fail (✗) = IRC-opt R→P contains NONE, or computation failed
  N/A  (—) = No IRC/structures data (incomplete)
"""

import csv
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from validate_benchmark import (
    run_bond_summary,
    find_input_endpoints,
    find_irc_endpoints,
    filter_nonmetal,
    judge_reaction,
    fmt_bc,
    bond_changes_empty,
    parse_timing,
    fmt_time,
    find_rp_inputs_for_tsonly,
    check_ts_imag_freq,
    CONDITIONS,
)

# Default to the shipped benchmark inputs under <repo>/examples/benchmark.
# Override with PDB2REACTION_BENCHMARK_DIR if your results live elsewhere.
BASE = Path(
    os.environ.get(
        "PDB2REACTION_BENCHMARK_DIR",
        str(Path(__file__).resolve().parent.parent / "examples" / "benchmark"),
    )
)
BACKENDS = ["uma_s1p1", "uma_s1p2", "uma_m1p1", "mace", "orb"]
CONDITION_OUT_FILES = {
    "default": "all.out",
    "nrp": "all_nrp.out",
    "tsonly": "all_tsonly.out",
    "ana": "all_ana.out",
}


def get_run_dir(backend, pdb, model):
    if "/" in model:
        parts = model.split("/")
        return BASE / f"p2r_cluster_{backend}" / pdb / parts[0] / parts[1]
    return BASE / f"p2r_cluster_{backend}" / pdb / model


def analyze(backend, pdb, model, condition="default"):
    """Analyze one (backend, pdb, model, condition). Returns dict with judge + details."""
    run_dir = get_run_dir(backend, pdb, model)
    if not run_dir.exists():
        return {"judge": "✗", "reason": "NO_DIR"}

    result_dir_name = CONDITIONS.get(condition, "result")
    out_file_name = CONDITION_OUT_FILES.get(condition, "all.out")

    # Check for error in .out
    all_out = run_dir / out_file_name
    if all_out.exists():
        tail = all_out.read_text()[-500:]
        if "SVD did not converge" in tail or "LinAlgError" in tail:
            return {"judge": "✗", "reason": "SVD_FAIL"}

    # Cluster systems are single-step reactions — multi-segment paths are failures
    if condition != "tsonly":
        import json as _json
        if condition == "tsonly":
            summary_json = run_dir / result_dir_name / "tsopt_single" / "summary.json"
        else:
            summary_json = run_dir / result_dir_name / "summary.json"
            if not summary_json.exists():
                summary_json = run_dir / result_dir_name / "path_opt" / "summary.json"
            if not summary_json.exists():
                summary_json = run_dir / result_dir_name / "path_search" / "summary.json"
        if summary_json.exists():
            try:
                with open(summary_json) as _f:
                    _data = _json.load(_f)
                _segs = _data.get("segments", []) if _data else []
                if len(_segs) > 1:
                    return {"judge": "✗", "reason": f"SPLIT_{len(_segs)}_SEGS"}
            except Exception:
                pass

    # Input endpoints
    if condition == "tsonly":
        inp_R, inp_P = find_rp_inputs_for_tsonly(run_dir)
    else:
        inp_R, inp_P = find_input_endpoints(run_dir, result_dir_name)
    if not inp_R or not inp_P:
        return {"judge": "✗", "reason": "NO_INPUT"}

    # Input R→P (non-metal)
    bc_input = run_bond_summary(inp_R, inp_P)
    bc_input_nm = filter_nonmetal(bc_input)

    # IRC-optimized structures
    irc_R, irc_P, has_irc = find_irc_endpoints(run_dir, result_dir_name)
    if not has_irc:
        if condition == "tsonly":
            summary = run_dir / result_dir_name / "tsopt_single" / "summary.json"
        else:
            summary = run_dir / result_dir_name / "summary.json"
        if not summary.exists():
            return {"judge": "—", "reason": "INCOMPLETE", "input_bc": fmt_bc(bc_input_nm)}
        return {"judge": "—", "reason": "NO_STRUCTURES", "input_bc": fmt_bc(bc_input_nm)}

    # Two-way assignment: try normal vs swapped and pick closer to input R/P
    from validate_benchmark import assign_irc_to_input
    try:
        irc_R2, irc_P2, _swapped = assign_irc_to_input(inp_R, inp_P, irc_R, irc_P)
        irc_R, irc_P = irc_R2, irc_P2
    except Exception:
        pass

    # IRC-opt R→P (non-metal)
    bc_irc = run_bond_summary(irc_R, irc_P)
    bc_irc_nm = filter_nonmetal(bc_irc)

    # Input R → opt R (non-metal)
    bc_r_shift = run_bond_summary(inp_R, irc_R)
    bc_r_shift_nm = filter_nonmetal(bc_r_shift)

    # Input P → opt P (non-metal)
    bc_p_shift = run_bond_summary(inp_P, irc_P)
    bc_p_shift_nm = filter_nonmetal(bc_p_shift)

    # Judge reaction reproduction
    judge, detail = judge_reaction(bc_input_nm, bc_irc_nm)

    r_ok = "OK" if (bc_r_shift_nm and bond_changes_empty(bc_r_shift_nm)) else "SHIFT"
    p_ok = "OK" if (bc_p_shift_nm and bond_changes_empty(bc_p_shift_nm)) else "SHIFT"

    if r_ok == "SHIFT" or p_ok == "SHIFT":
        inp_bonds = set(bc_input_nm["formed"]) | set(bc_input_nm["broken"])
        r_shift_bonds = set(bc_r_shift_nm["formed"] + bc_r_shift_nm["broken"]) if bc_r_shift_nm else set()
        p_shift_bonds = set(bc_p_shift_nm["formed"] + bc_p_shift_nm["broken"]) if bc_p_shift_nm else set()
        if r_shift_bonds & inp_bonds or p_shift_bonds & inp_bonds:
            judge = "✗"

    # TS imaginary mode count (informational; does not affect judge).
    n_imag = check_ts_imag_freq(run_dir, result_dir_name)

    return {
        "judge": judge,
        "input_bc": fmt_bc(bc_input_nm),
        "irc_bc": fmt_bc(bc_irc_nm),
        "r_shift": r_ok,
        "p_shift": p_ok,
        "r_shift_bc": fmt_bc(bc_r_shift_nm),
        "p_shift_bc": fmt_bc(bc_p_shift_nm),
        "ts_n_imag": n_imag if n_imag is not None else "N/A",
    }


def main():
    # Collect systems from CSV (or discover from filesystem)
    csv_path = BASE / "results_all.csv"
    if csv_path.exists():
        with open(csv_path) as f:
            rows = list(csv.DictReader(f))
        systems = []
        seen = set()
        for r in rows:
            key = (r["pdb"], r["model"])
            if key not in seen:
                seen.add(key)
                systems.append(key)
        lit_map = {}
        for r in rows:
            if r["lit_barrier"] != "N/A":
                lit_map[(r["pdb"], r["model"])] = r["lit_barrier"]
    else:
        print("Warning: results_all.csv not found. Run validate_benchmark.py first.", file=sys.stderr)
        return

    cond_list = list(CONDITIONS.keys())

    # === Table 1: O/△/X/- summary per condition ===
    for cond in cond_list:
        print(f"## Table 1 [{cond}]: Qualitative reaction reproduction (✓/✗/—)")
        print()
        print("| PDB  | Model/Step | 文献 | s1p1 | s1p2 | m1p1 | MACE | ORB |")
        print("|------|-----------|------|------|------|------|------|-----|")

        prev_pdb = ""
        all_v = {be: [] for be in BACKENDS}
        for pdb, model in systems:
            lit = lit_map.get((pdb, model), "")
            vals = []
            for be in BACKENDS:
                res = analyze(be, pdb, model, condition=cond)
                j = res["judge"]
                vals.append(j)
                all_v[be].append(j)
            sep = pdb if pdb != prev_pdb else ""
            prev_pdb = pdb
            m = model.replace("model_", "").replace("model/", "")
            print(f"| {sep:4s} | {m:9s} | {lit:>4s} | {vals[0]:>4s} | {vals[1]:>4s} | {vals[2]:>4s} | {vals[3]:>4s} | {vals[4]:>3s} |")

        print("|------|-----------|------|------|------|------|------|-----|")
        for be in BACKENDS:
            v = all_v[be]
            o, t, x, d = v.count("✓"), v.count("~"), v.count("✗"), v.count("—")
            print(f"| {be:>21s} || ✓={o} ~={t} ✗={x} —={d} / {len(v)} |")
        print()

    # === Table 2: Detailed bond changes (default condition only) ===
    print()
    print("## Table 2 [default]: Detailed bond changes (non-metal only)")
    print()
    all_results = {}
    for pdb, model in systems:
        for be in BACKENDS:
            all_results[(be, pdb, model)] = analyze(be, pdb, model, condition="default")
    for pdb, model in systems:
        m = model.replace("model_", "").replace("model/", "")
        print(f"### {pdb}/{m}")
        for be in BACKENDS:
            res = all_results.get((be, pdb, model), {})
            if "input_bc" in res:
                print(f"  Input R→P: {res['input_bc']}")
                break

        for be in BACKENDS:
            res = all_results.get((be, pdb, model), {})
            j = res["judge"]
            irc_bc = res.get("irc_bc", "N/A")
            r_shift = res.get("r_shift", "N/A")
            p_shift = res.get("p_shift", "N/A")
            r_shift_bc = res.get("r_shift_bc", "")
            p_shift_bc = res.get("p_shift_bc", "")
            reason = res.get("reason", "")
            n_imag = res.get("ts_n_imag", "N/A")

            detail = f"IRC R→P: {irc_bc}" if irc_bc != "N/A" else reason
            opt_info = ""
            if r_shift != "N/A":
                opt_info = f" | optR={r_shift}"
                if r_shift == "SHIFT":
                    opt_info += f"({r_shift_bc})"
                opt_info += f" optP={p_shift}"
                if p_shift == "SHIFT":
                    opt_info += f"({p_shift_bc})"
            nimag_info = f" | nImag={n_imag}" if n_imag != "N/A" else ""

            print(f"  {be:12s} [{j}] {detail}{opt_info}{nimag_info}")
        print()

    # === Table 3: Timing ===
    print()
    print_timing_table(systems)


def print_timing_table(systems):
    """Print Table 3: Pipeline wall-clock time per system per backend per condition."""
    cond_list = list(CONDITIONS.keys())

    for cond in cond_list:
        out_file = CONDITION_OUT_FILES[cond]
        print(f"## Table 3 [{cond}]: Pipeline wall-clock time")
        print()
        print("| PDB  | Model/Step | s1p1 | s1p2 | m1p1 | MACE | ORB |")
        print("|------|-----------|------|------|------|------|-----|")

        prev_pdb = ""
        backend_times = {be: [] for be in BACKENDS}

        for pdb, model in systems:
            vals = []
            for be in BACKENDS:
                run_dir = get_run_dir(be, pdb, model)
                all_out = run_dir / out_file
                if all_out.exists():
                    timing = parse_timing(str(all_out))
                    total = timing.get("total_sec")
                    vals.append(fmt_time(total))
                    if total is not None:
                        backend_times[be].append(total)
                else:
                    vals.append("—")
            sep = pdb if pdb != prev_pdb else ""
            prev_pdb = pdb
            m = model.replace("model_", "").replace("model/", "")
            print(f"| {sep:4s} | {m:9s} | {vals[0]:>6s} | {vals[1]:>6s} | {vals[2]:>6s} | {vals[3]:>6s} | {vals[4]:>5s} |")

        print("|------|-----------|------|------|------|------|-----|")
        avgs = []
        for be in BACKENDS:
            t = backend_times[be]
            avgs.append(fmt_time(np.mean(t)) if t else "—")
        print(f"| {'Avg':>15s} || {avgs[0]:>6s} | {avgs[1]:>6s} | {avgs[2]:>6s} | {avgs[3]:>6s} | {avgs[4]:>5s} |")
        meds = []
        for be in BACKENDS:
            t = backend_times[be]
            meds.append(fmt_time(np.median(t)) if t else "—")
        print(f"| {'Median':>15s} || {meds[0]:>6s} | {meds[1]:>6s} | {meds[2]:>6s} | {meds[3]:>6s} | {meds[4]:>5s} |")
        print()

        # Stage breakdown (DFT separated)
        print(f"### Stage breakdown (average seconds) [{cond}]")
        print()
        print("| Backend | Path Opt | TS Opt | IRC | Freq | DFT | Total |")
        print("|---------|----------|--------|-----|------|-----|-------|")
        for be in BACKENDS:
            stage_totals = {"path_opt_sec": [], "ts_opt_sec": [], "irc_sec": [],
                            "freq_total_sec": [], "dft_total_sec": [], "total_sec": []}
            for pdb, model in systems:
                run_dir = get_run_dir(be, pdb, model)
                all_out = run_dir / out_file
                if all_out.exists():
                    timing = parse_timing(str(all_out))
                    for key in stage_totals:
                        if timing[key] is not None:
                            stage_totals[key].append(timing[key])
            row = [be]
            for key in ("path_opt_sec", "ts_opt_sec", "irc_sec", "freq_total_sec", "dft_total_sec", "total_sec"):
                vals = stage_totals[key]
                row.append(fmt_time(np.mean(vals)) if vals else "—")
            print(f"| {row[0]:7s} | {row[1]:>8s} | {row[2]:>6s} | {row[3]:>5s} | {row[4]:>6s} | {row[5]:>5s} | {row[6]:>7s} |")
        print()


if __name__ == "__main__":
    main()
