#!/usr/bin/env bash
# pdb2reaction smoke tests — GPU required.
# Coverage: all subcommands, input formats, --ref-pdb, solvent, dry-run, utilities.
#
# This script assumes the calling environment already has:
#   - a Python with `pdb2reaction` installed and importable
#   - CUDA available
#   - a writable scratch copy of this directory (artefacts land in `testNN*`)
# It does NOT activate conda, load modules, or contain HPC scheduler
# directives. Copy the fixtures to scratch, then run from that copy:
#
#   cp -a tests/smoke /path/to/scratch/p2r-smoke
#   cd /path/to/scratch/p2r-smoke
#   bash run.sh
#
# If you need an HPC scheduler wrapper, keep that in
# your own out-of-tree submission script and have it invoke this file as
# the body.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel 2>/dev/null || true)"
if [[ -n "$REPO_ROOT" ]]; then
  echo "ERROR: copy tests/smoke to writable scratch and run that copy; do not run the repository script in place." >&2
  exit 2
fi

# Deterministic run gate: pin PYTHONHASHSEED so Set / dict iteration order
# does not leak hash randomisation into the produced output. Combined with the
# UMA backend's deterministic-algorithms wiring this removes the remaining
# I/O-side non-determinism between two consecutive runs.
export PYTHONHASHSEED=0
# Reduce CUDA allocator fragmentation across the 40+ stage processes.
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"

# Every capability in this lane is required. Missing prerequisites block the
# release lane before case 1; they are never converted into a passing skip.
command -v xtb >/dev/null 2>&1 || {
  echo "[smoke] BLOCKED: xtb is required by the solvent cell" >&2
  exit 2
}
python - <<'PY'
from importlib.metadata import version
from pathlib import Path
import subprocess
import sys

import pdb2reaction
import torch

installed = version("pdb2reaction")
if pdb2reaction.__version__ != installed:
    raise SystemExit(
        "[smoke] BLOCKED: distribution/module version mismatch: "
        f"metadata={installed}, module={pdb2reaction.__version__}"
    )
cli_version = subprocess.run(
    [sys.executable, "-m", "pdb2reaction", "--version"],
    check=True,
    capture_output=True,
    text=True,
).stdout.strip().split()[-1]
if cli_version != installed:
    raise SystemExit(
        "[smoke] BLOCKED: distribution/CLI version mismatch: "
        f"metadata={installed}, cli={cli_version}"
    )
if not torch.cuda.is_available():
    raise SystemExit("[smoke] BLOCKED: CUDA is required by this lane")
print(
    "[smoke] package pdb2reaction "
    f"version={installed} source={Path(pdb2reaction.__file__).resolve()}"
)
PY
pdb2reaction() { python -m pdb2reaction "$@"; }
python assert_tr_cuda_parity.py

# Clean only artifacts authored by this harness. The digit-qualified glob must
# not be widened to `test*`, which would also match the repository's tests/.
rm -rf -- test[0-9]* p_complex_model.pdb model_r.pdb model_multi.pdb model_from_cif.pdb model_from_cif.cif r_complex_elem.pdb r_complex_fixalt.pdb

P_COMPLEX_MODEL_FREEZE_ATOMS="1,32"

# --- Subcommand tests (individual) ---

# test1: opt (grad / lbfgs)
pdb2reaction opt -i r.pdb -q -1 --opt-mode grad --max-cycles 5 --thresh gau_loose --dump --out-dir test1 > test1.out 2>&1

# test2: opt (hess / rfo)
pdb2reaction opt -i r.pdb -q -1 --opt-mode hess --max-cycles 3 --thresh gau_loose --out-dir test2 > test2.out 2>&1

# test3: tsopt (hess)
pdb2reaction tsopt -i ts.pdb -q 0 --max-cycles 5 --thresh gau_loose --out-dir test3 > test3.out 2>&1

# test4: freq
pdb2reaction freq -i r.pdb -q -1 --out-dir test4 > test4.out 2>&1

# test5: irc
pdb2reaction irc -i ts.pdb -q 0 --max-cycles 3 --never-stop --out-dir test5 > test5.out 2>&1

# test6: dft (lightweight: hf/sto-3g, cpu)
pdb2reaction dft -i h2.gjf --func-basis 'hf/sto-3g' --grid-level 0 --conv-tol 1e-5 --max-cycle 40 --engine cpu --out-dir test6 > test6.out 2>&1

# test7: scan (1D)
pdb2reaction scan -i r.pdb -q -1 --scan-lists "[(1,5,1.4)]" --max-step-size 2.0 --relax-max-cycles 3 --no-preopt --no-endopt --out-dir test7 > test7.out 2>&1

# test8: scan2d (extract model first)
pdb2reaction extract -i p_complex.pdb -c 'PRE' --ligand-charge 'PRE:-2' -r 5.0 --no-exclude-backbone --include-h2o -o p_complex_model.pdb
pdb2reaction scan2d -i p_complex_model.pdb --ligand-charge 'PRE:-2' --freeze-atoms "$P_COMPLEX_MODEL_FREEZE_ATOMS" --scan-lists "[('PRE 8 C3','PRE 8 O1\'',1.4,1.8),('PRE 8 C1','PRE 8 C8',3.0,3.4)]" --max-step-size 2.0 --relax-max-cycles 100 --thresh gau_loose --out-dir test8 > test8.out 2>&1

# test9: scan3d
pdb2reaction scan3d -i p_complex_model.pdb --ligand-charge 'PRE:-2' --freeze-atoms "$P_COMPLEX_MODEL_FREEZE_ATOMS" --scan-lists "[('PRE 8 C3','PRE 8 O1\'',1.4,1.8),('PRE 8 C1','PRE 8 C8',3.0,3.4),('PRE 8 C1','PRE 8 C7',1.4,1.6)]" --max-step-size 2.0 --relax-max-cycles 100 --thresh gau_loose --out-dir test9 > test9.out 2>&1

# test10: path-opt (gsm)
pdb2reaction path-opt -i r.pdb p.pdb -q -1 --max-nodes 5 --max-cycles 5 --no-preopt --no-climb --out-dir test10 > test10.out 2>&1

# test11: path-search
pdb2reaction path-search -i r.pdb p.pdb -q -1 --max-nodes 5 --max-cycles 5 --out-dir test11 > test11.out 2>&1

# --- Input format tests (all command) ---

# test12: all (pdb+pdb, --no-refine-path for single-pass path-opt)
pdb2reaction -i r.pdb p.pdb -q -1 --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --out-dir test12 > test12.out 2>&1

# test13: all (xyz+xyz, --ref-pdb for PDB conversion, --no-refine-path for single-pass path-opt)
pdb2reaction -i r.xyz p.xyz -q -1 --ref-pdb r.pdb --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --out-dir test13 > test13.out 2>&1

# test14: all (gjf+gjf, --no-refine-path for single-pass path-opt)
pdb2reaction -i r.gjf p.gjf --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --out-dir test14 > test14.out 2>&1

# test15: all (scan-lists, pdb)
pdb2reaction -i r.pdb -q -1 --scan-lists "[(1,5,1.4)]" --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --out-dir test15 > test15.out 2>&1

# test16: all (scan-lists, xyz)
pdb2reaction -i r.xyz -q -1 --scan-lists "[(1,5,1.4)]" --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --out-dir test16 > test16.out 2>&1

# --- Complex system tests ---

# test17: opt (complex, ligand-charge)
pdb2reaction opt -i p_complex.pdb --ligand-charge 'PRE:-2' --max-cycles 3 --thresh gau_loose --out-dir test17 > test17.out 2>&1

# test18: all (complex, extract + scan-lists)
pdb2reaction -i p_complex.pdb -c 'PRE' --ligand-charge 'PRE:-2' --scan-lists "[('PRE 8 C3','PRE 8 O1\'',1.4),('PRE 8 C1','PRE 8 C8',3.3)]" -r 5.0 --no-exclude-backbone --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --out-dir test18 > test18.out 2>&1

# test19: all (complex, multi-input, --no-refine-path for single-pass path-opt)
pdb2reaction -i r_complex.pdb p_complex.pdb -c 'PRE' --ligand-charge 'PRE:-2' -r 5.0 --no-exclude-backbone --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --out-dir test19 > test19.out 2>&1

# --- TSOPT-only mode ---

# test20: all (ts input, --tsopt + --thermo)
pdb2reaction -i ts.pdb -q 0 --tsopt --opt-mode-post grad --thermo --irc-never-stop --max-cycles 100 --thresh gau --thresh-post gau_loose --out-dir test20 > test20.out 2>&1
python assert_release_result.py all test20 --require-thermo >> test20.out 2>&1

# test21: all (ts input, --tsopt, opt-mode hess)
pdb2reaction -i ts.pdb -q 0 --tsopt --opt-mode-post hess --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --out-dir test21 > test21.out 2>&1

# --- MEP mode ---

# test22: all (pdb+pdb, mep-mode dmf, --no-refine-path for single-pass path-opt)
pdb2reaction -i r.pdb p.pdb -q -1 --mep-mode dmf --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --out-dir test22 > test22.out 2>&1

# --- TSOPT for complex systems ---

# test23: tsopt (complex, hess)
pdb2reaction tsopt -i ts_complex.pdb -q -1 --opt-mode hess --max-cycles 3 --thresh gau_loose --out-dir test23 > test23.out 2>&1

# test24: tsopt (complex, grad)
pdb2reaction tsopt -i ts_complex.pdb -q -1 --opt-mode grad --max-cycles 100 --thresh gau --out-dir test24 > test24.out 2>&1

# --- Dry-run validation ---

# test25: opt --dry-run
pdb2reaction opt -i r.pdb -q -1 --dry-run --out-dir test25 > test25.out 2>&1

# test26: tsopt --dry-run
pdb2reaction tsopt -i ts.pdb -q 0 --dry-run --out-dir test26 > test26.out 2>&1

# test27: freq --dry-run
pdb2reaction freq -i r.pdb -q -1 --dry-run --out-dir test27 > test27.out 2>&1

# test28: dft --dry-run
pdb2reaction dft -i r.pdb -q -1 --dry-run --out-dir test28 > test28.out 2>&1

# test29: path-search --dry-run
pdb2reaction path-search -i r.pdb p.pdb -q -1 --dry-run --out-dir test29 > test29.out 2>&1

# test30: irc --dry-run
pdb2reaction irc -i ts.pdb -q 0 --dry-run --out-dir test30 > test30.out 2>&1

# --- Utility subcommands ---

# test31: extract
pdb2reaction extract -i r_complex.pdb -c 'PRE' --ligand-charge 'PRE:-2' -r 5.0 --no-exclude-backbone -o model_r.pdb > test31.out 2>&1

# test32: add-elem-info
pdb2reaction add-elem-info -i r_complex.pdb -o r_complex_elem.pdb > test32.out 2>&1

# test33: trj2fig
pdb2reaction trj2fig -i test1/optimization_trj.xyz -o test33.png > test33.out 2>&1

# test34: energy-diagram
pdb2reaction energy-diagram -i "[0, 12.5, 4.3, 18.7, -1.2]" -o test34.png > test34.out 2>&1

# --- Bond-summary & fix-altloc ---

# test35: bond-summary (two PDB structures)
pdb2reaction bond-summary -i r.pdb p.pdb > test35.out 2>&1

# test36: fix-altloc
pdb2reaction fix-altloc -i r_complex.pdb -o r_complex_fixalt.pdb > test36.out 2>&1

# --- YAML scan spec ---

# test37: scan (1D, YAML spec file)
pdb2reaction scan -i r.pdb -q -1 -s scan_spec.yaml --max-step-size 2.0 --relax-max-cycles 3 --no-preopt --no-endopt --out-dir test37 > test37.out 2>&1

# test38: scan2d (YAML spec file)
pdb2reaction scan2d -i p_complex_model.pdb --ligand-charge 'PRE:-2' --freeze-atoms "$P_COMPLEX_MODEL_FREEZE_ATOMS" -s scan2d_spec.yaml --max-step-size 2.0 --relax-max-cycles 100 --thresh gau_loose --out-dir test38 > test38.out 2>&1

# --- Solvent correction (requires xTB) ---

# test39: opt (solvent water; xTB is a lane preflight requirement)
pdb2reaction opt -i r.pdb -q -1 --opt-mode grad --max-cycles 3 --thresh gau_loose --solvent water --out-dir test39 > test39.out 2>&1

# --- dist-freeze ---

# test40: opt --dist-freeze --dry-run (inline 3-tuple + 2-tuple)
pdb2reaction opt -i r.pdb -q -1 --dist-freeze "[(1,2,1.5),(3,4)]" --dry-run --out-dir test40 > test40.out 2>&1

# --- refine-path ---

# test41: all (pdb+pdb, --refine-path = recursive path_search opt-in; the default is now single-pass path-opt)
pdb2reaction -i r.pdb p.pdb -q -1 --refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --out-dir test41 > test41.out 2>&1

# --- Polish-train new CLI flags (A1 + W3 + B4 wires; all opt-in, defaults preserve Table 1 numerics) ---

# test42: tsopt --opt-mode trim (A1 Helgaker trust-region image-min TS opt)
pdb2reaction tsopt -i ts.pdb -q 0 --opt-mode trim --max-cycles 5 --thresh gau_loose --out-dir test42 > test42.out 2>&1

# test43: tsopt --opt-mode rsprfo (A1 Banerjee restricted-step P-RFO TS opt)
pdb2reaction tsopt -i ts.pdb -q 0 --opt-mode rsprfo --max-cycles 5 --thresh gau_loose --out-dir test43 > test43.out 2>&1

# test45: irc --irc-pos-def (Sella backport 1: PSD-Hessian convergence guard)
pdb2reaction irc -i ts.pdb -q 0 --max-cycles 3 --irc-pos-def --out-dir test45 > test45.out 2>&1

# test46: opt --print-every 3 (W3a debug throttle, no behavior change)
pdb2reaction opt -i r.pdb -q -1 --opt-mode hess --max-cycles 5 --thresh gau_loose --print-every 3 --out-dir test46 > test46.out 2>&1

# --- Determinism gate ---

# test47: `all` pipeline determinism GATE (`--deterministic`).
# Runs the full pipeline twice with identical inputs / args + `--deterministic`
# and REQUIRES the two runs to be bit-identical. Default (non-deterministic) GPU
# runs carry ~ULP scatter/atomic non-determinism and are not asserted here;
# `--deterministic` enables torch deterministic algorithms and MUST be
# bit-reproducible, so any drift is a real regression and fails the smoke.
det_args="-i r.pdb p.pdb -q -1 --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --deterministic"
pdb2reaction $det_args --out-dir test47_a > test47_a.out 2>&1
pdb2reaction $det_args --out-dir test47_b > test47_b.out 2>&1
find test47_a -type f \( -name "*.pdb" -o -name "*.xyz" \) -printf '%P\n' | LC_ALL=C sort > test47_a.manifest
find test47_b -type f \( -name "*.pdb" -o -name "*.xyz" \) -printf '%P\n' | LC_ALL=C sort > test47_b.manifest
if ! cmp -s test47_a.manifest test47_b.manifest; then
  echo "[smoke] FAIL test47: deterministic runs produced different file manifests" > test47.out
  comm -3 test47_a.manifest test47_b.manifest >> test47.out
  cat test47.out
  exit 1
fi
total=$(wc -l < test47_a.manifest)
if [ "$total" -eq 0 ]; then
  echo "[smoke] FAIL test47: deterministic gate found no PDB/XYZ artifacts" > test47.out
  cat test47.out
  exit 1
fi
drifted=0
while IFS= read -r rel; do
  if ! cmp -s "test47_a/$rel" "test47_b/$rel"; then
    drifted=$((drifted + 1))
    echo "DRIFT: $rel" >> test47.out
  fi
done < test47_a.manifest
echo "[det_check] compared $total PDB/XYZ files; $drifted differ" >> test47.out
if [ "$drifted" -ne 0 ]; then
  echo "[smoke] FAIL test47: --deterministic runs differ" >> test47.out
  cat test47.out
  exit 1
fi

# --- --coord-type CLI plumbing (throttled, fast) ---

# test48: `all --coord-type cart` — explicit cart (== default), verifies CLI plumbing.
pdb2reaction all -i r.pdb p.pdb -q -1 --coord-type cart --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --out-dir test48 > test48.out 2>&1

# test49: `all --coord-type dlc` — DLC propagated to child opt / tsopt / path-opt stages.
pdb2reaction all -i r.pdb p.pdb -q -1 --coord-type dlc --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --out-dir test49 > test49.out 2>&1

# test50: `sp` (single-point) — energy + forces.
pdb2reaction sp -i r.pdb -q -1 --out-dir test50 > test50.out 2>&1

# test51: `sp --hess` — energy + forces + Hessian (UMA analytical).
pdb2reaction sp -i r.pdb -q -1 --hess --out-dir test51 > test51.out 2>&1

# --- Full-pipeline, NO max-cycles throttle (release-gate runs) ---
# These exercise the canonical `all` flow with default convergence thresholds
# and the production-realistic optimizer cycle budgets. Each run takes
# substantially longer (~30-60 min) than the throttled tests above; they are
# the "does the pipeline actually finish on a real input" gate.

# test52: full `all` cart — default thresh, no max-cycles cap.
pdb2reaction all -i r.pdb p.pdb -q -1 --out-dir test52 > test52.out 2>&1

# test53: `all` dlc — verifies the DLC code-path lights up end-to-end.
# Capped at max-cycles 5 + thresh gau_loose + --no-tsopt/thermo/dft because
# DLC GSM on this system needs many cycles to converge with default `gau`
# thresh, and the post-stages (TS / IRC / DFT) depend on a converged HEI
# from the MEP — silently broken structure handoff otherwise. test52 (cart)
# keeps the no-cap default-behaviour check.
pdb2reaction all -i r.pdb p.pdb -q -1 --coord-type dlc --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --out-dir test53 > test53.out 2>&1

# --- Per-stage internal-coordinate code paths ---
# Each test scoped at max-cycles 3 + thresh gau_loose to exercise the
# --coord-type code path without long convergence. Frequency analysis remains
# Cartesian because the only thing an internal-coordinate flag would add is
# an unrelated B-matrix construction:
#   - freq: vibrational analysis uses cart Hessian regardless of
#     --coord-type; the only thing the flag adds is an SVD on the internal
#     B-matrix that occasionally fails on cuSOLVER. freq + cart is already
#     covered by test4.

# test53a: opt --coord-type dlc
pdb2reaction opt -i r.pdb -q -1 --coord-type dlc --max-cycles 3 --thresh gau_loose --out-dir test53a_opt_dlc > test53a_opt_dlc.out 2>&1

# test53b-e: Hessian TS optimization must accept every documented coordinate
# system. The DLC case includes frozen atoms and therefore exercises the
# partial-Cartesian-Hessian -> internal-coordinate handoff regression.
pdb2reaction tsopt -i ts.pdb -q 0 --opt-mode hess --coord-type cart --max-cycles 2 --thresh gau_loose --out-dir test53b_ts_cart > test53b_ts_cart.out 2>&1
pdb2reaction tsopt -i ts.pdb -q 0 --opt-mode hess --coord-type redund --max-cycles 2 --thresh gau_loose --out-dir test53c_ts_redund > test53c_ts_redund.out 2>&1
pdb2reaction tsopt -i ts.pdb -q 0 --opt-mode hess --coord-type dlc --freeze-atoms 1,3,5 --max-cycles 2 --thresh gau_loose --out-dir test53d_ts_dlc_freeze > test53d_ts_dlc_freeze.out 2>&1
pdb2reaction tsopt -i ts.pdb -q 0 --opt-mode hess --coord-type tric --max-cycles 2 --thresh gau_loose --out-dir test53e_ts_tric > test53e_ts_tric.out 2>&1

# test53d: scan --coord-type dlc (1D)
pdb2reaction scan -i r.pdb -q -1 --coord-type dlc --scan-lists "[(1,5,1.4)]" --max-step-size 2.0 --relax-max-cycles 3 --no-preopt --no-endopt --out-dir test53d_scan_dlc > test53d_scan_dlc.out 2>&1

# test53g: opt --coord-type redund
pdb2reaction opt -i r.pdb -q -1 --coord-type redund --max-cycles 3 --thresh gau_loose --out-dir test53g_opt_redund > test53g_opt_redund.out 2>&1

# test53j: scan --coord-type redund
pdb2reaction scan -i r.pdb -q -1 --coord-type redund --scan-lists "[(1,5,1.4)]" --max-step-size 2.0 --relax-max-cycles 3 --no-preopt --no-endopt --out-dir test53j_scan_redund > test53j_scan_redund.out 2>&1

# test53k: opt --coord-type tric
pdb2reaction opt -i r.pdb -q -1 --coord-type tric --max-cycles 3 --thresh gau_loose --out-dir test53k_opt_tric > test53k_opt_tric.out 2>&1

# --- Multi-mode flag code-path verify (single-stage) ---

# test53m: opt --precision fp64 (alternate precision)
pdb2reaction opt -i r.pdb -q -1 --precision fp64 --max-cycles 3 --thresh gau_loose --out-dir test53m_opt_fp64 > test53m_opt_fp64.out 2>&1

# test53n: opt --precision fp32 (explicit fp32 dispatch; default is fp32, this pins the explicit path alongside test53m fp64)
pdb2reaction opt -i r.pdb -q -1 --precision fp32 --max-cycles 3 --thresh gau_loose --out-dir test53n_opt_fp32 > test53n_opt_fp32.out 2>&1

# test54: full `all` with `--backend orb` — exercises the non-default MLIP backend.
pdb2reaction all -i r.pdb p.pdb -q -1 --backend orb --out-dir test54 > test54.out 2>&1
python assert_release_result.py provenance test54 --expected-backend orb --expected-model orb_v3_conservative_omol --expected-precision fp64 >> test54.out 2>&1

# ---- Coverage-gap regression (subcommand-specific code paths; coverage audit 2026-06-05) ----
# test55: opt --dist-freeze + --bias-k (harmonic restraint actually applied at runtime, not dry-run)
pdb2reaction opt -i r.pdb -q -1 --dist-freeze "[(1,2,1.5)]" --bias-k 150 --max-cycles 3 --thresh gau_loose --out-dir test55_opt_biask > test55_opt_biask.out 2>&1

# test56: opt --freeze-atoms (explicit user-frozen DOF, distinct from --freeze-links)
pdb2reaction opt -i r.pdb -q -1 --freeze-atoms '1,3,5' --max-cycles 3 --thresh gau_loose --no-flatten --out-dir test56_opt_freeze > test56_opt_freeze.out 2>&1

# test57: freq --hessian-calc-mode Analytical (workflow analytical path)
pdb2reaction freq -i r.pdb -q -1 --hessian-calc-mode Analytical --max-write 3 --out-dir test57_freq_anahess > test57_freq_anahess.out 2>&1

# test58: irc --hessian-calc-mode analytical (IRC-initiating Hessian path)
pdb2reaction irc -i ts.pdb -q 0 --hessian-calc-mode analytical --max-cycles 2 --out-dir test58_irc_anahess > test58_irc_anahess.out 2>&1

# test59: scan --opt-mode hess (RFO relaxation vs default LBFGS)
pdb2reaction scan -i r.pdb -q -1 --opt-mode hess --scan-lists "[(1,5,1.4)]" --max-step-size 2.0 --relax-max-cycles 2 --no-preopt --no-endopt --out-dir test59_scan_hess > test59_scan_hess.out 2>&1

# test60: scan2d --opt-mode hess (RFO per-grid relaxation)
pdb2reaction scan2d -i p_complex_model.pdb --ligand-charge 'PRE:-2' --freeze-atoms "$P_COMPLEX_MODEL_FREEZE_ATOMS" --scan-lists "[('PRE 8 C3','PRE 8 O1\'',1.4,1.8),('PRE 8 C1','PRE 8 C8',3.0,3.4)]" --opt-mode hess --max-step-size 2.0 --relax-max-cycles 3 --thresh gau_loose --out-dir test60_scan2d_hess > test60_scan2d_hess.out 2>&1

# test61: scan3d --opt-mode hess (RFO per-grid relaxation)
pdb2reaction scan3d -i p_complex_model.pdb --ligand-charge 'PRE:-2' --freeze-atoms "$P_COMPLEX_MODEL_FREEZE_ATOMS" --scan-lists "[('PRE 8 C3','PRE 8 O1\'',1.4,1.8),('PRE 8 C1','PRE 8 C8',3.0,3.4),('PRE 8 C1','PRE 8 C7',1.4,1.6)]" --opt-mode hess --max-step-size 2.0 --relax-max-cycles 3 --thresh gau_loose --out-dir test61_scan3d_hess > test61_scan3d_hess.out 2>&1

# test62: scan3d --precision fp64 (fp64 backend dispatch)
pdb2reaction scan3d -i p_complex_model.pdb --ligand-charge 'PRE:-2' --freeze-atoms "$P_COMPLEX_MODEL_FREEZE_ATOMS" --scan-lists "[('PRE 8 C3','PRE 8 O1\'',1.4,1.8),('PRE 8 C1','PRE 8 C8',3.0,3.4),('PRE 8 C1','PRE 8 C7',1.4,1.6)]" --precision fp64 --max-step-size 2.0 --relax-max-cycles 3 --thresh gau_loose --out-dir test62_scan3d_fp64 > test62_scan3d_fp64.out 2>&1

# test63: path-opt --mep-mode dmf (Direct Max Flux at the subcommand level)
pdb2reaction path-opt -i r.pdb p.pdb -q -1 --mep-mode dmf --max-nodes 5 --max-cycles 3 --no-preopt --no-climb --out-dir test63_pathopt_dmf > test63_pathopt_dmf.out 2>&1

# test64: path-opt --coord-type dlc (p2r keeps DLC for pure-MLIP)
pdb2reaction path-opt -i r.pdb p.pdb -q -1 --coord-type dlc --max-nodes 5 --max-cycles 3 --no-preopt --no-climb --out-dir test64_pathopt_dlc > test64_pathopt_dlc.out 2>&1

# test65: path-search --mep-mode dmf
pdb2reaction path-search -i r.pdb p.pdb -q -1 --mep-mode dmf --max-nodes 5 --max-cycles 3 --no-preopt --no-endopt --out-dir test65_ps_dmf > test65_ps_dmf.out 2>&1

# test66: path-search --opt-mode hess (RFO single-structure preopt; keep preopt ON)
pdb2reaction path-search -i r.pdb p.pdb -q -1 --opt-mode hess --workers 1 --max-nodes 5 --max-cycles 3 --no-endopt --out-dir test66_ps_hess > test66_ps_hess.out 2>&1

# test67: required positive MEP -> TSopt -> IRC -> thermo handoff.
# --tsopt-max-cycles must cover the opt-in --flatten repair (which draws from this
# global budget); this system converges to n_imag=1 with no flatten iterations, but
# 200 keeps the budget above flatten_max_iter (=50, ~135 cycles) so a soft spectator
# mode could never be starved short of a clean saddle. Product default is 10000.
pdb2reaction all -i r.pdb p.pdb -q -1 --tsopt --thermo --flatten --irc-never-stop --max-cycles 5 --tsopt-max-cycles 200 --thresh gau_loose --thresh-post gau --out-dir test67_all_tsopt > test67_all_tsopt.out 2>&1
python assert_release_result.py all test67_all_tsopt --require-thermo >> test67_all_tsopt.out 2>&1
grep -Fq '[irc] Reusing cached TS Hessian from tsopt.' test67_all_tsopt.out || { echo '[smoke] FAIL test67: IRC did not report cached TS Hessian reuse' >> test67_all_tsopt.out; exit 1; }

# test68: all --scan-lists (single-PDB scan->path mode of `all`)
pdb2reaction all -i r.pdb -q -1 --scan-lists "[(1,5,1.4)]" --max-cycles 5 --thresh gau_loose --no-tsopt --no-thermo --no-dft --out-dir test68_all_scan > test68_all_scan.out 2>&1

# test69: opt --coord-type dlc on extracted complex model with explicit frozen boundary atoms
pdb2reaction opt -i p_complex_model.pdb --ligand-charge 'PRE:-2' --coord-type dlc --freeze-atoms "$P_COMPLEX_MODEL_FREEZE_ATOMS" --max-cycles 3 --thresh gau_loose --out-dir test69_opt_complex_freeze_dlc > test69_opt_complex_freeze_dlc.out 2>&1
python check_frozen_atoms.py p_complex_model.pdb test69_opt_complex_freeze_dlc/final_geometry.pdb "$P_COMPLEX_MODEL_FREEZE_ATOMS" test69 >> test69_opt_complex_freeze_dlc.out 2>&1

# test70: scan --coord-type dlc on extracted complex model with explicit frozen boundary atoms.
# Keep preopt enabled: this is the regression path that can trigger internal-coordinate rebuilds.
# Use a small non-reactive target so this tests coordinate integrity, not chemistry.
pdb2reaction scan -i p_complex_model.pdb --ligand-charge 'PRE:-2' --coord-type dlc --freeze-atoms "$P_COMPLEX_MODEL_FREEZE_ATOMS" --scan-lists "[('PRE 8 C3','PRE 8 O1\'',3.2)]" --max-step-size 0.2 --relax-max-cycles 3 --no-endopt --out-dir test70_scan_complex_freeze_dlc > test70_scan_complex_freeze_dlc.out 2>&1
python check_frozen_atoms.py p_complex_model.pdb test70_scan_complex_freeze_dlc/stage_01/result.pdb "$P_COMPLEX_MODEL_FREEZE_ATOMS" test70 >> test70_scan_complex_freeze_dlc.out 2>&1
if grep -q "Covalent-bond changes (start vs final): Yes" test70_scan_complex_freeze_dlc.out; then
  echo "[bond-check] test70: unexpected covalent-bond changes in non-reactive DLC+freeze scan" >> test70_scan_complex_freeze_dlc.out
  exit 1
fi

# --- refine-path opt-in (recursive path_search) extra coverage ---
# Since the `all` default is now single-pass path-opt, exercise the recursive
# path_search opt-in (`--refine-path`) across more input modes than test41.

# test71: all --scan-lists --refine-path (single-PDB scan -> recursive path_search)
pdb2reaction all -i r.pdb -q -1 --scan-lists "[(1,5,1.4)]" --refine-path --max-cycles 5 --thresh gau_loose --no-tsopt --no-thermo --no-dft --out-dir test71_rp_scan > test71_rp_scan.out 2>&1

# test72: all (complex multi-input) --refine-path (recursive path_search on the extracted model)
pdb2reaction all -i r_complex.pdb p_complex.pdb -c 'PRE' --ligand-charge 'PRE:-2' -r 5.0 --no-exclude-backbone --refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --out-dir test72_rp_complex > test72_rp_complex.out 2>&1

# test73: extract MULTI-INPUT via space-separated '-i a.pdb b.pdb' (one flag, two paths).
# Regression guard: a single -i with several space-separated paths must NOT drop the 2nd input.
# A single -o yields one multi-MODEL PDB, so both endpoints must appear (-> exactly 2 MODEL records).
pdb2reaction extract -i r_complex.pdb p_complex.pdb -c 'PRE' --ligand-charge 'PRE:-2' -r 5.0 --no-exclude-backbone -o model_multi.pdb > test73_multi_extract.out 2>&1
n_models=$(grep -c '^MODEL' model_multi.pdb 2>/dev/null)
if [ "${n_models:-0}" -ne 2 ]; then
  echo "[extract-multi] test73: space-separated '-i a b' yielded ${n_models:-0} MODEL records (expected 2); the 2nd input was dropped" >> test73_multi_extract.out
  exit 1
fi

# test74: --backend-model routing — the override must reach the resolved config.
# dry-run + show-config (no model download, fast): proves --backend-model is honored.
pdb2reaction opt -i r.pdb -q -1 --backend-model uma-m-1p1 --show-config --dry-run --out-dir test74_backend_model > test74_backend_model.out 2>&1
grep -Eq '^\[backend\] uma \(uma-m-1p1, fp32\)$' test74_backend_model.out || { echo "[smoke] FAIL test74: non-default backend model missing from resolved runtime summary" >> test74_backend_model.out; exit 1; }

# test75: mmCIF input crosses the internal PDB bridge for a real MLIP run and
# restores the original long chain / large residue identifiers in public CIF.
pdb2reaction opt -i r_complex.cif -q 0 --max-cycles 1 --thresh gau_loose --out-dir test75_opt_cif > test75_opt_cif.out 2>&1
test -s test75_opt_cif/final_geometry.pdb || { echo "[smoke] FAIL test75: final PDB missing" >> test75_opt_cif.out; exit 1; }
test -s test75_opt_cif/final_geometry.cif || { echo "[smoke] FAIL test75: final CIF missing" >> test75_opt_cif.out; exit 1; }
grep -q 'LONG_CHAIN' test75_opt_cif/final_geometry.cif || { echo "[smoke] FAIL test75: auth chain was not restored" >> test75_opt_cif.out; exit 1; }
grep -q '10001' test75_opt_cif/final_geometry.cif || { echo "[smoke] FAIL test75: auth residue number was not restored" >> test75_opt_cif.out; exit 1; }

# test76: chain + residue-name + residue-number selection remains exact after
# mmCIF normalization, and extract emits both its internal PDB and public CIF.
pdb2reaction extract -i r_complex.cif -c 'LONG_CHAIN:PRE:10001' -r 0.1 --no-add-linkh -o model_from_cif.pdb -v 0 > test76_extract_cif.out 2>&1
test -s model_from_cif.pdb || { echo "[smoke] FAIL test76: extracted PDB missing" >> test76_extract_cif.out; exit 1; }
test -s model_from_cif.cif || { echo "[smoke] FAIL test76: extracted CIF missing" >> test76_extract_cif.out; exit 1; }
python assert_cif_selection.py model_from_cif.pdb model_from_cif.cif >> test76_extract_cif.out 2>&1

# test77: runtime YAML precedence reaches a real calculation (config-only
# charge/model/precision, with CLI max-cycles overriding the YAML budget).
pdb2reaction opt -i r.pdb --config runtime_config.yaml --show-config --max-cycles 1 --out-json --out-dir test77_config > test77_config.out 2>&1
python assert_release_result.py opt-config test77_config --expected-charge -1 --expected-model uma-s-1p2 --expected-max-cycles 1 --expected-precision fp64 >> test77_config.out 2>&1

# test78: a user ASE calculator that implements energy/forces only supports
# finite-difference Hessian consumers.
pdb2reaction sp -i r.pdb -q -1 --calc-file harmonic_calc.py --hess --hessian-calc-mode FiniteDifference --out-json --out-dir test78_custom > test78_custom.out 2>&1
python assert_release_result.py sp-hessian test78_custom >> test78_custom.out 2>&1

# test79: explicit Analytical + workers>1 is a required rejection contract.
rc=0
pdb2reaction sp -i r.pdb -q -1 --hess --hessian-calc-mode Analytical --workers 2 --out-dir test79_workers > test79_workers.out 2>&1 || rc=$?
if [ "$rc" -ne 1 ] || ! grep -Fq "Analytical Hessian cannot be combined with UMA workers>1: the parallel predictor exposes no autograd model. Use workers=1 or select hessian_calc_mode='FiniteDifference'." test79_workers.out || [ -e test79_workers/hessian.npy ]; then
  echo "[smoke] FAIL test79: Analytical + workers>1 was not rejected exactly" >> test79_workers.out
  exit 1
fi

# test80: force the never-stop energy-convergence guard to fire, then require
# both IRC branches to continue past it.
pdb2reaction irc -i ts.pdb -q 0 --config never_stop_config.yaml --never-stop --max-cycles 2 --out-json --out-dir test80_irc_never_stop > test80_irc_never_stop.out 2>&1
python assert_release_result.py irc-never-stop test80_irc_never_stop >> test80_irc_never_stop.out 2>&1

# Numerical analytical-vs-FD agreement for every backend installed in the
# default strict environment. MACE/AIMNet2 use this same required wrapper in
# their dependency-isolated cluster environments.
bash run_backend_hessian.sh uma orb > backend_hessian.out 2>&1

echo "[smoke] PASS: required core-GPU and structure-I/O lane completed with zero skips."
