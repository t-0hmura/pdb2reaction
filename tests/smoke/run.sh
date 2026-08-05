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
pdb2reaction path-opt -i r.pdb p.pdb -q -1 --max-nodes 5 --max-cycles 5 --thresh-gsm gau_loose --no-preopt --no-climb --out-dir test10 > test10.out 2>&1

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

# test39: opt (solvent water)
pdb2reaction opt -i r.pdb -q -1 --opt-mode grad --max-cycles 3 --thresh gau_loose --solvent water --out-dir test39 > test39.out 2>&1

# --- dist-freeze ---

# test40: opt --dist-freeze --dry-run (inline 3-tuple + 2-tuple)
pdb2reaction opt -i r.pdb -q -1 --dist-freeze "[(1,2,1.5),(3,4)]" --dry-run --out-dir test40 > test40.out 2>&1

# --- refine-path ---

# test41: all (pdb+pdb, --refine-path = recursive path_search opt-in; the default is now single-pass path-opt)
pdb2reaction -i r.pdb p.pdb -q -1 --refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --out-dir test41 > test41.out 2>&1

# --- Opt-in TS and IRC methods ---

# test42: tsopt --opt-mode trim (Helgaker trust-region image-min TS opt)
pdb2reaction tsopt -i ts.pdb -q 0 --opt-mode trim --max-cycles 5 --thresh gau_loose --out-dir test42 > test42.out 2>&1

# test43: tsopt --opt-mode rsprfo (Banerjee restricted-step P-RFO TS opt)
pdb2reaction tsopt -i ts.pdb -q 0 --opt-mode rsprfo --max-cycles 5 --thresh gau_loose --out-dir test43 > test43.out 2>&1

# test44: irc --irc-pos-def (PSD-Hessian convergence guard)
pdb2reaction irc -i ts.pdb -q 0 --max-cycles 3 --irc-pos-def --out-dir test44 > test44.out 2>&1

# test45: opt --print-every 3 (print-throttling check, no behavior change)
pdb2reaction opt -i r.pdb -q -1 --opt-mode hess --max-cycles 5 --thresh gau_loose --print-every 3 --out-dir test45 > test45.out 2>&1

# --- Determinism gate ---

# test46: `all` pipeline determinism GATE (`--deterministic`).
# Runs the full pipeline twice with identical inputs / args + `--deterministic`
# and REQUIRES their geometry artifacts (*.pdb / *.xyz) to be bit-identical --
# that manifest is what the gate compares. Default (non-deterministic) GPU
# runs carry ~ULP scatter/atomic non-determinism and are not asserted here;
# `--deterministic` enables torch deterministic algorithms and MUST be
# bit-reproducible, so any drift is a real regression and fails the smoke.
det_args="-i r.pdb p.pdb -q -1 --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --deterministic"
pdb2reaction $det_args --out-dir test46_a > test46_a.out 2>&1
pdb2reaction $det_args --out-dir test46_b > test46_b.out 2>&1
find test46_a -type f \( -name "*.pdb" -o -name "*.xyz" \) -printf '%P\n' | LC_ALL=C sort > test46_a.manifest
find test46_b -type f \( -name "*.pdb" -o -name "*.xyz" \) -printf '%P\n' | LC_ALL=C sort > test46_b.manifest
if ! cmp -s test46_a.manifest test46_b.manifest; then
  echo "[smoke] FAIL test46: deterministic runs produced different file manifests" > test46.out
  comm -3 test46_a.manifest test46_b.manifest >> test46.out
  cat test46.out
  exit 1
fi
total=$(wc -l < test46_a.manifest)
if [ "$total" -eq 0 ]; then
  echo "[smoke] FAIL test46: deterministic gate found no PDB/XYZ artifacts" > test46.out
  cat test46.out
  exit 1
fi
drifted=0
while IFS= read -r rel; do
  if ! cmp -s "test46_a/$rel" "test46_b/$rel"; then
    drifted=$((drifted + 1))
    echo "DRIFT: $rel" >> test46.out
  fi
done < test46_a.manifest
echo "[det_check] compared $total PDB/XYZ files; $drifted differ" >> test46.out
if [ "$drifted" -ne 0 ]; then
  echo "[smoke] FAIL test46: --deterministic runs differ" >> test46.out
  cat test46.out
  exit 1
fi

# --- --coord-type CLI plumbing, plus the single-point lanes (throttled, fast) ---

# test47: `all --coord-type cart` — explicit cart (== default), verifies CLI plumbing.
pdb2reaction all -i r.pdb p.pdb -q -1 --coord-type cart --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --out-dir test47 > test47.out 2>&1

# test48: `all --coord-type dlc` — DLC propagated to child opt / tsopt / path-opt stages.
pdb2reaction all -i r.pdb p.pdb -q -1 --coord-type dlc --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --out-dir test48 > test48.out 2>&1

# test49: `sp` (single-point) — energy + forces.
pdb2reaction sp -i r.pdb -q -1 --out-dir test49 > test49.out 2>&1

# test50: `sp --hess` — energy + forces + Hessian (default FiniteDifference).
pdb2reaction sp -i r.pdb -q -1 --hess --out-dir test50 > test50.out 2>&1

# --- Full-pipeline release-gate runs ---
# test51 is unthrottled: the canonical `all` flow with default convergence
# thresholds and no optimizer-cycle cap. It is the full-pipeline completion
# gate. test52 remains throttled because it only exercises the DLC path end to
# end.

# test51: full `all` cart — default thresh, no max-cycles cap.
pdb2reaction all -i r.pdb p.pdb -q -1 --out-dir test51 > test51.out 2>&1

# test52: `all` dlc — verifies the DLC code-path lights up end-to-end.
# Capped at max-cycles 5 + thresh gau_loose + --no-tsopt/thermo/dft because
# DLC GSM on this system needs many cycles to converge with default `gau`
# thresh, and the post-stages (TS / IRC / DFT) depend on a converged HEI
# from the MEP — silently broken structure handoff otherwise. test51 (cart)
# keeps the no-cap default-behaviour check.
pdb2reaction all -i r.pdb p.pdb -q -1 --coord-type dlc --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --out-dir test52 > test52.out 2>&1

# --- Per-stage internal-coordinate code paths ---
# Each test scoped at a 2-3 cycle cap + thresh gau_loose to exercise the
# --coord-type code path without long convergence. Frequency analysis remains
# Cartesian because the only thing an internal-coordinate flag would add is
# an unrelated B-matrix construction:
#   - freq: vibrational analysis uses cart Hessian regardless of
#     --coord-type; the only thing the flag adds is an SVD on the internal
#     B-matrix that occasionally fails on cuSOLVER. freq + cart is already
#     covered by test4.

# test52a: opt --coord-type dlc
pdb2reaction opt -i r.pdb -q -1 --coord-type dlc --max-cycles 3 --thresh gau_loose --out-dir test52a_opt_dlc > test52a_opt_dlc.out 2>&1

# test52b-e: Hessian TS optimization must accept every documented coordinate
# system. The DLC case includes frozen atoms and therefore exercises the
# partial-Cartesian-Hessian -> internal-coordinate handoff regression.
pdb2reaction tsopt -i ts.pdb -q 0 --opt-mode hess --coord-type cart --max-cycles 2 --thresh gau_loose --out-dir test52b_ts_cart > test52b_ts_cart.out 2>&1
pdb2reaction tsopt -i ts.pdb -q 0 --opt-mode hess --coord-type redund --max-cycles 2 --thresh gau_loose --out-dir test52c_ts_redund > test52c_ts_redund.out 2>&1
pdb2reaction tsopt -i ts.pdb -q 0 --opt-mode hess --coord-type dlc --freeze-atoms 1,3,5 --max-cycles 2 --thresh gau_loose --out-dir test52d_ts_dlc_freeze > test52d_ts_dlc_freeze.out 2>&1
pdb2reaction tsopt -i ts.pdb -q 0 --opt-mode hess --coord-type tric --max-cycles 2 --thresh gau_loose --out-dir test52e_ts_tric > test52e_ts_tric.out 2>&1

# test52d_scan_dlc: scan --coord-type dlc (1D).  The bare `test52d` label belongs
# to the tsopt lane above; this one is named after its own out-dir.
pdb2reaction scan -i r.pdb -q -1 --coord-type dlc --scan-lists "[(1,5,1.4)]" --max-step-size 2.0 --relax-max-cycles 3 --no-preopt --no-endopt --out-dir test52d_scan_dlc > test52d_scan_dlc.out 2>&1

# test52g: opt --coord-type redund
pdb2reaction opt -i r.pdb -q -1 --coord-type redund --max-cycles 3 --thresh gau_loose --out-dir test52g_opt_redund > test52g_opt_redund.out 2>&1

# test52j: scan --coord-type redund
pdb2reaction scan -i r.pdb -q -1 --coord-type redund --scan-lists "[(1,5,1.4)]" --max-step-size 2.0 --relax-max-cycles 3 --no-preopt --no-endopt --out-dir test52j_scan_redund > test52j_scan_redund.out 2>&1

# test52k: opt --coord-type tric
pdb2reaction opt -i r.pdb -q -1 --coord-type tric --max-cycles 3 --thresh gau_loose --out-dir test52k_opt_tric > test52k_opt_tric.out 2>&1

# --- Multi-mode flag code-path verify (single-stage) ---

# test52m: opt --precision fp64 (alternate precision)
pdb2reaction opt -i r.pdb -q -1 --precision fp64 --max-cycles 3 --thresh gau_loose --out-dir test52m_opt_fp64 > test52m_opt_fp64.out 2>&1

# test52n: opt --precision fp32 (explicit fp32 dispatch; default is fp32, this pins the explicit path alongside test52m fp64)
pdb2reaction opt -i r.pdb -q -1 --precision fp32 --max-cycles 3 --thresh gau_loose --out-dir test52n_opt_fp32 > test52n_opt_fp32.out 2>&1

# test53: full `all` with `--backend orb` — exercises the non-default MLIP backend.
pdb2reaction all -i r.pdb p.pdb -q -1 --backend orb --out-dir test53 > test53.out 2>&1
python assert_release_result.py provenance test53 --expected-backend orb --expected-model orb_v3_conservative_omol --expected-precision fp64 >> test53.out 2>&1

# ---- Subcommand-specific regression coverage ----
# test54: opt --dist-freeze + --bias-k (harmonic restraint actually applied at runtime, not dry-run)
pdb2reaction opt -i r.pdb -q -1 --dist-freeze "[(1,2,1.5)]" --bias-k 150 --max-cycles 3 --thresh gau_loose --out-dir test54_opt_biask > test54_opt_biask.out 2>&1

# test55: opt --freeze-atoms (explicit user-frozen DOF, distinct from --freeze-links)
pdb2reaction opt -i r.pdb -q -1 --freeze-atoms '1,3,5' --max-cycles 3 --thresh gau_loose --no-flatten --out-dir test55_opt_freeze > test55_opt_freeze.out 2>&1

# test56: freq --hessian-calc-mode Analytical (workflow analytical path)
pdb2reaction freq -i r.pdb -q -1 --hessian-calc-mode Analytical --max-write 3 --out-dir test56_freq_anahess > test56_freq_anahess.out 2>&1

# test57: irc --hessian-calc-mode analytical (IRC-initiating Hessian path)
pdb2reaction irc -i ts.pdb -q 0 --hessian-calc-mode analytical --max-cycles 2 --out-dir test57_irc_anahess > test57_irc_anahess.out 2>&1

# test58: scan --opt-mode hess (RFO relaxation vs default LBFGS)
pdb2reaction scan -i r.pdb -q -1 --opt-mode hess --scan-lists "[(1,5,1.4)]" --max-step-size 2.0 --relax-max-cycles 2 --no-preopt --no-endopt --out-dir test58_scan_hess > test58_scan_hess.out 2>&1

# test59: scan2d --opt-mode hess (RFO per-grid relaxation).  Start from a
# converged test8 point and use the smallest genuine 2D grid (2 × 2).  The
# starting geometry is whatever test8 relaxed to, so a corner of this tiny grid
# may stay off its target distance and leave either no eligible points or too
# few points to interpolate.  A full 2 × 2 and either controlled outcome prove
# the RFO per-grid branch ran, while any other failure is fatal.
rc=0
pdb2reaction scan2d -i test8/grid/point_i140_j300.pdb --ligand-charge 'PRE:-2' --freeze-atoms "$P_COMPLEX_MODEL_FREEZE_ATOMS" --scan-lists "[('PRE 8 C3','PRE 8 O1\'',1.40,1.41),('PRE 8 C1','PRE 8 C8',3.00,3.01)]" --opt-mode hess --max-step-size 2.0 --relax-max-cycles 100 --thresh gau_loose --out-dir test59_scan2d_hess > test59_scan2d_hess.out 2>&1 || rc=$?
test -s test59_scan2d_hess/surface.csv || { echo "[smoke] FAIL test59: surface.csv missing" >> test59_scan2d_hess.out; exit 1; }
grep -Fq '[hessian] Completed FiniteDifference Hessian:' test59_scan2d_hess.out || { echo "[smoke] FAIL test59: Hessian branch did not run" >> test59_scan2d_hess.out; exit 1; }
if grep -Fq 'Traceback (most recent call last)' test59_scan2d_hess.out; then
  echo "[smoke] FAIL test59: unexpected traceback" >> test59_scan2d_hess.out
  exit 1
fi
if [ "$rc" -eq 0 ]; then
  test -s test59_scan2d_hess/scan2d_map.png || { echo "[smoke] FAIL test59: successful run omitted the 2D plot" >> test59_scan2d_hess.out; exit 1; }
elif [ "$rc" -eq 1 ]; then
  if grep -Fq 'No converged finite grid point with a written geometry is available' test59_scan2d_hess.out; then
    grep -Fq '[plot] No finite data for plotting.' test59_scan2d_hess.out \
      || { echo "[smoke] FAIL test59: no-data exit omitted the plotting diagnostic" >> test59_scan2d_hess.out; exit 1; }
  else
    grep -Fq '[plot] ERROR: A 2D energy surface requires at least three non-collinear' test59_scan2d_hess.out \
      || { echo "[smoke] FAIL test59: exit 1 was not a controlled insufficient-support outcome" >> test59_scan2d_hess.out; exit 1; }
  fi
else
  echo "[smoke] FAIL test59: unexpected exit status $rc" >> test59_scan2d_hess.out
  exit 1
fi

# test60: scan3d --opt-mode hess (RFO per-grid relaxation).  A deliberately
# short budget may either converge or end with the controlled no-eligible-data
# outcome; both prove the Hessian branch ran, while any other failure is fatal.
rc=0
pdb2reaction scan3d -i r.pdb -q -1 --scan-lists "[(1,5,2.20,2.21),(1,6,1.75,1.76),(2,3,1.65,1.66)]" --opt-mode hess --max-step-size 0.1 --relax-max-cycles 3 --no-preopt --thresh gau_loose --out-dir test60_scan3d_hess > test60_scan3d_hess.out 2>&1 || rc=$?
test -s test60_scan3d_hess/surface.csv || { echo "[smoke] FAIL test60: surface.csv missing" >> test60_scan3d_hess.out; exit 1; }
grep -Fq '[hessian] Completed FiniteDifference Hessian:' test60_scan3d_hess.out || { echo "[smoke] FAIL test60: Hessian branch did not run" >> test60_scan3d_hess.out; exit 1; }
if grep -Fq 'Traceback (most recent call last)' test60_scan3d_hess.out; then
  echo "[smoke] FAIL test60: unexpected traceback" >> test60_scan3d_hess.out
  exit 1
fi
if [ "$rc" -eq 0 ]; then
  test -s test60_scan3d_hess/scan3d_density.html || { echo "[smoke] FAIL test60: successful run omitted the 3D plot" >> test60_scan3d_hess.out; exit 1; }
elif [ "$rc" -eq 1 ]; then
  grep -Fq 'No converged finite grid point with a written geometry is available' test60_scan3d_hess.out \
    && grep -Fq '[plot] No finite data for plotting.' test60_scan3d_hess.out \
    || { echo "[smoke] FAIL test60: exit 1 was not the controlled nonconvergence outcome" >> test60_scan3d_hess.out; exit 1; }
else
  echo "[smoke] FAIL test60: unexpected exit status $rc" >> test60_scan3d_hess.out
  exit 1
fi

# test61: scan3d --precision fp64 option plumbing.  Runtime fp64 dispatch is
# exercised by test52m; this pins the scan3d-specific Click/config path.
pdb2reaction scan3d -i r.pdb -q -1 --scan-lists "[(1,5,2.20,2.21),(1,6,1.75,1.76),(2,3,1.65,1.66)]" --precision fp64 --dry-run --out-dir test61_scan3d_fp64 > test61_scan3d_fp64.out 2>&1
grep -Fq '[scan3d] --dry-run: input, charge/spin parity, and --scan-lists parse OK.' test61_scan3d_fp64.out || { echo "[smoke] FAIL test61: fp64 option did not reach scan3d dry-run" >> test61_scan3d_fp64.out; exit 1; }

# test62: path-opt --mep-mode dmf (Direct Max Flux at the subcommand level)
pdb2reaction path-opt -i r.pdb p.pdb -q -1 --mep-mode dmf --max-nodes 5 --max-cycles 3 --thresh-dmf middle --no-preopt --no-climb --out-dir test62_pathopt_dmf > test62_pathopt_dmf.out 2>&1

# test63: path-opt --coord-type dlc (p2r keeps DLC for pure-MLIP)
pdb2reaction path-opt -i r.pdb p.pdb -q -1 --coord-type dlc --max-nodes 5 --max-cycles 3 --no-preopt --no-climb --out-dir test63_pathopt_dlc > test63_pathopt_dlc.out 2>&1

# test64: path-search --mep-mode dmf
pdb2reaction path-search -i r.pdb p.pdb -q -1 --mep-mode dmf --max-nodes 5 --max-cycles 3 --no-preopt --out-dir test64_ps_dmf > test64_ps_dmf.out 2>&1

# test65: path-search --opt-mode hess (RFO single-structure preopt; keep preopt ON)
pdb2reaction path-search -i r.pdb p.pdb -q -1 --opt-mode hess --workers 1 --max-nodes 5 --max-cycles 3 --out-dir test65_ps_hess > test65_ps_hess.out 2>&1

# test66: required positive MEP -> TSopt -> IRC -> thermo handoff.
# --tsopt-max-cycles must cover the opt-in --flatten repair (which draws from this
# global budget); this system converges to n_imag=1 with no flatten iterations, but
# 200 keeps the budget above flatten_max_iter (=50, ~135 cycles) so a soft spectator
# mode could never be starved short of a clean saddle. Product default is 10000.
pdb2reaction all -i r.pdb p.pdb -q -1 --tsopt --thermo --flatten --irc-never-stop --max-cycles 5 --tsopt-max-cycles 200 --thresh gau_loose --thresh-post gau --out-dir test66_all_tsopt > test66_all_tsopt.out 2>&1
python assert_release_result.py all test66_all_tsopt --require-thermo >> test66_all_tsopt.out 2>&1
grep -Fq '[irc] Reusing cached TS Hessian from tsopt.' test66_all_tsopt.out || { echo '[smoke] FAIL test66: IRC did not report cached TS Hessian reuse' >> test66_all_tsopt.out; exit 1; }

# test67: all --scan-lists (single-PDB scan->path mode of `all`)
pdb2reaction all -i r.pdb -q -1 --scan-lists "[(1,5,1.4)]" --max-cycles 5 --thresh gau_loose --no-tsopt --no-thermo --no-dft --out-dir test67_all_scan > test67_all_scan.out 2>&1

# test68: opt --coord-type dlc on extracted complex model with explicit frozen boundary atoms
pdb2reaction opt -i p_complex_model.pdb --ligand-charge 'PRE:-2' --coord-type dlc --freeze-atoms "$P_COMPLEX_MODEL_FREEZE_ATOMS" --max-cycles 3 --thresh gau_loose --out-dir test68_opt_complex_freeze_dlc > test68_opt_complex_freeze_dlc.out 2>&1
python check_frozen_atoms.py p_complex_model.pdb test68_opt_complex_freeze_dlc/final_geometry.pdb "$P_COMPLEX_MODEL_FREEZE_ATOMS" test68 >> test68_opt_complex_freeze_dlc.out 2>&1

# test69: scan --coord-type dlc on extracted complex model with explicit frozen boundary atoms.
# Keep preopt enabled: this is the regression path that can trigger internal-coordinate rebuilds.
# Use a small non-reactive target so this tests coordinate integrity, not chemistry.
pdb2reaction scan -i p_complex_model.pdb --ligand-charge 'PRE:-2' --coord-type dlc --freeze-atoms "$P_COMPLEX_MODEL_FREEZE_ATOMS" --scan-lists "[('PRE 8 C3','PRE 8 O1\'',3.2)]" --max-step-size 0.2 --relax-max-cycles 3 --no-endopt --out-dir test69_scan_complex_freeze_dlc > test69_scan_complex_freeze_dlc.out 2>&1
python check_frozen_atoms.py p_complex_model.pdb test69_scan_complex_freeze_dlc/stage_01/result.pdb "$P_COMPLEX_MODEL_FREEZE_ATOMS" test69 >> test69_scan_complex_freeze_dlc.out 2>&1
if grep -q "Covalent-bond changes (start vs final): Yes" test69_scan_complex_freeze_dlc.out; then
  echo "[bond-check] test69: unexpected covalent-bond changes in non-reactive DLC+freeze scan" >> test69_scan_complex_freeze_dlc.out
  exit 1
fi

# --- refine-path opt-in (recursive path_search) extra coverage ---
# Since the `all` default is now single-pass path-opt, exercise the recursive
# path_search opt-in (`--refine-path`) across more input modes than test41.

# test70: all --scan-lists --refine-path (single-PDB scan -> recursive path_search)
pdb2reaction all -i r.pdb -q -1 --scan-lists "[(1,5,1.4)]" --refine-path --max-cycles 5 --thresh gau_loose --no-tsopt --no-thermo --no-dft --out-dir test70_rp_scan > test70_rp_scan.out 2>&1

# test71: all (complex multi-input) --refine-path (recursive path_search on the extracted model)
pdb2reaction all -i r_complex.pdb p_complex.pdb -c 'PRE' --ligand-charge 'PRE:-2' -r 5.0 --no-exclude-backbone --refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --out-dir test71_rp_complex > test71_rp_complex.out 2>&1

# test72: extract MULTI-INPUT via space-separated '-i a.pdb b.pdb' (one flag, two paths).
# Regression guard: a single -i with several space-separated paths must NOT drop the 2nd input.
# A single -o yields one multi-MODEL PDB, so both endpoints must appear (-> exactly 2 MODEL records).
pdb2reaction extract -i r_complex.pdb p_complex.pdb -c 'PRE' --ligand-charge 'PRE:-2' -r 5.0 --no-exclude-backbone -o model_multi.pdb > test72_multi_extract.out 2>&1
n_models=$(grep -c '^MODEL' model_multi.pdb 2>/dev/null)
if [ "${n_models:-0}" -ne 2 ]; then
  echo "[extract-multi] test72: space-separated '-i a b' yielded ${n_models:-0} MODEL records (expected 2); the 2nd input was dropped" >> test72_multi_extract.out
  exit 1
fi

# test73: --backend-model routing — the override must reach the resolved config.
# dry-run + show-config (no model download, fast): proves --backend-model is honored.
pdb2reaction opt -i r.pdb -q -1 --backend-model uma-m-1p1 --show-config --dry-run --out-dir test73_backend_model > test73_backend_model.out 2>&1
grep -Eq '^\[backend\] uma \(uma-m-1p1, fp32\)$' test73_backend_model.out || { echo "[smoke] FAIL test73: non-default backend model missing from resolved runtime summary" >> test73_backend_model.out; exit 1; }

# test74: mmCIF input crosses the internal PDB bridge for a real MLIP run and
# restores the original long chain / large residue identifiers in public CIF.
pdb2reaction opt -i r_complex.cif -q 0 --max-cycles 1 --thresh gau_loose --out-dir test74_opt_cif > test74_opt_cif.out 2>&1
test -s test74_opt_cif/final_geometry.pdb || { echo "[smoke] FAIL test74: final PDB missing" >> test74_opt_cif.out; exit 1; }
test -s test74_opt_cif/final_geometry.cif || { echo "[smoke] FAIL test74: final CIF missing" >> test74_opt_cif.out; exit 1; }
grep -q 'LONG_CHAIN' test74_opt_cif/final_geometry.cif || { echo "[smoke] FAIL test74: auth chain was not restored" >> test74_opt_cif.out; exit 1; }
grep -q '10001' test74_opt_cif/final_geometry.cif || { echo "[smoke] FAIL test74: auth residue number was not restored" >> test74_opt_cif.out; exit 1; }

# test75: chain + residue-name + residue-number selection remains exact after
# mmCIF normalization, and extract emits both its internal PDB and public CIF.
pdb2reaction extract -i r_complex.cif -c 'LONG_CHAIN:PRE:10001' -r 0.1 --no-add-linkh -o model_from_cif.pdb -v 0 > test75_extract_cif.out 2>&1
test -s model_from_cif.pdb || { echo "[smoke] FAIL test75: extracted PDB missing" >> test75_extract_cif.out; exit 1; }
test -s model_from_cif.cif || { echo "[smoke] FAIL test75: extracted CIF missing" >> test75_extract_cif.out; exit 1; }
python assert_cif_selection.py model_from_cif.pdb model_from_cif.cif >> test75_extract_cif.out 2>&1

# test76: runtime YAML precedence reaches a real calculation (config-only
# charge/model/precision, with CLI max-cycles overriding the YAML budget).
pdb2reaction opt -i r.pdb --config runtime_config.yaml --show-config --max-cycles 1 --out-json --out-dir test76_config > test76_config.out 2>&1
python assert_release_result.py opt-config test76_config --expected-charge -1 --expected-model uma-s-1p2 --expected-max-cycles 1 --expected-precision fp64 >> test76_config.out 2>&1

# test77: a user ASE calculator that implements energy/forces only supports
# finite-difference Hessian consumers.
pdb2reaction sp -i r.pdb -q -1 --calc-file harmonic_calc.py --hess --hessian-calc-mode FiniteDifference --out-json --out-dir test77_custom > test77_custom.out 2>&1
python assert_release_result.py sp-hessian test77_custom >> test77_custom.out 2>&1

# test78: explicit Analytical + workers>1 is a required rejection contract.
rc=0
pdb2reaction sp -i r.pdb -q -1 --hess --hessian-calc-mode Analytical --workers 2 --out-dir test78_workers > test78_workers.out 2>&1 || rc=$?
if [ "$rc" -ne 1 ] || ! grep -Fq "Analytical Hessian cannot be combined with UMA workers>1: the parallel predictor exposes no autograd model. Use workers=1 or select hessian_calc_mode='FiniteDifference'." test78_workers.out || [ -e test78_workers/hessian.npy ]; then
  echo "[smoke] FAIL test78: Analytical + workers>1 was not rejected exactly" >> test78_workers.out
  exit 1
fi

# test79: force the never-stop energy-convergence guard to fire, then require
# both IRC branches to continue past it.
pdb2reaction irc -i ts.pdb -q 0 --config never_stop_config.yaml --never-stop --max-cycles 2 --out-json --out-dir test79_irc_never_stop > test79_irc_never_stop.out 2>&1
python assert_release_result.py irc-never-stop test79_irc_never_stop >> test79_irc_never_stop.out 2>&1
if grep -Fq '[irc] IRC stopped after only a few frames' test79_irc_never_stop.out; then
  echo '[smoke] FAIL test79: cycle-cap completion was reported as early IRC termination' >> test79_irc_never_stop.out
  exit 1
fi

# Numerical analytical-vs-FD agreement for every backend installed in the
# default strict environment. MACE/AIMNet2 use this same required wrapper in
# their dependency-isolated cluster environments.
bash run_backend_hessian.sh uma orb > backend_hessian.out 2>&1

echo "[smoke] PASS: required core-GPU and structure-I/O lane completed with zero skips."
