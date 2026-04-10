#!/usr/bin/env bash
#PBS -N smoke_p2r
#PBS -q default
#PBS -l nodes=1:ppn=8:gpus=1,mem=60GB,walltime=4:00:00
#PBS -o /dev/null
#PBS -e /dev/null
cd "${PBS_O_WORKDIR:-.}"
if [ -n "${PBS_O_WORKDIR:-}" ]; then
  . /home/apps/Modules/init/profile.sh
  module load cuda/12.9
  source /data2/tohmura/miniconda3/etc/profile.d/conda.sh
  conda activate p2r
  export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
fi
set -euo pipefail
# pdb2reaction smoke tests — GPU required
# Coverage: all subcommands, input formats, --ref-pdb, solvent, dry-run, utilities.

# Clean previous results
rm -rf test* p_complex_model.pdb model_r.pdb r_complex_elem.pdb r_complex_fixalt.pdb

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
pdb2reaction irc -i ts.pdb -q 0 --max-cycles 3 --out-dir test5 > test5.out 2>&1

# test6: dft (lightweight: hf/sto-3g, cpu)
pdb2reaction dft -i h2.gjf --func-basis 'hf/sto-3g' --grid-level 0 --conv-tol 1e-5 --max-cycle 40 --engine cpu --out-dir test6 > test6.out 2>&1

# test7: scan (1D)
pdb2reaction scan -i r.pdb -q -1 --scan-lists "[(1,5,1.4)]" --max-step-size 2.0 --relax-max-cycles 3 --no-preopt --no-endopt --out-dir test7 > test7.out 2>&1

# test8: scan2d (extract model first)
pdb2reaction extract -i p_complex.pdb -c 'PRE' --ligand-charge 'PRE:-2' -r 5.0 --no-exclude-backbone --include-h2o -o p_complex_model.pdb
pdb2reaction scan2d -i p_complex_model.pdb --ligand-charge 'PRE:-2' --scan-lists "[('PRE 8 C3','PRE 8 O1\'',1.4,1.8),('PRE 8 C1','PRE 8 C8',3.0,3.4)]" --max-step-size 2.0 --relax-max-cycles 100 --thresh gau_loose --out-dir test8 > test8.out 2>&1

# test9: scan3d
pdb2reaction scan3d -i p_complex_model.pdb --ligand-charge 'PRE:-2' --scan-lists "[('PRE 8 C3','PRE 8 O1\'',1.4,1.8),('PRE 8 C1','PRE 8 C8',3.0,3.4),('PRE 8 C1','PRE 8 C7',1.4,1.6)]" --max-step-size 2.0 --relax-max-cycles 100 --thresh gau_loose --out-dir test9 > test9.out 2>&1

# test10: path-opt (gsm)
pdb2reaction path-opt -i r.pdb p.pdb -q -1 --max-nodes 5 --max-cycles 5 --no-preopt --no-climb --out-dir test10 > test10.out 2>&1

# test11: path-search
pdb2reaction path-search -i r.pdb p.pdb -q -1 --max-nodes 5 --max-cycles 5 --out-dir test11 > test11.out 2>&1

# --- Input format tests (all command) ---

# test12: all (pdb+pdb, --no-refine-path for single-pass path-opt)
pdb2reaction -i r.pdb p.pdb -q -1 --refine-path false --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --out-dir test12 > test12.out 2>&1

# test13: all (xyz+xyz, --ref-pdb for PDB conversion, --no-refine-path for single-pass path-opt)
pdb2reaction -i r.xyz p.xyz -q -1 --ref-pdb r.pdb --refine-path false --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --out-dir test13 > test13.out 2>&1

# test14: all (gjf+gjf, --no-refine-path for single-pass path-opt)
pdb2reaction -i r.gjf p.gjf --refine-path false --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --out-dir test14 > test14.out 2>&1

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
pdb2reaction -i r_complex.pdb p_complex.pdb -c 'PRE' --ligand-charge 'PRE:-2' -r 5.0 --no-exclude-backbone --refine-path false --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --out-dir test19 > test19.out 2>&1

# --- TSOPT-only mode ---

# test20: all (ts input, --tsopt + --thermo)
pdb2reaction -i ts.pdb -q 0 --tsopt --opt-mode-post grad --thermo --max-cycles 100 --thresh gau --thresh-post gau_loose --out-dir test20 > test20.out 2>&1

# test21: all (ts input, --tsopt, opt-mode hess)
pdb2reaction -i ts.pdb -q 0 --tsopt --opt-mode-post hess --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --out-dir test21 > test21.out 2>&1

# --- MEP mode ---

# test22: all (pdb+pdb, mep-mode dmf, --no-refine-path for single-pass path-opt)
pdb2reaction -i r.pdb p.pdb -q -1 --mep-mode dmf --refine-path false --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --out-dir test22 > test22.out 2>&1

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
pdb2reaction scan2d -i p_complex_model.pdb --ligand-charge 'PRE:-2' -s scan2d_spec.yaml --max-step-size 2.0 --relax-max-cycles 100 --thresh gau_loose --out-dir test38 > test38.out 2>&1

# --- Solvent correction (requires xTB) ---

# test39: opt (solvent water) — skip if xtb is not available
if command -v xtb &>/dev/null; then
  pdb2reaction opt -i r.pdb -q -1 --opt-mode grad --max-cycles 3 --thresh gau_loose --solvent water --out-dir test39 > test39.out 2>&1
else
  echo "SKIP test39: xtb not found" > test39.out
fi

# --- dist-freeze ---

# test40: opt --dist-freeze --dry-run (inline 3-tuple + 2-tuple)
pdb2reaction opt -i r.pdb -q -1 --dist-freeze "[(1,2,1.5),(3,4)]" --dry-run --out-dir test40 > test40.out 2>&1

# --- refine-path ---

# test41: all (pdb+pdb, --refine-path — now the default, kept explicit for clarity)
pdb2reaction -i r.pdb p.pdb -q -1 --refine-path true --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --out-dir test41 > test41.out 2>&1
