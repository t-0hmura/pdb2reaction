# Basic opt command tests
pdb2reaction opt -i r.pdb -q -1 --dump True --out-dir 'test1' > test1.out 2>&1
pdb2reaction opt -i r.pdb -q -1 --dump True --out-dir 'test2' --args-yaml input.yaml > test2.out 2>&1

# Input format tests (pdb/xyz/gjf) - test3 keeps --tsopt/--thermo for verification
pdb2reaction -i r.pdb p.pdb -q -1 --out-dir 'test3' --tsopt True --thermo True > test3.out 2>&1
pdb2reaction -i r.xyz p.xyz -q -1 --out-dir 'test4' > test4.out 2>&1
pdb2reaction -i r.gjf p.gjf -q -1 --out-dir 'test5' > test5.out 2>&1
pdb2reaction -i r.gjf p.gjf --out-dir 'test6' > test6.out 2>&1

# Scan-lists tests (--tsopt/--thermo removed for efficiency)
pdb2reaction -i r.pdb -q -1 --out-dir 'test7' --scan-lists "[(1,5,1.4)]" > test7.out 2>&1
pdb2reaction -i r.gjf --out-dir 'test8' --scan-lists "[(1,5,1.4)]" > test8.out 2>&1
pdb2reaction -i r.xyz -q -1 --out-dir 'test9' --scan-lists "[(1,5,1.4)]" > test9.out 2>&1

# Complex system tests
pdb2reaction opt -i p_complex.pdb --ligand-charge 'PRE:-2' --dump True --out-dir 'test10' > test10.out 2>&1
pdb2reaction -i p_complex.pdb -c 'PRE' --ligand-charge 'PRE:-2' --scan-lists "[('PRE 8 C3','PRE 8 O1\'',1.4),('PRE 8 C1','PRE 8 C8',3.3)]" -r 5.0 --exclude-backbone False --out-dir 'test11' > test11.out 2>&1
pdb2reaction -i r_complex.pdb p_complex.pdb -c 'PRE' --ligand-charge 'PRE:-2' -r 5.0 --exclude-backbone False --out-dir 'test12' > test12.out 2>&1

# 2D/3D scan tests
pdb2reaction extract -i p_complex.pdb -c 'PRE' --ligand-charge 'PRE:-2' -r 5.0 --exclude-backbone False --include-H2O True -o p_complex_pocket.pdb
pdb2reaction scan2d -i p_complex_pocket.pdb --ligand-charge 'PRE:-2' --scan-lists "[('PRE 8 C3','PRE 8 O1\'',1.4,1.8),('PRE 8 C1','PRE 8 C8',3.0,3.4)]" --out-dir 'test13' > test13.out 2>&1
pdb2reaction scan3d -i p_complex_pocket.pdb --ligand-charge 'PRE:-2' --scan-lists "[('PRE 8 C3','PRE 8 O1\'',1.4,1.8),('PRE 8 C1','PRE 8 C8',3.0,3.4),('PRE 8 C1','PRE 8 C7',1.4,1.6)]" --out-dir 'test14' > test14.out 2>&1

# Multiple scan-lists and MEP-mode tests
pdb2reaction -i p_complex.pdb -c 'PRE' --ligand-charge 'PRE:-2' --scan-lists "[('PRE 8 C3','PRE 8 O1\'',1.4),('PRE 8 C1','PRE 8 C8',3.3)]" "[('PRE 8 C1','PRE 8 C7',3.5)]" -r 5.0 --exclude-backbone False --out-dir 'test15' > test15.out 2>&1
pdb2reaction -i r_complex.pdb p_complex.pdb -c 'PRE' --ligand-charge 'PRE:-2' -r 5.0 --exclude-backbone False --mep-mode dmf --out-dir 'test16' > test16.out 2>&1

# TS input tests - test17 keeps --tsopt/--thermo for verification
pdb2reaction -i ts.pdb -q 0 --out-dir 'test17' --tsopt True --thermo True > test17.out 2>&1
pdb2reaction -i ts.pdb -q 0 --out-dir 'test18' --tsopt True --opt-mode light > test18.out 2>&1
pdb2reaction -i ts.gjf --out-dir 'test19' --tsopt True > test19.out 2>&1

# Additional subcommand tests
pdb2reaction dft -i h2.gjf --out-dir 'test20' --engine cpu > test20.out 2>&1
pdb2reaction freq -i r.pdb -q -1 --out-dir 'test21' > test21.out 2>&1
pdb2reaction tsopt -i ts.pdb -q 0 --out-dir 'test22' > test22.out 2>&1

# TSOPT for complex systems
pdb2reaction -i ts_complex.pdb -q -1 --out-dir 'test23' --tsopt True --opt-mode heavy > test23.out 2>&1
pdb2reaction -i ts_complex.pdb -q -1 --out-dir 'test24' --tsopt True --opt-mode light > test24.out 2>&1
pdb2reaction tsopt -i ts_complex.pdb -q -1 --out-dir 'test25' --opt-mode light > test25.out 2>&1