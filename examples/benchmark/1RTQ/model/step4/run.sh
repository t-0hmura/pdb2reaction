#!/bin/bash
pdb2reaction all -i 8.IM3.xyz 10.P.xyz -q 0 --config config.yaml --tsopt --thermo --dft --dft-func-basis "wb97m-v/def2-svp" --out-dir result > all.out 2>&1
pdb2reaction all -i 8.IM3.xyz 10.P.xyz -q 0 --config config.yaml --tsopt --thermo --dft --dft-func-basis "wb97m-v/def2-svp" --no-refine-path --out-dir result_nrp > all_nrp.out 2>&1
pdb2reaction all -i 9.TS4.xyz -q 0 --config config.yaml --tsopt --thermo --dft --dft-func-basis "wb97m-v/def2-svp" --out-dir result_tsonly > all_tsonly.out 2>&1
pdb2reaction all -i 8.IM3.xyz 10.P.xyz -q 0 --config config.yaml --tsopt --thermo --dft --dft-func-basis "wb97m-v/def2-svp" --no-refine-path --out-dir result_ana --hessian-calc-mode Analytical > all_ana.out 2>&1
