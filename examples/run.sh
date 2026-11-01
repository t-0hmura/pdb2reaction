#!/usr/bin/env bash

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

pdb2reaction all -i "${script_dir}/1.R.pdb" "${script_dir}/3.P.pdb" -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --tsopt --thermo --out-dir result_mep

pdb2reaction all -i "${script_dir}/1.R.pdb" -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --scan-lists '[("CS1 SAM 320","GPP 321 C7",1.60)]' '[("GPP 321 H11","GLU 186 OE2",0.90)]' --tsopt --thermo --out-dir result_scan

# COMT catechol O-methyltransferase (examples/comt): endpoint mode with a 4.0 Å pocket.
pdb2reaction all -i "${script_dir}/comt/1.R.pdb" "${script_dir}/comt/3.P.pdb" -c 'CAT,SAM,MG' -l 'CAT:-1,SAM:0' -r 4.0 --tsopt --thermo --out-dir result_comt
