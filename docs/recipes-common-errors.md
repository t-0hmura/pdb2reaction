# Common Error Recipes

Use this page when you know the symptom but do not know which subcommand page to open first.
For full details, keep [Troubleshooting](troubleshooting.md) open in parallel.

## Quick routing

| Symptom | Start here | Then read |
| --- | --- | --- |
| Missing element columns / extraction aborts | `add-elem-info` on the original PDB | [Troubleshooting](troubleshooting.md) |
| "Charge is required" errors | Set `-q/--charge` or `-l/--ligand-charge` explicitly | [Troubleshooting](troubleshooting.md) |
| Energies/states look wrong after a run | Re-check charge/multiplicity policy in CLI conventions | [Troubleshooting](troubleshooting.md) |
| DMF mode import errors (`cyipopt`) | Run `conda install -c conda-forge cyipopt` | [Troubleshooting](troubleshooting.md) |
| TSOPT/IRC does not converge | For LBFGS/Dimer: adjust `max_step`. For RFO/RS-I-RFO: adjust `trust_radius`/`trust_min`/`trust_max`. Increase cycles, validate TS quality first | [Troubleshooting](troubleshooting.md) |
| CUDA/GPU runtime mismatch | Verify `torch.cuda.is_available()` and CUDA build pairing | [Troubleshooting](troubleshooting.md) |
| Plot export failures | Run `plotly_get_chrome -y` to install headless Chrome | [Troubleshooting](troubleshooting.md) |

## Recipe 1: Extraction fails before MEP starts

Signal:  
- Errors mention missing element symbols, atom-count mismatch, or empty active site models (binding pockets).  

First checks:  
- Confirm all inputs are prepared by the same workflow and atom ordering is consistent.  
- Ensure element columns are present before running `extract` or `all`.  

Typical fix path:  
- Repair elements with `pdb2reaction add-elem-info -i input.pdb -o input_fixed.pdb` -> rerun extraction -> confirm active site model size (`--radius`) and residue inclusion (`--selected-resn`).  

## Recipe 2: Charge/spin validation fails

Signal:  
- Errors mention unresolved charge, especially for non-`.gjf` inputs.  

First checks: 
- Ensure net charge and multiplicity are physically correct for the target state.  
- If using residue maps, validate each residue key in `--ligand-charge`.  
- Verify the resolution rules in [CLI Conventions](cli-conventions.md) when results look physically inconsistent.  

Typical fix path:  
- Prefer explicit `-q/--charge` or `-l/--ligand-charge` and `-m` for critical runs, then retry scan/path/tsopt.   

## Recipe 3: Environment blockers

Signal:  
- DMF import failures, CUDA mismatch, or unavailable plotting backends.

First checks:  
- Confirm optional dependencies are installed in the active env.
- Validate GPU visibility and PyTorch CUDA compatibility.

Typical fix path:
- Repair environment first, verify with `pdb2reaction --version` and `python -c "import torch; print(torch.cuda.is_available())"`, then rerun with `--dry-run` before full execution.

## Recipe 4: Convergence and post-processing failures

Signal:  
- TSOPT stalls, IRC branches look unstable, or MEP refinement stops unexpectedly.  

First checks:  
- Confirm TS candidate quality: exactly one imaginary frequency with |ν| ≥ 100 cm⁻¹, and the corresponding imaginary mode shows displacement along the reaction coordinate.  
- For LBFGS/Dimer: reduce `max_step`. For RFO/RS-I-RFO: reduce `trust_radius`/`trust_min`/`trust_max`. Increase cycle limits.  

Typical fix path:
- Run a smaller diagnostic case, tune thresholds/step sizes, then scale back up.
