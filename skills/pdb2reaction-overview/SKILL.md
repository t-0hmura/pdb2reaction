---
name: pdb2reaction-overview
description: "Orientation, TS strategy, and output reading for pdb2reaction, a PDB-native toolkit for MLIP reaction-path calculations on enzyme active-site clusters. SKILL.md first picks the `all` mode (endpoint MEP, scan, or TS-only) from the available structures, then covers stage-by-stage runs, how to judge each stage, and where the source code lives; ts-strategy.md covers precision, routes to a TS candidate and retries when n_imag is wrong, product-start scans, staged vs concerted scans, and controlled comparisons; outputs.md covers summary.json, R/TS/P paths, bond changes, energy diagrams, and failed runs. TRIGGER on first-touch questions, choosing an all mode or workflow, building or debugging a TS candidate, reading summary.json, extracting barriers or Gibbs energies, or locating code. SKIP for one subcommand (pdb2reaction-cli), structure files, charge, or cluster building (pdb2reaction-model-setup), install or CUDA (pdb2reaction-install), and job scripts (pdb2reaction-hpc)."
---

# pdb2reaction

`pdb2reaction all` picks its mode from the inputs: two or more structures in reaction order → endpoint MEP (`all-endpoint-mep.md`); one structure with `-s` → scan (`all-scan-list.md`); one TS candidate with `--tsopt` → TS-only (`all-ts-only.md`). Add `-c` to cut the cluster and `--tsopt --thermo` for the TS, IRC, and Gibbs energies; run the stages one by one when you want to judge each result first.

## Pick an all mode

| You have | Mode | Read |
|---|---|---|
| Two or more structures in reaction order (R, any intermediates, P) | Endpoint MEP | [all-endpoint-mep.md](../pdb2reaction-cli/all-endpoint-mep.md) |
| One structure (R) and the distances to drive | Scan (`-s`) | [all-scan-list.md](../pdb2reaction-cli/all-scan-list.md) |
| One TS candidate | TS-only (`--tsopt`) | [all-ts-only.md](../pdb2reaction-cli/all-ts-only.md) |

```bash
# R and P (put intermediates between them, in order)
pdb2reaction all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --tsopt --thermo -o result_mep
# R only: one -s, then one literal per stage
pdb2reaction all -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -s '[("SAM,320,CS1","GPP,321,C7",1.60)]' '[("GPP,321,H11","GLU,186,OE2",0.90)]' \
    --tsopt --thermo -o result_scan
# A TS candidate from another code or an earlier run
pdb2reaction all -i ts_guess.pdb -l 'SAM:1,GPP:-3' --tsopt --thermo -o result_ts
```

- `-c` cuts the active-site cluster around the named residues; without it, the input is used as the cluster. How to choose the residues, radius, and boundary: [pdb2reaction-model-setup](../pdb2reaction-model-setup/SKILL.md).
- `--tsopt` adds TS optimization, IRC, and endpoint optimization; `--thermo` adds frequencies and Gibbs energies for R, TS, and P; `--dft` adds DFT single points on them. `--thermo` and `--dft` need `--tsopt`.
- Each neighbouring pair of inputs becomes one segment. Without intermediates, `--refine-path` splits the MEP where bonds change; `n_segments` can then exceed 1, and each extra segment is a candidate step to check with TS optimization and IRC.
- DFT//MLIP: `--dft` (with `--thermo`) evaluates R, TS, and P with DFT on the MLIP geometries; to run the single points yourself, see [dft.md](../pdb2reaction-cli/dft.md).
- With two or more structures, `-s` is an error. One structure with both `-s` and `--tsopt` runs the scan mode. One structure with neither is an error.

## Run stage by stage and judge each stage

Run the stages as separate commands when you want to check each result before spending GPU time on the next. Pass the same `-q`/`-l`/`-m` and `-b` to every stage, and add `--out-json` so each stage writes the `result.json` read below (standalone commands default to `--no-out-json`). Commands are on the pages under [pdb2reaction-cli](../pdb2reaction-cli/SKILL.md).

| Stage | Command | Role |
|---|---|---|
| Cluster | `extract` | Cuts the active-site cluster, adds cap H, sums the charge |
| MEP | `path-opt` (default in `all`) or `path-search` (`all --refine-path`) | Path between neighbouring structures; its highest-energy image (HEI) is the TS candidate |
| TS | `tsopt` | TS optimization (RS-P-RFO by default, Dimer as an alternative) and n_imag |
| IRC | `irc`, then `opt` on both ends | Follows the reaction mode both ways; the optimized ends become R and P |
| Thermo | `freq` | Frequencies and QRRHO Gibbs energies (`all --thermo`) |
| DFT | `dft` | Single points on R, TS, and P, ωB97M-V/def2-SVP by default (`all --dft`) |

**MEP.** In `path-opt` `result.json`, `optimization_status` must be `"converged"`; `"completed"` only means the run returned. Look at `final_geometries_trj.xyz`, its energy profile, and `hei.pdb`, check that both ends have the same atoms in the same order, and run `bond-summary` on the end pair. For R → IM → P, run one `path-opt` per neighbouring pair. With `path-search`, read `summary.json`, check each segment's `bond_changes`, and start each TS from that segment's `hei_seg_NN.pdb`.

**TS.** In `tsopt` `result.json`, `optimization_status` is `"converged"`, `hessian_status` is `"completed"`, `saddle_validation` is `"first_order"` (`n_imaginary_modes` = 1), and the imaginary mode moves the reacting atoms; the console prints `[tsopt] Converged (n_imag=1).` A successful TS optimization gives one imaginary mode along the reaction coordinate. A run that stops at max cycles without converging computes no Hessian, so it reports no n_imag. A run stopped on an energy plateau (`--stop-plateau`) always computes the Hessian and reports n_imag. `freq` on the TS is optional (all modes, thermochemistry).

**IRC.** `irc` writes `finished_first.xyz` and `finished_last.xyz`. IRC has no pass/fail verdict of its own: `completed` means it returned, and `*_integration_converged` is a diagnostic. Optimize both ends with `opt` and require `optimization_status` `"converged"` for each. Then decide which end is R and which is P by comparing bonds and coordinates with the MEP ends, not from `first`/`last` or from energy. Even if the IRC does not converge, the result is usable when the endpoint optimizations reach the intended R and P. A first-order TS alone does not show the intended reaction.

**Thermo.** The thermochemistry values of R, TS, and P must be finite. n_imag of R and P is a diagnostic, not this gate. Label R and P from the IRC assignment.

**DFT.** Each `result.json` shows `converged: true`. Label the end energies from the IRC assignment, then draw the profile with `energy-diagram` ([outputs.md](outputs.md#energy-diagrams)).

Pitfalls:

- `-l` reads residue names, so it is rejected on bare `.xyz`/`.gjf`. Give a stage a `.pdb` or `.cif` (stages write one when the input had residues), or pass `-q`, or keep `-l` and add `--ref-pdb` with the cluster PDB.
- Standalone `irc` does not write `reactant.pdb`/`product.pdb`; use `all` for the `segments/seg_NN/` layout and automatic R/P orientation.
- In `all`, IRC starts only after the TS converged, its final PHVA (partial Hessian vibrational analysis) finished, and a negative mode was chosen. n_imag = 0, non-convergence, or a failed PHVA stops the segment before IRC and keeps the TS files. A converged TS with n_imag ≥ 2 still runs IRC as a diagnostic (the log says `this is not first-order TS certification`); that IRC is not a TS check.
- After a walltime stop, rerun `all` with the same MEP settings and `--resume-segment N` ([all.md](../pdb2reaction-cli/all.md)), or continue with the stage commands. On any status other than `success`, read `summary.log` and then the stage outputs under `segments/seg_NN/` before retrying.
- A large dense Hessian can exceed GPU memory. Freezing a justified boundary (PHVA) or `--hessian-calc-mode FiniteDifference` lowers the peak, but the Hessian of the moving atoms stays dense.

## What it does

`pdb2reaction` runs MLIP reaction-path calculations on enzyme active-site cluster models. With `-c`, `all` cuts the cluster from a protein–ligand PDB; without it, the PDB/mmCIF/XYZ/GJF model is used as is. It optimizes the endpoints, searches the minimum-energy path (MEP), and stops at the MEP's highest-energy image, a TS candidate; `--tsopt`, `--thermo`, and `--dft` add the later stages. Each stage is also its own subcommand.

- **PDB-native setup**: a residue-aware extractor cuts the cluster, sums residue and ligand charges, and adds cap H at carbon cut points, with no manual atom mapping.
- **Bundled pysisyphus fork**: optimizers, TS search (RS-P-RFO by default, Dimer as an alternative), and IRC have GPU code paths for the MLIP backends; which operations stay on the GPU depends on the backend.
- **Bond-change path splitting**: with `--refine-path`, when R and P differ by more than one step, the path search finds bond changes along the MEP and splits it into narrower segments. These segments are candidates; check each HEI with TS optimization, n_imag, and IRC.

## When to use it, and when not

Use it for:

- Cluster-model enzyme mechanisms, one step or several: `all`.
- Checking a TS candidate with IRC and thermochemistry: TS-only mode, or `tsopt` → `irc` → `freq`.
- DFT//MLIP barriers: DFT on the IRC-refined R, TS, and P (`all --dft`); a TS single point alone is not a barrier.
- A single-point energy on any geometry: `sp` (MLIP energy and forces, `--hess` for a Hessian) or `dft`.

Not for:

- QM methods that PySCF/GPU4PySCF does not provide; use a dedicated QM code (plain DFT runs with `-b dft`).
- Explicit-solvent QM/MM with a force-field environment; pdb2reaction uses cluster models only.
- Free-energy simulations such as umbrella sampling or metadynamics.

## Quick check

```bash
pdb2reaction --version
pdb2reaction --help               # subcommands
pdb2reaction all --help           # main flags
pdb2reaction all --help-advanced  # every flag
```

If `pdb2reaction` is not on PATH, start with [pdb2reaction-install](../pdb2reaction-install/SKILL.md).

## Backend choice

The default is `-b uma`; `orb`, `mace`, `aimnet2`, and `dft` (PySCF/GPU4PySCF) are the alternatives, and [pdb2reaction-install](../pdb2reaction-install/SKILL.md) covers installing and choosing them.

## Where the code lives

Find the installed package with `python -c "import pdb2reaction, os; print(os.path.dirname(pdb2reaction.__file__))"`.

- `pdb2reaction/cli/`: command-line entry point, shared options, `--help-advanced`.
- `pdb2reaction/workflows/`: one module per stage subcommand (`all.py`, `extract.py`, `path_opt.py`, `path_search.py`, `scan.py`, `tsopt.py`, `irc.py`, `freq.py`, `dft.py`, …).
- `pdb2reaction/domain/`: chemistry helpers (bond changes, bond summary, element repair, residue tables).
- `pdb2reaction/backends/`: calculator adapters for UMA, ORB, MACE, AIMNet2, and PySCF DFT.
- `pdb2reaction/io/`: summary writer, energy diagrams, trajectory plots, charge, Hessian cache, altloc fix.
- `pdb2reaction/core/`: `defaults.py` (most calculation defaults; some stay in the command modules, so check live `--help`) and shared utilities.
- `pdb2reaction/mcp/`: MCP server.
- `pysisyphus/` (optimizers, TS search, IRC) and `thermoanalysis/` (thermochemistry): bundled forks, installed as separate top-level packages. A numerical change there needs a regression test and the matching benchmark.

Imports run `cli` → `workflows` → `domain`/`backends`/`io` → `core`; the import-graph check keeps `core` and `domain` from importing `workflows`. More: [docs/architecture.md](../../docs/architecture.md), [CONTRIBUTING.md](../../CONTRIBUTING.md), [check_engineering_markers.py](../../.github/scripts/check_engineering_markers.py), [check_import_graph.py](../../.github/scripts/check_import_graph.py).

## Where to go next

- [ts-strategy.md](ts-strategy.md): studying a mechanism (hypothesis, precision, TS candidates, splitting the reaction, wrong n_imag, a TS that does not come out, comparisons, barriers).
- [outputs.md](outputs.md): `summary.json` and the output tree.
- [pdb2reaction-cli](../pdb2reaction-cli/SKILL.md): running and judging each subcommand.
- [pdb2reaction-model-setup](../pdb2reaction-model-setup/SKILL.md): file formats, residue and atom selectors, charge and multiplicity, and building, trimming, and enlarging the cluster.
- [pdb2reaction-install](../pdb2reaction-install/SKILL.md): the package, backends, CUDA, and checking the environment.
- [pdb2reaction-hpc](../pdb2reaction-hpc/SKILL.md): job scripts.
- [pdb2reaction-mcp](../pdb2reaction-mcp/SKILL.md): MCP tools.
- [colab-local-gpu-runtime](../colab-local-gpu-runtime/SKILL.md): a Colab local runtime.
