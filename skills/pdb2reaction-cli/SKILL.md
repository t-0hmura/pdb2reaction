---
name: pdb2reaction-cli
description: "Task-level guide to the 18 pdb2reaction subcommands: which command to run, a minimal working invocation, how to judge success, common pitfalls and recovery, and what to run next. SKILL.md is a one-line input-to-output cheatsheet with shared conventions (charge, spin, backend, output, dry run); all.md and its three mode pages, extract.md (with frozen atoms), opt.md, path.md (path-opt and path-search), scan.md (scan, scan2d, scan3d), tsopt.md, irc.md, freq.md, dft.md, and utilities.md (sp, fix-altloc, add-elem-info, bond-summary, trj2fig, energy-diagram) hold the details. Full flag lists come from --help-advanced and the generated CLI reference. TRIGGER on a question about a specific subcommand, a shell invocation, or a run that failed. SKIP for install, HPC, output-schema, or structure-format questions, and for choosing what goes into the cluster (pdb2reaction-model-setup)."
---

# pdb2reaction CLI

Pick the command from the table, run its minimal command, then read its page
for how to judge success. Charge and multiplicity must be chemically right for
every run: the examples use a neutral singlet `-q 0 -m 1`, so replace those
values with the verified charge (or `-l 'RES:Q,...'` for PDB/mmCIF) and set
`-m` for open-shell systems.

| sub | role | minimal command | primary output | page |
|---|---|---|---|---|
| `all` | Extraction, MEP or staged scan, then TS/IRC, freq, DFT as requested | `pdb2reaction all -i 1.R.pdb 3.P.pdb -q 0 -m 1 --tsopt --thermo -o out` | `out/summary.json`, `out/segments/seg_NN/{reactant,ts,product}.*` | [all.md](all.md), [endpoint MEP](all-endpoint-mep.md), [scan](all-scan-list.md), [TS-only](all-ts-only.md) |
| `extract` | Active-site cluster cut | `pdb2reaction extract -i raw.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' -r 2.6 -o cluster.pdb` | `cluster.pdb` (`-o` is a file path, not a directory) | [extract.md](extract.md) |
| `path-opt` | Single-pass MEP between two structures | `pdb2reaction path-opt -i 1.R.pdb 2.P.pdb -q 0 -m 1 -o out` | `out/final_geometries_trj.xyz`, `out/hei.xyz` | [path.md](path.md) |
| `path-search` | Recursive MEP split at bond changes | `pdb2reaction path-search -i 1.R.pdb 3.P.pdb -q 0 -m 1 -o out` | `out/mep_trj.xyz`, `out/hei_seg_NN.xyz`, `out/summary.json` | [path.md](path.md) |
| `opt` | Geometry minimization (L-BFGS / RFO) | `pdb2reaction opt -i geom.pdb -q 0 -m 1 -o out` | `out/final_geometry.xyz` | [opt.md](opt.md) |
| `tsopt` | TS optimization (RS-P-RFO / Dimer) | `pdb2reaction tsopt -i ts.xyz -q 0 -m 1 -o out` | `out/final_geometry.xyz`; confirm with n_imag = 1, the mode, and IRC | [tsopt.md](tsopt.md) |
| `irc` | IRC from a TS | `pdb2reaction irc -i ts.xyz -q 0 -m 1 -o out` | `out/{forward,backward,finished}_irc_trj.xyz` | [irc.md](irc.md) |
| `freq` | Hessian and QRRHO thermochemistry | `pdb2reaction freq -i geom.xyz -q 0 -m 1 -o out` | `out/frequencies_cm-1.txt`; `thermoanalysis.yaml` with `--dump` | [freq.md](freq.md) |
| `dft` | Single-point DFT (PySCF / GPU4PySCF) | `pdb2reaction dft -i geom.pdb -q 0 -m 1 --func-basis 'wb97m-v/def2-tzvpd' -o out` | `out/result.yaml` | [dft.md](dft.md) |
| `scan` | Staged scan under restraints | `pdb2reaction scan -i 1.R.pdb -q 0 -m 1 -s '[(a,b,1.6)]' -o out` | `out/scan_trj.xyz`, `stage_NN/result.xyz` | [scan.md](scan.md) |
| `scan2d` | 2D distance grid | `pdb2reaction scan2d -i 1.R.pdb -q 0 -m 1 -s '[(a,b,1.3,3.1),(c,d,1.2,3.2)]' -o out` | `out/surface.csv`, `out/scan2d_map.png` | [scan.md](scan.md) |
| `scan3d` | 3D distance grid | `pdb2reaction scan3d -i 1.R.pdb -q 0 -m 1 -s '[(a,b,L,H),(c,d,L,H),(e,f,L,H)]' -o out` | `out/surface.csv`, `out/scan3d_density.html` | [scan.md](scan.md) |
| `sp` | Single-point energy and forces | `pdb2reaction sp -i geom.pdb -q 0 -m 1 -o out` | energy on stdout, `out/forces.npy`; `hessian.npy` with `--hess` | [utilities.md](utilities.md) |
| `trj2fig` | Energy profile from an XYZ trajectory | `pdb2reaction trj2fig -i trj.xyz` | `energy.png` | [utilities.md](utilities.md) |
| `energy-diagram` | Diagram from energy values | `pdb2reaction energy-diagram -i "[0.0, 21.5, -0.7]" --label-x "['R','TS','P']"` | `energy_diagram.png` | [utilities.md](utilities.md) |
| `add-elem-info` | Fill the PDB element column | `pdb2reaction add-elem-info -i raw.pdb -o fixed.pdb` | `fixed.pdb` | [utilities.md](utilities.md) |
| `fix-altloc` | Keep one alternate location per residue | `pdb2reaction fix-altloc -i raw.pdb -o fixed.pdb` | `fixed.pdb` | [utilities.md](utilities.md) |
| `bond-summary` | Bond changes between consecutive structures | `pdb2reaction bond-summary -i reactant.pdb -i product.pdb` | text on stdout; JSON with `--json` | [utilities.md](utilities.md) |

mmCIF input and very large PDB input are handled as `.pdb` inside the run, so
the `.pdb` names above still apply; the outputs also include a `.cif` that
keeps the original chain IDs and residue numbers.

## Common conventions

- `-i, --input`: input structure(s). Geometry commands read `.pdb`, `.cif` / `.mmcif`, `.xyz`, and `.gjf`.
- `-q, --charge`: total charge. `-l, --ligand-charge 'RES1:Q1,RES2:Q2'`: charges of non-standard residues in a PDB/mmCIF, from which the total is derived. `-m, --multiplicity`: 2S+1, default 1.
- Charge order: explicit `-q` wins over the `-l` derivation, then the YAML value. This holds for `all -c/--center` too: extraction derives the cluster charge, an explicit `-q` sets the total, and a mismatch is reported as a warning. A bare XYZ needs `-q` (or `--ref-pdb` with `-l`); a GJF supplies its header value.
- `-b, --backend`: `uma` (default), `orb`, `mace`, `aimnet2`, or `dft`.
- `-o, --out-dir`: output directory; each command has its own default (`all` writes to `./result_all/`). `extract`, `add-elem-info`, `fix-altloc`, `trj2fig`, and `energy-diagram` take output file paths with `-o` instead.
- `--config FILE`: YAML applied on top of the built-in defaults and below explicit CLI flags.
- `--show-config`: prints the configuration and **continues** with the run. `all`, `path-search`, and `sp` print the settings after merging defaults, YAML, and CLI; the other commands that accept it print the loaded YAML and its top-level keys.
- `--dry-run`: checks options and inputs, then exits before any MLIP or DFT stage. `all -c/--center --dry-run` also runs extraction to check the derived charge and electron parity.
- `--ref-pdb FILE`: a PDB/mmCIF that gives residue names and topology to XYZ/GJF inputs while their coordinates are kept.
- `--solvent NAME`: for MLIP backends, an expensive xTB correction `E_xTB(solvent) - E_xTB(vacuum)`, meant mainly for small molecules in solution; with `-b dft` and the `dft` command, PySCF PCM/SMD (`--solvent-model`). `none` turns it off. Install notes: [xTB solvent correction](../pdb2reaction-install-backends/backends.md#xtb-solvent-correction).

## Frozen atoms

Cap hydrogens, `--freeze-links`, `--freeze-atoms`, and YAML `geom.freeze_atoms` are combined into one frozen set; choose a chemically justified boundary and inspect it rather than reusing one freeze set for every cluster. See [extract.md](extract.md#freeze-atoms-at-the-cluster-boundary).

## Cross-cutting pitfalls

- **`--scan-lists` syntax error.** Each value is a Python literal. Wrap it in single quotes and use double quotes inside; do not confuse a backtick (`` ` ``) with a backslash (`\`).
- **Wrong charge with no error.** Settle protonation and oxidation states and read the per-residue charge breakdown. `--dry-run` prints the charge and checks electron parity, but it cannot prove the chemistry is right.
- **Backend left to the default.** Pass `-b` explicitly in production scripts so a change of default cannot reroute the run.
- **YAML value seems ignored.** The order is built-in defaults, then `--config`, then explicit CLI flags; a value also given on the command line overrides the YAML.
- **Flag not in `--help`.** Advanced flags are listed only by `--help-advanced`; pin the package version when a workflow is shared.
- **Out of memory in the Hessian step.** Keep the default `--hessian-calc-mode FiniteDifference` rather than `Analytical`, which builds an autograd graph; use a justified frozen boundary (PHVA) or a smaller model. The Hessian of the active atoms stays dense, so measure memory on a representative system.
- **UMA with `--uma-workers` above 1 and an explicit `Analytical` Hessian.** This raises `BackendError`; the requested method is never changed silently. Use `--uma-workers 1` for `Analytical`, or keep `FiniteDifference`. ORB, MACE, and AIMNet2 ignore the worker flags, and all four built-in backends implement analytical Hessians.

## Where flags and defaults live

- `pdb2reaction <sub> --help-advanced` lists every flag with its default; `docs/reference/commands/` holds the same text.
- `pdb2reaction.core.defaults` holds the shared defaults as `*_KW` dictionaries and `OUT_DIR_*` paths; some command-only defaults live in the command module, so check `--help-advanced` too:

```bash
python -c "import pdb2reaction.core.defaults as d; print(sorted(n for n in dir(d) if not n.startswith('_')))"
python -c "import pdb2reaction.core.defaults as d; print(d.RSIRFO_KW)"   # or LBFGS_KW, IRC_KW, UMA_CALC_KW
```

## Next step

- [pdb2reaction-overview](../pdb2reaction-overview/SKILL.md): which `all` mode to use, and how to run and judge the stages one by one.
- [Reading outputs](../pdb2reaction-overview/outputs.md): `summary.json` keys, R/TS/P paths, bond changes, energy diagrams.
- [TS strategy](../pdb2reaction-overview/ts-strategy.md): wrong n_imag, or no TS.
- [pdb2reaction-structure-io](../pdb2reaction-structure-io/SKILL.md): input formats, charge and multiplicity.
- [pdb2reaction-model-setup](../pdb2reaction-model-setup/SKILL.md): what goes into the cluster.
- [pdb2reaction-install-backends](../pdb2reaction-install-backends/SKILL.md): install and backends.
- [pdb2reaction-hpc](../pdb2reaction-hpc/SKILL.md): running on PBS or SLURM.
