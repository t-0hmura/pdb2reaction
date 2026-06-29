# Glossary

## Reaction Path & Optimization

| Term | Full Name | Description |
|------|-----------|-------------|
| **MEP** | Minimum Energy Path | The lowest-energy pathway on a potential energy surface (PES) connecting reactants to products through a transition state. |
| **TS** | Transition State | A first-order saddle point on the potential energy surface — a stationary point with exactly one direction of negative curvature (one imaginary frequency) along the reaction coordinate. |
| **IRC** | Intrinsic Reaction Coordinate | Classically defined as the mass-weighted steepest-descent path from a TS toward reactants and products, used to validate TS connectivity. In pdb2reaction the EulerPC integrator runs in **unweighted Cartesian coordinates** (`irc.md`: `--step-size` is in Bohr of unweighted Cartesian); `geom.coord_type` is forced to `cart`. |
| **GSM** | Growing String Method | A string-based method that grows images from endpoints and optimizes them to approximate an MEP. |
| **DMF** | Direct Max Flux | A chain-of-states method for optimizing an MEP by maximizing flux along the pathway. In pdb2reaction it is selected with `--mep-mode dmf`. |
| **HEI** | Highest-Energy Image | The image along an MEP with maximum energy; often used as a TS guess. |
| **Image** | — | A single geometry (one "node") along a chain-of-states path. |
| **Segment** | — | An MEP between two adjacent endpoints (e.g., R → I1, I1 → I2, …). |
| **Reactive segment** | — | A segment in which covalent bond changes are detected between the endpoints. Only reactive segments proceed to TS optimization. |
| **Bridge segment** | — | A segment connecting two non-adjacent intermediates that still contains unresolved bond changes; `path-search` recursively subdivides bridge segments until all reactive regions are isolated. |
| **Kink** | — | A region along an MEP where no covalent bond change is detected but a geometric distortion persists. `path-search` inserts linearly interpolated nodes and optimizes them individually rather than running a full string calculation. |
| **PES** | Potential Energy Surface | A hypersurface of energy as a function of atomic coordinates. |

## Optimization Algorithms

| Term | Full Name | Description |
|------|-----------|-------------|
| **BFGS** | Broyden-Fletcher-Goldfarb-Shanno | A quasi-Newton Hessian update scheme (`hessian_update: bfgs`). |
| **L-BFGS** | Limited-memory BFGS | A quasi-Newton optimization algorithm that approximates the Hessian using a limited history of gradients. Used in `opt --opt-mode grad`. |
| **RFO** | Rational Function Optimization | A trust-region optimization method that uses explicit Hessian information. Used in `opt --opt-mode hess`. |
| **RS-I-RFO** | Restricted-Step Image-RFO | A variant of RFO for saddle point (TS) optimization that follows one negative eigenvalue. Default for `tsopt --opt-mode hess`. |
| **Dimer** | Dimer Method | A TS optimization method that estimates the lowest curvature mode without computing the full Hessian. Used in `tsopt --opt-mode grad`. pdb2reaction uses a Hessian-Guided Dimer variant that periodically evaluates an exact active-subspace Hessian to update the dimer direction. |
| **Bofill** | Bofill Update | A Hessian update scheme that blends SR1 (symmetric rank-one) and PSB (Powell-symmetric-Broyden) updates, well suited to saddle-point searches. Selected via `hessian_update: bofill` and used by RS-I-RFO and the Dimer flatten loop. |
| **SR1** | Symmetric Rank-One | A rank-one Hessian update scheme; one of the two components blended by Bofill. |
| **PSB** | Powell-Symmetric-Broyden | A symmetric Hessian update scheme; the second component blended by Bofill. |
| **EulerPC** | Euler Predictor-Corrector | An integration scheme for IRC calculations: a predictor step along the gradient direction followed by a corrector step that refines the path. |
| **PHVA** | Partial Hessian Vibrational Analysis | Vibrational analysis performed only on the active (non-frozen) degrees of freedom. Automatically applied when `freeze_atoms` is set. |
| **Active DOF** | Active Degrees of Freedom | The 3N Cartesian coordinates of atoms not listed in `freeze_atoms`. PHVA, partial-Hessian TS optimization, and the analytical Hessian path operate only on this active subspace; frozen atoms contribute neither rows nor columns to the reduced Hessian. |
| **DLC** | Delocalized Internal Coordinates | A redundant internal coordinate system constructed from interatomic distances, angles, and dihedrals. Available via `coord_type: dlc` (default is `cart` = Cartesian). |

## Machine Learning & Calculators

| Term | Full Name | Description |
|------|-----------|-------------|
| **MLIP** | Machine Learning Interatomic Potential | A model (often neural-network-based) that predicts energies and forces from atomic structures, trained on quantum-mechanical data. |
| **UMA** | Universal Models for Atoms | Meta's family of pretrained MLIPs used as the default calculator backend in pdb2reaction. |
| **ORB** | ORB Models | Orbital Materials' MLIP backend. Selected with `-b orb`. |
| **MACE** | MACE | Equivariant message-passing MLIP. Selected with `-b mace`. |
| **AIMNet2** | Atoms-in-Molecules Neural Network Potential, 2nd generation | Charge-aware neural-network potential (Anstine et al., *Chem. Sci.* 2025); selected with `-b aimnet2`. |
| **xTB** | Extended Tight Binding | A semi-empirical quantum chemistry method. In pdb2reaction, used for implicit solvation correction via `--solvent`. |
| **fairchem** | — | Meta's open-source foundation-model toolkit that ships the UMA family of checkpoints. pdb2reaction depends on `fairchem-core` to load UMA predictors. |
| **ASE** | Atomic Simulation Environment | Python framework providing the Calculator API used by all MLIP backends in pdb2reaction (Larsen et al., *J. Phys. Condens. Matter* 2017). |
| **task_name** | — | UMA task tag recorded in each inference batch (YAML: `calc.task_name`, default `omol`). Selects the UMA task/preset that a checkpoint was trained for. |
| **Analytical Hessian** | — | Exact evaluation of the Hessian matrix via automatic differentiation; faster than finite differences but requires more VRAM. Selected with `--hessian-calc-mode Analytical`. |
| **Finite Difference** | — | Approximation of the Hessian via finite nuclear displacements; slower but more memory-efficient. Selected with `--hessian-calc-mode FiniteDifference` (default). |

## Quantum Chemistry

| Term | Full Name | Description |
|------|-----------|-------------|
| **QM** | Quantum Mechanics | First-principles electronic structure calculations (DFT, HF, post-HF, etc.). |
| **DFT** | Density Functional Theory | A quantum-mechanical method that models electronic structure via electron density functionals. |
| **DFT//MLIP** | — | Composite-method notation: DFT single-point energies evaluated at MLIP-optimized geometries. Combines MLIP geometry/dynamics with DFT energy refinement. The `//` separator follows the standard quantum-chemistry convention "energy-level // geometry-level". |
| **Hessian** | — | The matrix of second derivatives of energy with respect to atomic coordinates. Eigenvalues yield vibrational frequencies; eigenvectors yield vibrational modes (displacement vectors). Used for vibrational analysis and TS optimization. |
| **SP** | Single Point | A calculation at a fixed geometry (no optimization); often used for higher-level energy refinement. |
| **Spin Multiplicity** | — | 2S+1, where S is total spin. Singlet = 1, doublet = 2, triplet = 3, etc. Specified with `-m/--multiplicity` (default: 1). |
| **ALPB** | Analytical Linearized Poisson-Boltzmann | An implicit-solvent model available via xTB (`--solvent-model alpb`, default). |
| **CPCM-X** | Extended Conductor-like Polarizable Continuum Solvation Model | An implicit-solvent model available via xTB (`--solvent-model cpcmx`); "X" stands for "eXtended", obtained by coupling CPCM with COSMO-RS-style σ-profiles and an SMD-style non-electrostatic term to give realistic solvation free energies for arbitrary solvents (Stahn, Ehlert, Grimme, *J. Phys. Chem. A* 2023). |
| **cyipopt** | — | Python bindings for the IPOPT interior-point optimizer. Required by the DMF (`--mep-mode dmf`) path refinement pipeline. |
| **IPOPT** | Interior Point OPTimizer | Open-source nonlinear constrained optimizer (Wächter & Biegler 2006) used by the DMF path-refinement solver via `cyipopt` bindings. |
| **SCF** | Self-Consistent Field | Iterative procedure that converges the electronic wavefunction in DFT/HF; controlled in `pdb2reaction dft` by `--max-cycle` and `--conv-tol`. |

## Structural Biology & Active Site Model Extraction

| Term | Full Name | Description |
|------|-----------|-------------|
| **PDB** | Protein Data Bank | A file format and database for macromolecular 3D structures. |
| **XYZ** | — | A simple text format listing atomic symbols and Cartesian coordinates. |
| **GJF** | Gaussian Job File | An input format for Gaussian; pdb2reaction reads charge/multiplicity and coordinates from these files. |
| **Active Site Model** | Active Site Model (Binding Pocket) | The extraction region around the substrate(s), defined by `-c/--center` and `-r/--radius`. pdb2reaction uses "active site model" and "cluster model" interchangeably when describing the downstream calculation input; strictly, "active site model" refers to the geometric selection and "cluster model" to the same selection after it has been capped with cap hydrogens and prepared for QM/MLIP calculation. |
| **Cluster Model** | — | The QM/MLIP computational subsystem obtained by taking the extracted active site model and capping severed covalent bonds with hydrogen atoms (cap hydrogens). |
| **Cap Hydrogen** | — | A hydrogen atom added to cap severed bonds when extracting an active site model from a larger structure. |
| **Backbone** | — | The main chain of a protein (N–Cα–C–O atoms). Can be excluded during active site model extraction with `--exclude-backbone`. |

## Thermochemistry

| Term | Full Name | Description |
|------|-----------|-------------|
| **ZPE** | Zero-Point Energy | The vibrational energy at 0 K; a quantum correction to the electronic energy. |
| **Gibbs Energy** | Gibbs Free Energy (G) | G = H − TS; includes thermal and entropic contributions. |
| **Enthalpy** | (H) | H = E + PV; total heat content at constant pressure. |
| **Entropy** | (S) | A measure of disorder; contributes −TS to Gibbs energy. |
| **QRRHO** | Quasi-Rigid-Rotor Harmonic Oscillator | A thermochemical approximation incorporating Grimme's correction for low-frequency vibrations. Automatically applied in `freq`. |

## Units & Constants

| Term | Description |
|------|-------------|
| **Hartree** | Atomic unit of energy; 1 Hartree ≈ 627.5 kcal/mol ≈ 27.21 eV. |
| **RMSD** | Root-Mean-Square Deviation; used as the segment stitch / bridge similarity metric in `path-search` (`stitch_rmsd_thresh`, `bridge_rmsd_thresh`). |
| **MAE** | Mean Absolute Error; used in benchmark / regression reports. |
| **CI** | Confidence Interval (statistical context). Distinct from quantum-chemistry "configuration interaction"; pdb2reaction uses CI only in benchmark statistics. |
| **kcal/mol** | Kilocalories per mole; a common unit for reaction energetics. |
| **kJ/mol** | Kilojoules per mole; 1 kcal/mol ≈ 4.184 kJ/mol. |
| **eV** | Electron volt; 1 eV ≈ 23.06 kcal/mol. |
| **Bohr** | Atomic unit of length; 1 Bohr ≈ 0.529 Å. |
| **Angstrom (Å)** | 10⁻¹⁰ m; standard unit for interatomic distances. |
| **cm⁻¹** | Reciprocal centimeters (wavenumber); the standard unit for vibrational frequencies. Imaginary frequencies appear as negative values. |
| **Imaginary Frequency** | A vibrational frequency corresponding to a negative eigenvalue of the Hessian. A TS has exactly one (first-order saddle point). Reported as a negative cm⁻¹ value. The detection threshold is `hessian_dimer.neg_freq_thresh_cm` (default 5 cm⁻¹). |

(frequency-thresholds)=
### Frequency thresholds: 5 cm⁻¹ imaginary detection vs 100 cm⁻¹ QRRHO rotor cutoff

Two unrelated cm⁻¹ thresholds appear in `pdb2reaction`. They act on different mode populations and serve different purposes:

| Threshold | Role | Source |
|-----------|------|--------|
| **5 cm⁻¹** | *Imaginary-mode detection cutoff.* Negative eigenvalues with magnitude below this are not counted as imaginary (treated as rigid-body / numerical noise). | `pdb2reaction/core/defaults.py` as `hessian_dimer.neg_freq_thresh_cm = 5.0`; tunable via YAML. |
| **100 cm⁻¹** | *QRRHO rotor cutoff* (Grimme). Positive low-frequency vibrations are damped between harmonic-oscillator and free-rotor entropy in `freq` thermochemistry; it changes only entropy / Gibbs free energy. | `thermoanalysis/config.py` as `ROTOR_CUT_DEFAULT = 100.0`. |

## CLI Conventions

| Term | Description |
|------|-------------|
| **Boolean option** | CLI flags that accept toggle form (`--flag` / `--no-flag`) or value form (`True`/`False`, `yes`/`no`, `1`/`0`). Example: `--tsopt`. |
| **Residue selector** | A specification like `'SAM,GPP'` (names) or `'A:123,B:456'` (chain:ID). |
| **Atom selector** | A specification like `'TYR,285,CA'` identifying a specific atom by residue name, number, and atom name. |

## See Also

- [Getting Started](getting-started.md) — installation and a first run
- [Installation](installation.md) — setup and dependencies
- [all](all.md) — how active site model extraction, MEP search, and post-processing fit together
- [Common Error Recipes](recipes-common-errors.md) — symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — common errors and fixes
- [YAML Reference](yaml-reference.md) — configuration file format
- [MLIP Calculator](uma-pysis.md) — machine learning potential details
