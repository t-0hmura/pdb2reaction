# Glossary

Each abbreviation, method name, and unit used in the docs is defined here in one line, grouped by field. Options, output fields, and status values are on each command page and in [JSON Output](json-output.md).

## Reaction Path & Optimization

| Term | Full Name | Description |
|------|-----------|-------------|
| **MEP** | Minimum Energy Path | The lowest-energy pathway on a potential energy surface (PES) connecting reactants to products through a transition state. |
| **TS** | Transition State | A first-order saddle point on the potential energy surface — a stationary point with exactly one direction of negative curvature (one imaginary frequency) along the reaction coordinate. |
| **n_imag** | Number of imaginary modes | The number of vibrational modes below the imaginary-mode criterion (ν < −5.00 cm⁻¹ by default); a TS has n_imag = 1, and `result.json` records it as `n_imaginary_modes` (`tsopt`) or `n_imaginary` (`freq`). |
| **IRC** | Intrinsic Reaction Coordinate | Classically defined as the mass-weighted steepest-descent path from a TS toward reactants and products, used to validate TS connectivity. In pdb2reaction the EulerPC integrator advances mass-weighted coordinates; `--step-size` is a length in Bohr in unweighted Cartesian coordinates. |
| **GSM** | Growing String Method | A string-based method that grows images from endpoints and optimizes them to approximate an MEP. |
| **DMF** | Direct Max Flux | A chain-of-states method for optimizing an MEP by maximizing flux along the pathway. In pdb2reaction it is selected with `--mep-mode dmf`. |
| **FB-ENM** | Flat-Bottom Elastic Network Model | The model that builds the initial path for DMF (Koda & Saito, *J. Chem. Theory Comput.* 2024); CFB-ENM is its correlated variant. |
| **HEI** | Highest-Energy Image | The image along an MEP with maximum energy; often used as a TS guess. HEI±1 are the images on either side of it. |
| **Image** | — | A single geometry (one "node") along a chain-of-states path. |
| **COS** | Chain-of-States | A path method that optimizes a whole chain of images together, such as GSM, the string optimizer (`stopt`), and DMF. |
| **Segment** | — | An MEP between two adjacent endpoints (e.g., R → I1, I1 → I2, …). |
| **Reactive segment** | — | A segment that holds a TS candidate: every MEP segment except bridges and kinks (normally one whose ends differ by a covalent bond change), or the input TS in [TS-only mode](quickstart-tsopt.md) (one structure with `all --tsopt`); `all` runs the requested TS optimization, IRC, thermochemistry, and DFT only on these. |
| **Bridge segment** | — | A short connecting MEP between two neighbouring segments: when `path-search` joins the segments into one path and the end of one does not coincide with the start of the next, with no bond change between them, it fills the gap with a bridge ([path-search → How it works](path-search.md#how-it-works), step 5). |
| **Kink** | — | A segment where only the conformation changes: the two structures that `path-search` optimizes on either side of the HEI (End1 and End2; [path-search → How it works](path-search.md#how-it-works), step 2) differ by no covalent bond change. `path-search` fills it with a few linearly interpolated nodes (`search.kink_max_nodes`, default 3) and optimizes each one instead of running a new GSM or DMF path. |
| **PES** | Potential Energy Surface | A hypersurface of energy as a function of atomic coordinates. |

## Optimization Algorithms

| Term | Full Name | Description |
|------|-----------|-------------|
| **BFGS** | Broyden-Fletcher-Goldfarb-Shanno | A quasi-Newton Hessian update scheme (`hessian_update: bfgs`). |
| **TS-BFGS** | Transition-State BFGS | A BFGS-type Hessian update that does not require positive curvature along the step, so the model Hessian can stay indefinite. Default for RFO minimization (`hessian_update: ts_bfgs`). |
| **L-BFGS** | Limited-memory BFGS | A quasi-Newton optimization algorithm that approximates the Hessian using a limited history of gradients. Used in `opt --opt-mode grad`. |
| **RFO** | Rational Function Optimization | A trust-region optimization method that uses explicit Hessian information. Used in `opt --opt-mode hess`. |
| **RS-I-RFO** | Restricted-Step Image-RFO | A variant of RFO for saddle point (TS) optimization that follows one negative eigenvalue. Selectable via `tsopt --opt-mode rsirfo` (the `hess` default is RS-P-RFO). |
| **RS-P-RFO** | Restricted-Step Partitioned RFO | The default TS optimizer of `tsopt` (`--opt-mode hess` or `rsprfo`): an RFO variant that maximizes the energy along the TS mode and minimizes it along all other modes (Banerjee et al. 1985). |
| **TRIM** | Trust-Region Image Minimization | A trust-region TS optimizer that reverses the sign of the Hessian eigenvalue and gradient along the TS mode and then minimizes (Helgaker 1991); selected with `tsopt --opt-mode trim`. |
| **Flatten** | — | With `--flatten` (`opt`, `tsopt`, and `all`; off by default), the structure is displaced along the extra imaginary modes and optimized again, until none is left for `opt` or one is left for `tsopt`, or the round limit is reached. |
| **Dimer** | Dimer Method | A TS optimization method that follows the lowest-curvature mode with gradient steps. pdb2reaction uses a Hessian-guided variant: an exact Hessian of the active atoms sets the dimer direction at the start and refreshes it every `hessian_dimer.update_interval_hessian` steps (default 500). Selected with `tsopt --opt-mode dimer` (`grad` is an alias). |
| **Bofill** | Bofill Update | A Hessian update scheme that blends SR1 (symmetric rank-one) and PSB (Powell-symmetric-Broyden) updates, well suited to saddle-point searches. The default `hessian_update` of the `rsirfo` section (RS-P-RFO / RS-I-RFO / TRIM) and of `irc`. |
| **SR1** | Symmetric Rank-One | A rank-one Hessian update scheme; one of the two components blended by Bofill. |
| **PSB** | Powell-Symmetric-Broyden | A symmetric Hessian update scheme; the second component blended by Bofill. |
| **EulerPC** | Euler Predictor-Corrector | An integration scheme for IRC calculations: a predictor step along the gradient direction followed by a corrector step that refines the path. |
| **PHVA** | Partial Hessian Vibrational Analysis | Vibrational analysis performed only on the active (non-frozen) degrees of freedom. Automatically applied when `freeze_atoms` is set. |
| **Active DOF** | Active Degrees of Freedom | The 3N Cartesian coordinates of atoms not listed in `freeze_atoms`. PHVA, partial-Hessian TS optimization, and the analytical Hessian path operate only on this active subspace; frozen atoms contribute neither rows nor columns to the reduced Hessian. |
| **DLC** | Delocalized Internal Coordinates | A non-redundant set of internal coordinates, each a linear combination of the primitive distances, angles, and dihedrals. Selected with `--coord-type dlc` (YAML `geom.coord_type: dlc`); the default is `cart` (Cartesian). |

## Machine Learning & Calculators

| Term | Full Name | Description |
|------|-----------|-------------|
| **MLIP** | Machine Learning Interatomic Potential | A model (often neural-network-based) that predicts energies and forces from atomic structures, trained on quantum-mechanical data. |
| **UMA** | Universal Models for Atoms | Meta's family of pretrained MLIPs used as the default calculator backend in pdb2reaction. |
| **ORB** | ORB Models | Orbital Materials' MLIP backend. Selected with `-b orb`. |
| **MACE** | MACE | Equivariant message-passing MLIP. Selected with `-b mace`. |
| **AIMNet2** | Atoms-in-Molecules Neural Network Potential, 2nd generation | Charge-aware neural-network potential (Anstine et al., *Chem. Sci.* 2025); selected with `-b aimnet2`. |
| **fairchem** | — | Meta's open-source foundation-model toolkit that ships the UMA family of checkpoints. pdb2reaction depends on `fairchem-core` to load UMA predictors. |
| **ASE** | Atomic Simulation Environment | Python framework providing the Calculator API used by all MLIP backends in pdb2reaction (Larsen et al., *J. Phys. Condens. Matter* 2017). |
| **task_name** | — | UMA task tag recorded in each inference batch (YAML: `calc.task_name`, default `omol`). Selects the UMA task/preset that a checkpoint was trained for. |
| **Analytical Hessian** | — | Automatic differentiation of the selected MLIP energy (up to floating-point/autograd behavior), avoiding finite-displacement truncation error. Runtime and accelerator-memory cost are backend/model/system dependent. Selected with `--hessian-calc-mode Analytical`. |
| **Finite Difference** | — | Approximation of the Hessian from finite nuclear displacements. It is the portable default and usually uses less peak accelerator memory, but runtime and displacement error depend on the setup. Selected with `--hessian-calc-mode FiniteDifference`. |

## Quantum Chemistry

| Term | Full Name | Description |
|------|-----------|-------------|
| **QM** | Quantum Mechanics | First-principles electronic structure calculations (DFT, HF, post-HF, etc.). |
| **DFT** | Density Functional Theory | A quantum-mechanical method that models electronic structure via electron density functionals. |
| **DFT//MLIP** | — | Composite-method notation: DFT single-point energies evaluated at MLIP-optimized geometries. Combines MLIP geometry/dynamics with a higher-level DFT energy evaluation. The `//` separator follows the standard quantum-chemistry convention "energy-level // geometry-level". |
| **MLIP_Gibbs**, **DFT//MLIP_Gibbs** | — | Method labels in `summary.json`, next to `MEP`, `MLIP`, and `DFT`: `MLIP_Gibbs` is the MLIP Gibbs energy, and `DFT//MLIP_Gibbs` is the DFT energy plus the MLIP thermal correction (with `-b dft`, `MLIP` and `MLIP_Gibbs` become `DFT` and `DFT_Gibbs`). |
| **Hessian** | — | The matrix of second derivatives of energy with respect to atomic coordinates. Eigenvalues yield vibrational frequencies; eigenvectors yield vibrational modes (displacement vectors). Used for vibrational analysis and TS optimization. |
| **SP** | Single Point | A calculation at a fixed geometry (no optimization); often used for a higher-level energy evaluation. |
| **Spin Multiplicity** | — | 2S+1, where S is total spin. Singlet = 1, doublet = 2, triplet = 3, etc. Specified with `-m/--multiplicity` (default: 1). |
| **cyipopt** | — | Python bindings for the IPOPT interior-point optimizer. Required by the DMF (`--mep-mode dmf`) path refinement pipeline. |
| **IPOPT** | Interior Point OPTimizer | Open-source nonlinear constrained optimizer (Wächter & Biegler 2006) used by the DMF path-refinement solver via `cyipopt` bindings. |
| **SCF** | Self-Consistent Field | Iterative procedure that converges the electronic wavefunction in DFT/HF; controlled in `pdb2reaction dft` by `--scf-max-cycles` and `--scf-tol`. |

## Structural Biology & Active Site Model Extraction

| Term | Full Name | Description |
|------|-----------|-------------|
| **PDB** | Protein Data Bank | A file format and database for macromolecular 3D structures. |
| **XYZ** | — | A simple text format listing atomic symbols and Cartesian coordinates. |
| **GJF** | Gaussian Job File | An input format for Gaussian; pdb2reaction reads charge/multiplicity and coordinates from these files. |
| **Active Site Model** | Active Site Model (Binding Pocket) | The region around the substrate(s) selected by `-c/--center` and `-r/--radius`. The docs use *active-site model* and *cluster model* for the same thing: the model that `extract` writes, with cap hydrogens on the severed bonds. |
| **Cluster Model** | — | The QM/MLIP computational subsystem obtained by taking the extracted region and capping severed covalent bonds with hydrogen atoms (cap hydrogens). |
| **Cap Hydrogen** | — | A hydrogen atom added to cap severed bonds when extracting an active site model from a larger structure. Also called a link hydrogen (`--add-linkh`, `--freeze-links`, `n_link_hydrogens`). |
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
| **kcal/mol** | Kilocalories per mole; a common unit for reaction energetics. |
| **kJ/mol** | Kilojoules per mole; 1 kcal/mol ≈ 4.184 kJ/mol. |
| **eV** | Electron volt; 1 eV ≈ 23.06 kcal/mol. |
| **Bohr** | Atomic unit of length; 1 Bohr ≈ 0.529 Å. |
| **Angstrom (Å)** | 10⁻¹⁰ m; standard unit for interatomic distances. |
| **cm⁻¹** | Reciprocal centimeters (wavenumber); the standard unit for vibrational frequencies. Imaginary frequencies appear as negative values. |
| **Imaginary Frequency** | A vibrational frequency corresponding to a negative eigenvalue of the Hessian. A TS has exactly one (first-order saddle point). Reported as a negative cm⁻¹ value. |

(frequency-thresholds)=
### Imaginary-mode criterion and QRRHO rotor cutoff

| Threshold | Role | Source |
|-----------|------|--------|
| **ν < −5.00 cm⁻¹** | Default imaginary-mode criterion. | Configurable: `freq.zero_cutoff_cm`. |
| **100 cm⁻¹** | *QRRHO rotor cutoff* (Grimme). Positive low-frequency vibrations are damped between harmonic-oscillator and free-rotor entropy in `freq` thermochemistry; it changes only entropy / Gibbs free energy. | Fixed (not configurable in pdb2reaction). |

## CLI Conventions

Boolean options, residue selectors, and atom selectors are described in [Common options and selectors](cli-conventions.md).

## Notes

* **Imaginary modes and negative signs**: every frequency is reported with its sign; the count of all negative values (`n_negative_modes` in `result.json`) is a separate diagnostic from n_imag and does not change convergence.

## See Also

- [Getting Started](getting-started.md) — the shortest run and which page to read next
- [Installation](installation.md) — setup and dependencies
- [all](all.md) — the full workflow
- [Troubleshooting](troubleshooting.md) — common errors and fixes
- [YAML Reference](yaml-reference.md) — configuration file format
- [MLIP Backends](backends.md) — machine learning potential details
