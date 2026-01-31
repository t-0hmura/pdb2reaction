# Glossary

This page provides definitions for abbreviations and technical terms used throughout the pdb2reaction documentation.

---

## Reaction Path & Optimization

| Term | Full Name | Description |
|------|-----------|-------------|
| **MEP** | Minimum Energy Path | The lowest-energy pathway connecting reactants to products through a transition state (on a potential energy surface). |
| **TS** | Transition State | A first-order saddle point on the potential energy surface, typically the highest-energy point along the reaction coordinate. |
| **IRC** | Intrinsic Reaction Coordinate | A mass-weighted steepest-descent path from a TS toward reactants and products. Often used to validate TS connectivity. |
| **GSM** | Growing String Method | A string-based method that grows images from endpoints and optimizes them to approximate an MEP. |
| **DMF** | Direct Max Flux | A chain-of-states method for optimizing an MEP by maximizing flux along the pathway. In pdb2reaction it is selected with `--mep-mode dmf`. |
| **NEB** | Nudged Elastic Band | A chain-of-states method that uses spring forces to maintain image spacing along a reaction path. |
| **HEI** | Highest-Energy Image | The image along an MEP with maximum energy; often used as a TS guess. |
| **Image** | — | A single geometry (one "node") along a chain-of-states path. |
| **Segment** | — | An MEP between two adjacent endpoints (e.g., R → I1, I1 → I2, …). |

---

## Optimization Algorithms

| Term | Full Name | Description |
|------|-----------|-------------|
| **L-BFGS** | Limited-memory BFGS | A quasi-Newton optimization algorithm that approximates the Hessian using a limited history of gradients. Used in `--opt-mode light`. |
| **RFO** | Rational Function Optimization | A trust-region optimization method that uses explicit Hessian information. Used in `--opt-mode heavy`. |
| **RS-I-RFO** | Restricted-Step Image-RFO | A variant of RFO for saddle point (TS) optimization that follows one negative eigenvalue. |
| **Dimer** | Dimer Method | A TS optimization method that estimates the lowest curvature mode without computing the full Hessian. Used in `--opt-mode light` for TSOPT. |

---

## Machine Learning & Calculators

| Term | Full Name | Description |
|------|-----------|-------------|
| **MLIP** | Machine Learning Interatomic Potential | A model (often neural-network-based) that predicts energies and forces from atomic structures, trained on quantum-mechanical data. |
| **UMA** | Universal Machine-learning potential for Atoms | Meta's family of pretrained MLIPs used as the default calculator backend in pdb2reaction. |
| **Analytical Hessian** | — | Computing the exact second derivatives of energy; faster but requires more VRAM. |
| **Finite Difference** | — | Approximating derivatives by small displacements; slower but more memory-efficient. |

---

## Quantum Chemistry

| Term | Full Name | Description |
|------|-----------|-------------|
| **QM** | Quantum Mechanics | First-principles electronic structure calculations (DFT, HF, post-HF, etc.). |
| **DFT** | Density Functional Theory | A quantum-mechanical method that models electronic structure via electron density functionals. |
| **Hessian** | — | The matrix of second derivatives of energy with respect to atomic coordinates; used for vibrational analysis and TS optimization. |
| **SP** | Single Point | A calculation at a fixed geometry (no optimization); often used for higher-level energy refinement. |
| **Spin Multiplicity** | — | 2S+1, where S is total spin. Singlet = 1, doublet = 2, triplet = 3, etc. |

---

## Structural Biology & Pocket Extraction

| Term | Full Name | Description |
|------|-----------|-------------|
| **PDB** | Protein Data Bank | A file format and database for macromolecular 3D structures. |
| **XYZ** | — | A simple text format listing atomic symbols and Cartesian coordinates. |
| **GJF** | Gaussian Job File | An input format for Gaussian; pdb2reaction reads charge/multiplicity and coordinates from these files. |
| **Pocket** | Active-site Pocket | A truncated structure around the substrate(s) used to reduce system size for MEP/TS search. Also called "cluster model". |
| **Cluster Model** | — | Synonym for pocket; a computationally tractable subset of the full enzyme–substrate complex. |
| **Link Hydrogen** | — | A hydrogen atom added to cap severed bonds when extracting a pocket from a larger structure. |
| **Backbone** | — | The main chain of a protein (N–Cα–C–O atoms). Can be excluded during pocket extraction with `--exclude-backbone True`. |

---

## Thermochemistry

| Term | Full Name | Description |
|------|-----------|-------------|
| **ZPE** | Zero-Point Energy | The vibrational energy at 0 K; a quantum correction to the electronic energy. |
| **Gibbs Energy** | Free Energy (G) | G = H − TS; includes thermal and entropic contributions. |
| **Enthalpy** | (H) | H = E + PV; total heat content at constant pressure. |
| **Entropy** | (S) | A measure of disorder; contributes −TS to Gibbs energy. |

---

## Units & Constants

| Term | Description |
|------|-------------|
| **Hartree** | Atomic unit of energy; 1 Hartree ≈ 627.5 kcal/mol ≈ 27.21 eV. |
| **kcal/mol** | Kilocalories per mole; a common unit for reaction energetics. |
| **kJ/mol** | Kilojoules per mole; 1 kcal/mol ≈ 4.184 kJ/mol. |
| **eV** | Electron volt; 1 eV ≈ 23.06 kcal/mol. |
| **Bohr** | Atomic unit of length; 1 Bohr ≈ 0.529 Å. |
| **Angstrom (Å)** | 10⁻¹⁰ meters; standard unit for atomic distances. |

---

## CLI Conventions

| Term | Description |
|------|-------------|
| **Boolean option** | CLI flags that take `True` or `False` (capitalized). Example: `--tsopt True`. |
| **Residue selector** | A specification like `'SAM,GPP'` (names) or `'A:123,B:456'` (chain:ID). |
| **Atom selector** | A specification like `'TYR,285,CA'` identifying a specific atom by residue name, number, and atom name. |

---

## See Also

- [Getting Started](getting-started.md) — installation and a first run
- [Concepts & Workflow](concepts.md) — how pocket extraction, MEP search, and post-processing fit together
- [Troubleshooting](troubleshooting.md) — common errors and fixes
- [YAML Reference](yaml-reference.md) — configuration file format
- [UMA Calculator](uma_pysis.md) — machine learning potential details
