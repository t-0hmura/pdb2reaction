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

---

## Machine Learning & Calculators

| Term | Full Name | Description |
|------|-----------|-------------|
| **MLIP** | Machine Learning Interatomic Potential | A model (often neural-network-based) that predicts energies and forces from atomic structures, trained on quantum-mechanical data. |
| **UMA** | Universal Machine-learning potential for Atoms | Meta's family of pretrained MLIPs used as the default calculator backend in pdb2reaction. |

---

## Quantum Chemistry

| Term | Full Name | Description |
|------|-----------|-------------|
| **QM** | Quantum Mechanics | First-principles electronic structure calculations (DFT, HF, post-HF, etc.). |
| **DFT** | Density Functional Theory | A quantum-mechanical method that models electronic structure via electron density functionals. |
| **Hessian** | — | The matrix of second derivatives of energy with respect to atomic coordinates; used for vibrational analysis and TS optimization. |

---

## Structural Biology

| Term | Full Name | Description |
|------|-----------|-------------|
| **PDB** | Protein Data Bank | A file format and database for macromolecular 3D structures. |
| **XYZ** | — | A simple text format listing atomic symbols and Cartesian coordinates. |
| **GJF** | Gaussian Job File | An input format for Gaussian; pdb2reaction reads charge/multiplicity and coordinates from these files. |

---

## Units & Constants

| Term | Description |
|------|-------------|
| **Hartree** | Atomic unit of energy; 1 Hartree ≈ 627.5 kcal/mol ≈ 27.21 eV. |
| **kcal/mol** | Kilocalories per mole; a common unit for reaction energetics. |
| **Bohr** | Atomic unit of length; 1 Bohr ≈ 0.529 Å. |
| **Angstrom (Å)** | 10⁻¹⁰ meters; standard unit for atomic distances. |

---

## See Also

- [Getting Started](getting-started.md) — installation and a first run
- [Concepts & Workflow](concepts.md) — how pocket extraction, MEP search, and post-processing fit together
- [Troubleshooting](troubleshooting.md) — common errors and fixes
- [YAML Reference](yaml-reference.md) — configuration file format
- [UMA Calculator](uma_pysis.md) — machine learning potential details
