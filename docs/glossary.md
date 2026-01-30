# Glossary

This page provides definitions for abbreviations and technical terms used throughout the pdb2reaction documentation.

---

## Reaction Path & Optimization

| Term | Full Name | Description |
|------|-----------|-------------|
| **MEP** | Minimum Energy Path | The lowest-energy pathway connecting reactants to products through a transition state. |
| **TS** | Transition State | A saddle point on the potential energy surface representing the highest-energy point along the reaction coordinate. |
| **IRC** | Intrinsic Reaction Coordinate | A mass-weighted steepest-descent path from the TS toward reactants and products. |
| **GSM** | Growing String Method | An algorithm for finding MEPs by iteratively growing a string of images from both ends toward the middle. |
| **NEB** | Nudged Elastic Band | A chain-of-states method that uses spring forces to maintain image spacing along a reaction path. |

---

## Machine Learning & Calculators

| Term | Full Name | Description |
|------|-----------|-------------|
| **MLIP** | Machine Learning Interatomic Potential | A neural-network-based energy/force predictor trained on quantum-mechanical data. |
| **UMA** | Universal Machine-learning potential for Atoms | Meta's family of pre-trained MLIPs used as the default calculator in pdb2reaction. |
| **DMF** | Distance Matrix Fingerprints | Structural descriptors based on interatomic distances, used for similarity checks in GSM. |

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
| **GJF** | Gaussian Job File | Input format for Gaussian; pdb2reaction reads charge/multiplicity and coordinates from these files. |

---

## Units & Constants

| Term | Description |
|------|-------------|
| **Hartree** | Atomic unit of energy; 1 Hartree ≈ 627.5 kcal/mol ≈ 27.21 eV. |
| **kcal/mol** | Kilocalories per mole; common unit for reaction energetics. |
| **Bohr** | Atomic unit of length; 1 Bohr ≈ 0.529 Å. |
| **Angstrom (Å)** | 10⁻¹⁰ meters; standard unit for atomic distances. |

---

## See Also

- [Getting Started](getting-started.md) — Installation and first run
- [YAML Reference](yaml-reference.md) — Configuration file format
- [UMA Calculator](uma_pysis.md) — Machine learning potential details
