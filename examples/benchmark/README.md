# pdb2reaction benchmark inputs (23 reaction steps from 6 enzymes)

This directory contains the raw inputs used for the validation campaign in
the `pdb2reaction` paper (*Automated Reaction-Path Modeling from PDB
Structures Using Machine-Learning Interatomic Potentials*, in preparation).

## What is here

```
examples/benchmark/
  1AH7/    Phospholipase C  (Zn3)
  1O9S/    Histone Lys methyltransferase  (SN2)
  1PWZ/    Haloalcohol dehalogenase  (Cl)
  1RTQ/    Aminopeptidase  (Zn2, 4-step)
  2E7Z/    Acetylene hydratase  (W, 4-step + alt)
  4OTA/    4-Oxalocrotonate tautomerase  (2-step x 2 models)
```

For every enzyme directory we provide:

- `README.md` — one-page description of the system, charges, DFT level of
  theory used in the original paper, and per-step literature barriers.
  This is the upstream Himo-group README, reproduced unchanged from the
  `dataset/Cluster/` working tree.
- `energy_profile.csv` — machine-readable digest of the literature energy
  profile (R -> IM1 -> TS1 -> ... -> P), copied verbatim from
  `dataset/Cluster/<pdb>/energy_profile.csv`. These numbers populate the
  `LITERATURE` dictionary in `scripts/validate_benchmark.py` under the
  left-input baseline convention described in the paper SI-G.
- per-model subdirectories (`model_A/`, `model_B/`, `model_Ba/`, `model/step1/`,
  etc.), each containing:
    - stationary-point XYZ files numbered in reaction order (`1.R.xyz`, optional
      `2.IM*.xyz` / `3.IM*.xyz`, `2.TS.xyz` or `3.TS.xyz`, final `*.P.xyz`).
      Single-step models follow the simpler `1.R.xyz`, `2.TS.xyz`, `3.P.xyz`
      convention; multistep models use an extended numbering with intermediate
      (`IM`) entries — see the actual files in each subdirectory.
    - `config.yaml` — pipeline configuration used for the benchmark run.
    - `run.sh` — driver script invoking one or more of the four pipeline
      conditions (`default`, `nrp`, `tsonly`, `ana`). Different cases ship with
      different conditions active by default; edit `run.sh` to enable/disable
      conditions for your reproduction, then execute `bash run.sh`.

We do **not** redistribute the original PDF reprints of the Himo-group
papers, the Kromann et al. 2016 PeerJ benchmark compilation, or the Himo
2017 JACS cluster-approach review; those remain copyrighted by their
publishers. The citations below give the primary sources.

## How to reproduce one case

```bash
cd examples/benchmark/1AH7/model_A
bash run.sh
```

The run produces a `result/` directory containing the full pipeline
output (MEP search, TS optimization, IRC, frequency analysis, and
optionally a DFT//MLIP single-point stage when `--dft` is active).

Paper timings were measured on a single NVIDIA GH200 120GB GPU (Miyabi
cluster). A consumer NVIDIA RTX 5080 reproduces the bezA application
profile with both transition-state Gibbs energies agreeing to within
~0.1 and ~0.6 kcal/mol respectively, and the intermediate and product
relative Gibbs energies differing by up to 2.2 kcal/mol (paper
Section 4.5 and Fig. 6B caption).

## Qualitative pass/fail judgment

After a case finishes, run `scripts/validate_benchmark.py` (a copy of
`validation2/validate_results.py` with the `VALIDATION_DIR` path adjusted
to this directory) to classify it into one of four display symbols that
match the paper's main-text qualitative table (SI-G, "Pass/Fail Judgment
Criteria"):

- **✓** (clean pass) — all non-metal covalent bond changes in the input
  R→P pair are reproduced in the IRC-relaxed endpoints *and* the
  converged TS has exactly one imaginary mode (`n_imag = 1`);
- **~** (topology pass) — the same bond-change criterion is satisfied,
  but the converged TS carries `n_imag ≥ 2` (reaction coordinate
  identified, TS spectrum not a clean harmonic saddle; typical for the
  Orb backend on this benchmark);
- **×** (fail) — the pipeline completed but the intended reaction is
  not reproduced (bond-change mismatch, no reactive segment detected,
  or shifted endpoints);
- **E** (workflow error) — the pipeline crashed before producing a
  validated TS (e.g. IRC `AssertionError`, SVD non-convergence,
  uncaught exception detected in `all.out`). The workflow-error display
  takes precedence over `×` whenever a crash signature is present,
  even if the validator's status code is otherwise `MULTI_SEGMENT` etc.

Raw validator status codes (`PERFECT`, `GOOD`, `WRONG_REACT`,
`PARTIAL_REACT`, `OPT_BROKEN`, `MULTI_SEGMENT`, `ANOMALOUS`,
`SVD_FAIL`, `ERROR`, `NOT_STARTED`, ...) are kept in the per-run CSV
for diagnostics; the mapping to the four display symbols above is
implemented once in `p2r_main/images/gen_data.py` of the manuscript
source and documented in the paper SI-G.

## Primary citations for the benchmark itself

The 23-step benchmark is curated from two literature sources. Local PDFs
exist under `dataset/Cluster/` in the working tree but are **not**
redistributed under the `pdb2reaction` repository; please cite the
originals.

1. **Kromann, J. C.; Christensen, A. S.; Cui, Q.; Jensen, J. H.**
   Towards a Barrier Height Benchmark Set for Biologically Relevant
   Systems.  *PeerJ* **2016**, *4*, e1994.
   DOI: [10.7717/peerj.1994](https://doi.org/10.7717/peerj.1994).
   (Local working-tree copy: `dataset/Cluster/peerj-1994.pdf`.)
   This is the benchmark compilation that we drew the six enzyme systems
   from; our `LITERATURE` dictionary in `scripts/validate_benchmark.py`
   is consistent with the B3LYP/6-311+G(2d,2p)//B3LYP/6-31G(d,p) barriers
   reported there, under the left-input baseline convention described in
   the paper SI-G.

2. **Himo, F.**  Recent Trends in Quantum Chemical Modeling of Enzymatic
   Reactions.  *J. Am. Chem. Soc.* **2017**, *139*, 6780–6786.
   DOI: [10.1021/jacs.7b02671](https://doi.org/10.1021/jacs.7b02671).
   (Local working-tree copy:
   `dataset/Cluster/the-quantum-chemical-cluster-approach-in-biocatalysis.pdf`.)
   The review article that frames the quantum-chemical cluster approach
   we validate `pdb2reaction` against.

## Per-enzyme primary citations

Each `<pdb>/README.md` lists the primary paper for that system. The
condensed table:

| PDB  | Enzyme                                  | Primary paper |
| :--- | :-------------------------------------- | :------------- |
| 1AH7 | Phospholipase C (Zn3)                   | Liao, Yu, Himo. *J. Phys. Chem. B* **2010**, *114*, 2533–2540. DOI: [10.1021/jp910992f](https://doi.org/10.1021/jp910992f). |
| 1O9S | Histone lysine methyltransferase (SN2)  | Georgieva, Himo. *J. Comput. Chem.* **2010**, *31*, 1707–1714. DOI: [10.1002/jcc.21458](https://doi.org/10.1002/jcc.21458). |
| 1PWZ | Haloalcohol dehalogenase                | Hopmann, Himo. *J. Chem. Theory Comput.* **2008**, *4*, 1129–1137. DOI: [10.1021/ct800135k](https://doi.org/10.1021/ct800135k). |
| 1RTQ | Aminopeptidase (Zn2, 4-step)            | Chen, Xu, Himo. *J. Phys. Chem. B* **2008**, *112*, 2494–2500. DOI: [10.1021/jp710035j](https://doi.org/10.1021/jp710035j). |
| 2E7Z | Acetylene hydratase (W, 4-step)         | Liao, Yu, Himo. *Proc. Natl. Acad. Sci. U.S.A.* **2010**, *107*, 22523–22527. DOI: [10.1073/pnas.1014060108](https://doi.org/10.1073/pnas.1014060108). |
| 4OTA | 4-Oxalocrotonate tautomerase (2-step)   | Sevastik, Himo. *Bioorg. Chem.* **2007**, *35*, 444–457. DOI: [10.1016/j.bioorg.2007.08.003](https://doi.org/10.1016/j.bioorg.2007.08.003). |

The individual `energy_profile.csv` files give the literature energies
for each state on the reported mechanism; each `<pdb>/README.md`
documents charges, models, and the exact computational details used in
the primary paper.

## License / redistribution note

The structured `energy_profile.csv` files and the stationary-point XYZ
inputs derive from geometries deposited by the Himo group in conjunction
with the publications above, and are redistributed here for
reproducibility of the `pdb2reaction` benchmark. The primary-source PDFs
are **not** redistributed; follow the DOI links above for access under
the respective publisher terms. When using this dataset please cite the
per-system primary reference, plus the Kromann 2016 benchmark paper and
the Himo 2017 cluster-approach review.
