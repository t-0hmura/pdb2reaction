## Test Data and Smoke Tests (Developer-Oriented)

The repository includes `tests/smoke/` with sample inputs and a convenience
script for quick smoke checks:

```
tests/smoke/
├── run.sh          # Batch of CLI smoke tests (redirects output to test*.out)
├── *.pdb/*.xyz/*.gjf
└── input.yaml      # Example YAML overrides
```

### Running the smoke tests

From the repository root:

```bash
bash tests/smoke/run.sh
```

Or from inside the `tests/smoke/` directory:

```bash
cd tests/smoke
bash run.sh
```

Notes:
- Outputs (`test1/`, `test2/`, …) and log files (`test1.out`, …) are created in
  `tests/smoke/` (the script changes into its own directory).
- Some tests use GPU UMA models; make sure CUDA, PyTorch, and UMA are available
  and that you have authenticated to Hugging Face (see Getting Started).
- The DFT test uses `--engine cpu` to avoid GPU4PySCF requirements.
