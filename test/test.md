## Test Data and Smoke Tests

The repository includes a small `test/` folder with sample inputs and a
convenience script for quick smoke checks:

```
pdb2reaction/test/
├── run.sh          # Batch of CLI smoke tests (redirects output to test*.out)
├── *.pdb/*.xyz/*.gjf
└── input.yaml      # Example YAML overrides
```

### Running the smoke tests

From the repository root:

```bash
bash pdb2reaction/test/run.sh
```

Or from inside the `test/` directory:

```bash
cd pdb2reaction/test
bash run.sh
```

Notes:
- Outputs (`test1/`, `test2/`, …) and log files (`test1.out`, …) are created in
  the **current working directory**.
- Some tests use GPU UMA models; make sure CUDA, PyTorch, and UMA are available
  and that you have authenticated to Hugging Face (see Getting Started).
- The DFT test uses `--engine cpu` to avoid GPU4PySCF requirements.