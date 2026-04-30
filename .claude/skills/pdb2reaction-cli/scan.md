# `pdb2reaction scan`

## Purpose

1D bond-length-driven scan with staged harmonic restraints and inter-stage
relaxation. Drive one or more bonds toward target distances, producing a
trajectory you can reuse to seed `path-search`. For most workflows prefer
`pdb2reaction all --scan-lists` (see `all-scan-list.md`); standalone
`scan` is for one-off exploration.

## Synopsis

```bash
pdb2reaction scan -i input.pdb \
    -s '[(idx_a, idx_b, target_A), ...]' \
    [-l 'RES:Q,...'] [-q / -m] \
    [-b uma|orb|mace|aimnet2] [-o ./result_scan/]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Reactant `.pdb` / `.xyz` / `.gjf` |
| `-s, --scan-lists` | str | required | Inline Python literal `'[(a,b,target),...]'`, or YAML/JSON spec path. **Pass multiple stages as space-separated literals after a single `-s`** — repeating `-s` is rejected. |
| `-q` / `-l` / `-m` | — | — | Charge / spin (common conventions) |
| `-b, --backend` | str | `uma` | MLIP backend |
| `--solvent` | str | `none` | xTB-ALPB solvent |
| `-o, --out-dir` | path | `./result_scan/` | Output directory |
| `--ref-pdb` | path | none | Residue context for XYZ/GJF inputs |
| `--config` / `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

The tuple grammar in `-s` accepts atom-index ints (`(1, 5, 1.4)`) or atom
specs (`("CS1 SAM 320", "C7 GPP 321", 1.60)`). Multiple stages chain
sequentially as space-separated literals after a **single** `-s`; each
stage starts from the previous stage's final geometry. Repeating `-s`
is rejected.

## Examples

### Single stage by atom name

```bash
pdb2reaction scan -i 1.R.pdb -l 'SAM:1,GPP:-3' \
    -s '[("CS1 SAM 320","C7 GPP 321",1.60)]' \
    -b uma -o result_scan
```

### Two sequential stages

```bash
pdb2reaction scan -i 1.R.pdb -l 'SAM:1,GPP:-3' \
    -s '[("CS1 SAM 320","C7 GPP 321",1.60)]' \
       '[("H11 GPP 321","OE2 GLU 186",0.90)]' \
    -b uma -o result_scan_staged
```

Each space-separated literal after a single `-s` is one stage; do **not** repeat `-s` (rejected with `repeated flags are not accepted`).

## Output

```
result_scan/
├── result.json                  # only when --out-json is passed
├── preopt/                      # only with --preopt
│   └── result.{xyz,pdb,gjf}     # pre-optimized starting geometry
├── stage_01/                    # per-stage relaxed snapshots
│   ├── result.xyz               # final geometry of stage
│   ├── scan_trj.xyz             # per-stage scan trajectory (always written)
│   └── scan_*.xyz               # intermediate steps (when opt.dump: true)
├── stage_02/                    # …
├── scan_trj.xyz                 # stitched scan trajectory across all stages
└── scan.pdb                     # PDB form (with --convert-files, default on)
```

`result.json` lists per-stage status, target distances, final energies,
and the stitched scan trajectory. Plot with `trj2fig.md`.

## Caveats

- `-s` is Python literal-eval. Quote with single quotes outside,
  double quotes inside. Atom-name strings use `"NAME RESNAME RESID"`
  with single spaces.
- Stage *k+1* starts from stage *k*'s final geometry; a diverged
  stage poisons downstream stages.
- For coupled multi-bond drives in one stage, put multiple tuples in
  one `-s` argument: `'[(a,b,1.6),(c,d,3.0)]'`.

## See also

- `scan2d.md`, `scan3d.md` — higher-dim analogs.
- `all-scan-list.md` — wraps `scan` inside the full pipeline.
- Defaults: `import pdb2reaction.defaults as d; print(d.BIAS_KW, d.BOND_KW, d.OUT_DIR_SCAN)`