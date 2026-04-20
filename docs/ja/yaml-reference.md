# YAML 設定リファレンス

```{tip}
**`all` 一気通貫コマンドをお探しの場合:** `all` コマンドはこのページに記載のすべてのセクションを読み取り、適切なサブコマンドへ転送します。どのセクションがどのステージで読み取られるかは [`all`](all.md)（特に "サブコマンド → YAML セクション" マッピング）を参照してください。
```

(ja-yaml-configuration-precedence)=
## 設定の優先順位

設定は以下の順序で解決されます（後のものが前のものを上書き）:

```
組み込みデフォルト  <  --config (YAML)  <  CLI フラグ
```

1. **組み込みデフォルト** — `pdb2reaction/defaults.py` に定義されたハードコード値。
2. **`--config`** — デフォルトを上書きする YAML ファイル（例: `--config my_settings.yaml`）。
3. **CLI フラグ** — コマンドラインで明示的に指定したオプション（例: `-q -1`, `--thresh gau_loose`）。*明示的に指定された*値のみが YAML を上書きし、CLI デフォルトのままのオプションは YAML の値を隠しません。

例: YAML で `charge: 0` を設定し、CLI で `-q -1` を渡した場合、電荷は `-1` になります。

この優先順位は `all`, `opt`, `tsopt`, `freq`, `irc`, `scan`, `scan2d`, `scan3d`, `path-opt`, `path-search`, `dft` に共通です。あわせて {ref}`CLI 規約: 設定の優先順位 <ja-configuration-precedence>` を参照してください。

(ja-common-cli-to-yaml-mapping)=
## 主要な CLI→YAML マッピング

| CLI フラグ | YAML キー | セクション |
|----------|----------|---------|
| `-q` / `--charge` | `charge` | `calc` |
| `-m` / `--multiplicity` | `spin` | `calc` |
| `-b` / `--backend` | `backend` | `calc` |
| `--model` | `model` | `calc` |
| `--solvent` | `solvent` | `calc` |
| `--device` | `device` | `calc` |
| `--thresh` | `thresh` | `opt` |
| `--max-cycles` | `max_cycles` | `opt` |
| `--dump` | `dump` | `opt` |
| `--opt-mode` | `opt_mode` | `opt` (tsopt) |
| `--freeze-atoms` | `freeze_atoms` | `geom` |
| `--coord-type` | `coord_type` | `geom` |
| `--temperature`（freq、`all --freq-temperature`） | `temperature` | `thermo` |
| `--pressure`（freq、`all --freq-pressure`） | `pressure_atm` | `thermo` |
| `--engine`（`dft` サブコマンド） / `--dft-engine`（`all` ラッパー） | `engine`（CLI 内部） | `dft` |

```{note}
**名前不一致 — `--pressure` vs `pressure_atm`.** CLI フラグは `--pressure`（単位は暗黙的に atm）、`thermo:` 配下の対応 YAML キーは `pressure_atm`（単位接尾辞付き）です。いずれも atm で扱い、内部で Pa に変換されます。
```

```{note}
**名前不一致 — `--engine` vs `--dft-engine`.** 単体の `dft` サブコマンドでは `--engine`（gpu / cpu）です。`pdb2reaction all` では他の engine 系オプションと衝突を避けるため、同じフラグが `--dft-engine` にリネームされます — {ref}`CLI 規約の --engine vs --dft-engine 節 <ja-engine-vs-dft-engine>` を参照してください。
```

### サブコマンド別の `--thresh` デフォルト

TS 最適化はより厳しい "baker" プリセットを、通常の極小化は "gau" プリセットを使うため、`--thresh` のデフォルトはサブコマンドごとに異なります。

| サブコマンド | デフォルト `--thresh` | 由来となるデフォルト辞書 |
|------------|---------------------|-----------------------|
| `opt` | `gau` | `OPT_BASE_KW`（→ `lbfgs` / `rfo`） |
| `tsopt`（Hessian Dimer） | `baker` | `HESSIAN_DIMER_KW`、内側の `LBFGS_TS_KW` |
| `tsopt`（RS-I-RFO） | `baker` | `RSIRFO_KW` |
| `scan`, `scan2d`, `scan3d` | `gau` | `OPT_BASE_KW` |
| `path-search`（各ステップの opt） | `gau` | `OPT_BASE_KW` |
| `path-opt` / StringOptimizer | `gau_loose` | `STOPT_KW` |
| `all`（pre-opt、post-opt 極小化） | `gau` | `OPT_BASE_KW` |
| `all`（post-opt TS 段階） | `baker` | `HESSIAN_DIMER_KW` / `RSIRFO_KW` |

受け付ける値: `gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`。実行ごとに `--thresh <preset>` または YAML の `opt.thresh` で上書きできます。

```{note}
**`--thresh` を持たないサブコマンド。** `irc`、`freq`、`dft` には `--thresh` が**ありません**:

- `irc` — 収束は `irc.rms_grad_thresh`、`irc.energy_thresh`、`irc.max_cycles` で制御されます（[`irc` セクション](#ja-irc-section) を参照）。IRC は予測子–修正子積分器に従うため、力ベース極小化用のプリセットファミリは適用されません。
- `freq` — 最適化ステップが無いため `--thresh` は存在しません。数値精度は `--hessian-calc-mode` と MLIP 自体の精度で決まります。
- `dft` — SCF 収束は `dft.conv_tol`（デフォルト `1e-9` Hartree）と `dft.max_cycle` で制御されます。`gau`/`baker` プリセットファミリは使用しません。[`dft` セクション](#ja-dft-section) を参照してください。
```

## 概要


| セクション | 説明 | 使用されるコマンド |
|---------|-------------|---------|
| [`geom`](#geom) | ジオメトリと座標設定 | all, opt, scan, scan2d, scan3d, tsopt, freq, irc, path-opt, path-search |
| [`calc`](#calc) | MLIP バックエンドの設定 | all, opt, scan, scan2d, scan3d, tsopt, freq, irc, path-opt, path-search |
| [`opt`](#opt) | 最適化の共通設定 | opt, scan, scan2d, scan3d, tsopt, path-opt, path-search |
| [`lbfgs`](#lbfgs) | L-BFGSの設定 | opt, scan, scan2d, scan3d, path-search |
| [`rfo`](#rfo) | RFOの設定 | opt, scan, scan2d, scan3d, path-search |
| [`gs`](#gs) | GSM（Growing String Method）設定 | path-opt, path-search |
| [`dmf`](#dmf) | DMF（Direct Max Flux）設定 | path-opt, path-search |
| [`stopt`](#stopt) | StringOptimizer 設定 | path-opt, path-search |
| [`irc`](#ja-irc-section) | IRC積分設定 | irc |
| [`freq`](#ja-freq-section) | 振動解析設定 | freq |
| [`thermo`](#thermo) | 熱化学設定 | freq |
| [`dft`](#ja-dft-section) | DFT計算設定 | dft |
| [`bias`](#bias) | 調和バイアス設定 | scan, scan2d, scan3d |
| [`bond`](#bond) | 結合変化検出設定 | scan, path-search |
| [`search`](#search) | 再帰的経路探索設定 | path-search |
| [`hessian_dimer`](#hessian_dimer) | Hessian Guided Dimer TS 最適化 | tsopt |
| [`rsirfo`](#rsirfo) | RS-I-RFO TS 最適化 | tsopt |

---

## 共通セクション

### `geom`

ジオメトリ読み込みと座標系の設定。

```yaml
geom:
 coord_type: cart # Coordinate type: "cart" (Cartesian) or "dlc" (delocalized internals)
 freeze_atoms: [] # 1-based indices of atoms to freeze; merged with CLI --freeze-links detection
```

**注記:**
- `freeze_atoms` は PDB 入力時の `--freeze-links` 検出原子とマージされます。
- 凍結原子は力がゼロ化され、ヘシアンの該当行/列もゼロ化されます。
- `irc` では `geom.coord_type` が YAML/CLI マージ後に `cart` へ強制されます。

---

### `calc`

MLIP バックエンドの設定。複数バックエンド（UMA, ORB, MACE, AIMNet2）と xTB 溶媒補正に対応。

```yaml
calc:
 backend: uma           # MLIP backend: "uma", "orb", "mace", or "aimnet2"
 charge: 0 # Total system charge (overridden by CLI -q)
 spin: 1 # Spin multiplicity 2S+1 (overridden by CLI -m)
 model: uma-s-1p1 # uma-s-1p1 | uma-m-1p1
 task_name: omol # Task tag recorded in UMA batches
 device: auto # Device: "cuda", "cpu", or "auto"
 max_neigh: null # Maximum neighbors for graph construction
 radius: null # Cutoff radius for neighbor search
 r_edges: false # Store radial edges
 workers: 1 # MLIP inference workers (workers>1 disables analytical Hessians; UMA backend only)
 workers_per_node: 1 # Workers per node for parallel predictor
 out_hess_torch: true # Return Hessian as torch.Tensor
 hessian_double: true # Assemble/return Hessian in float64
 # freeze_atoms: null # geom.freeze_atoms から継承されるため直接指定しない
 hessian_calc_mode: FiniteDifference # Hessian mode: "Analytical" or "FiniteDifference"
 return_partial_hessian: true  # active-DOF ブロックのヘシアンを返す
 print_timing: true # ヘシアン計算のタイミング内訳を表示
 print_vram: true # ヘシアン計算中の CUDA VRAM 使用量を表示 (UMA バックエンドのみ)
 # Solvent correction (xTB)
 solvent: none           # Implicit solvent name (e.g. "water", "methanol") or "none" to disable
 solvent_model: alpb     # xTB solvent model: "alpb" or "cpcmx"
 xtb_cmd: xtb            # Path to xTB executable
 xtb_acc: 0.2            # xTB accuracy parameter
```

**注記:**
- `backend` で MLIP エンジンを選択。UMA（デフォルト）は解析ヘシアンとマルチワーカー推論に対応。他のバックエンドは有限差分ヘシアンを使用。
- `workers` / `workers_per_node` は UMA バックエンドでのみ有効。
- `solvent` で xTB ベースの暗黙溶媒補正を有効化（デルタ補正方式）。`xtb` のインストールが必要。
- VRAM が十分な場合は `hessian_calc_mode: Analytical` を使用してください。
- `workers > 1` の場合、解析ヘシアンは無効化されます — `hessian_calc_mode: Analytical` を明示指定していても、`workers > 1` では **警告なしに有限差分へダウングレード**されます。デフォルトがそもそも `FiniteDifference` のため、通常問題になるのは `Analytical` を明示的に選択したときだけです。詳細は {ref}`MLIP Calculator のヘシアン評価モード <ja-hessian-evaluation>` を参照してください。
- 電荷/スピンは `.gjf` テンプレートがあればそれを継承します。
- `freq` はデフォルトで `calc.return_partial_hessian = true`（PHVA）を設定します（YAML で上書き可能）。
- IRC は `geom.coord_type = cart` と `calc.return_partial_hessian = true` を常に強制します（YAMLより優先、partial Hessian で active-DOF 処理）。
- `irc` では `calc.return_partial_hessian` が YAML/CLI マージ後に `true` へ強制されます。

---

### `opt`

L-BFGS/RFO で共通の単一構造最適化設定。

```yaml
opt:
 thresh: gau # Convergence preset: gau_loose, gau, gau_tight, gau_vtight, baker, never
 max_cycles: 10000 # Maximum optimizer iterations
 print_every: 100 # Logging stride
 min_step_norm: 1.0e-08 # Minimum step norm for acceptance
 assert_min_step: true # Stop if steps fall below threshold
 rms_force: null # Explicit RMS force target
 rms_force_only: false # Rely only on RMS force convergence
 max_force_only: false # Rely only on max force convergence
 force_only: false # Skip displacement checks
 converge_to_geom_rms_thresh: 0.05 # RMS threshold when converging to reference geometry
 overachieve_factor: 0.0 # Factor to tighten thresholds
 check_eigval_structure: false # Validate Hessian eigenstructure
 line_search: true # Enable line search
 energy_plateau: true # エネルギーがフラット化した場合にフォールバック収束を宣言（下記注記を参照）
 energy_plateau_thresh: 1.0e-04 # au (~0.06 kcal/mol); プラトー判定のレンジ閾値
 energy_plateau_window: 50 # プラトー判定に用いる直近ステップ数
 dump: false # Dump trajectory/restart data
 dump_restart: false # Dump restart checkpoints
 prefix: "" # Filename prefix
 out_dir: ./result_opt/ # Output directory
```

**エネルギープラトーによるフォールバック収束 (v0.3.5 新機能):**
`energy_plateau: true` の場合、直近 `energy_plateau_window` ステップのエネルギーレンジ
（max − min）が `energy_plateau_thresh`（デフォルト `1×10⁻⁴ au ≈ 0.06 kcal/mol`、50 ステップ）
を下回ると、オプティマイザーは収束したと判定します。MLIP の力のノイズフロア
（典型的には ~4×10⁻⁴ au）が力ベースの収束閾値（例: `baker` max_force = 3×10⁻⁴ au）を
上回る場合でも、エネルギーが明らかにフラット化していれば無駄なサイクルを消費せずに
停止できます。
chain-of-states（COS）オプティマイザー（`stopt`、`gs`、DMF など）はイメージごとの
エネルギー配列を保持するため、このフォールバックは**スキップ**されます。

**収束プリセット:**

| Preset | Max Force | RMS Force | Max Step | RMS Step |
|--------|-----------|-----------|----------|----------|
| `gau_loose` | 2.5e-3 | 1.7e-3 | 1.0e-2 | 6.7e-3 |
| `gau` | 4.5e-4 | 3.0e-4 | 1.8e-3 | 1.2e-3 |
| `gau_tight` | 1.5e-5 | 1.0e-5 | 6.0e-5 | 4.0e-5 |
| `gau_vtight` | 2.0e-6 | 1.0e-6 | 6.0e-6 | 4.0e-6 |
| `baker` | 3.0e-4 | 2.0e-4 | 3.0e-4 | 2.0e-4 |

---

### `lbfgs`

L-BFGSの設定（`opt` を拡張）。

```yaml
lbfgs:
 # Inherits all opt settings, plus:
 keep_last: 7 # History size for L-BFGS buffers
 beta: 1.0 # Initial damping beta
 gamma_mult: false # Multiplicative gamma update toggle
 max_step: 0.3 # Maximum step length
 control_step: true # Control step length adaptively
 double_damp: true # Double damping safeguard
 mu_reg: null # Regularization strength
 max_mu_reg_adaptions: 10 # Cap on mu adaptations
```

---

### `rfo`

RFO（Rational Function Optimizer）の設定（`opt` を拡張）。

```yaml
rfo:
 # Inherits all opt settings, plus:
 trust_radius: 0.10 # Trust-region radius
 trust_update: true # Enable trust-region updates
 trust_min: 0.0001 # Minimum trust radius
 trust_max: 0.10 # Maximum trust radius (bohr); v0.3.5 で 0.20 から 0.10 に変更
 max_energy_incr: null # Allowed energy increase per step
 hessian_update: bfgs # Hessian update scheme: bfgs, bofill, etc.
 hessian_init: calc # Hessian initialization: calc, unit, etc.
 hessian_recalc: 500 # Rebuild Hessian every N steps
 hessian_recalc_adapt: null # Adaptive Hessian rebuild factor
 small_eigval_thresh: 1.0e-08 # Eigenvalue threshold for stability
 alpha0: 1.0 # Initial micro step
 max_micro_cycles: 50 # Micro-iteration limit
 rfo_overlaps: false # Enable RFO overlaps
 gediis: false # Enable GEDIIS
 gdiis: true # Enable GDIIS
 gdiis_thresh: 0.0025 # GDIIS acceptance threshold
 gediis_thresh: 0.01 # GEDIIS acceptance threshold
 gdiis_test_direction: true # Test descent direction before DIIS
 adapt_step_func: true # Adaptive step scaling
```

---

## 経路最適化セクション

### `gs`

Growing String Method（GSM）の設定。

```yaml
gs:
 fix_first: true # Keep first endpoint fixed
 fix_last: true # Keep last endpoint fixed
 max_nodes: 20 # Maximum string nodes (internal images); GSM ではエンドポイント2点を加えた合計が総画像数
 perp_thresh: 0.005 # Perpendicular displacement threshold
 reparam_check: rms # Reparameterization check metric
 reparam_every: 1 # Reparameterization stride
 reparam_every_full: 1 # Full reparameterization stride
 param: equi # Parameterization scheme
 max_micro_cycles: 10 # Micro-iteration limit
 reset_dlc: true # Rebuild delocalized coordinates each step
 climb: true # Enable climbing image
 climb_rms: 0.0005 # Climbing RMS threshold
 climb_lanczos: true # Lanczos refinement for climbing
 climb_lanczos_rms: 0.0005 # Lanczos RMS threshold
 climb_fixed: false # Keep climbing image fixed
 scheduler: null # Optional scheduler backend
```

```{note}
**GSM と DMF で `max_nodes` の意味が異なります。** **GSM**（`mep-mode gs`）では `gs.max_nodes` は **内部画像数** を指し、エンドポイント2点が固定されるため総パス長は `max_nodes + 2` となります。**DMF**（`mep-mode dmf`）では CLI の `--max-nodes` は **可動画像数** を指し、エンドポイントの暗黙的な追加はありません。詳細は [`path-opt`](path-opt.md) を参照してください。
```

---

### `dmf`

Direct Max Flux（DMF）によるMEP最適化。

```{note}
**DMF の `--max-nodes` は「可動画像数」を意味します** — GSM と異なり、DMF は固定エンドポイント2点を加算しません。上記の `gs.max_nodes` の注記と比較してください。
```

```yaml
dmf:
 max_cycles: 300 # DMF/IPOPT の最大反復数（--max-cycles で上書き）
 correlated: true # Correlated DMF propagation
 sequential: true # Sequential DMF execution
 fbenm_only_endpoints: false # Run FB-ENM beyond endpoints
 fbenm_options:
 delta_scale: 0.2 # FB-ENM displacement scaling
 bond_scale: 1.25 # Bond cutoff scaling
 fix_planes: true # Enforce planar constraints
 cfbenm_options:
 bond_scale: 1.25 # CFB-ENM bond cutoff scaling
 corr0_scale: 1.1 # Correlation scale for corr0
 corr1_scale: 1.5 # Correlation scale for corr1
 corr2_scale: 1.6 # Correlation scale for corr2
 eps: 0.05 # Correlation epsilon
 pivotal: true # Pivotal residue handling
 single: true # Single-atom pivots
 remove_fourmembered: true # Prune four-membered rings
 dmf_options:
 remove_rotation_and_translation: false # Keep rigid-body motions
 mass_weighted: false # Toggle mass weighting
 parallel: false # Enable parallel DMF
 eps_vel: 0.01 # Velocity tolerance
 eps_rot: 0.01 # Rotational tolerance
 beta: 10.0 # Beta parameter for DMF
 update_teval: false # Update transition evaluation
 k_fix: 300.0 # Harmonic constant for restraints
```

---

### `search`

再帰的経路探索（path-searchのみ）。

```yaml
search:
 max_depth: 10 # Recursion depth limit
 stitch_rmsd_thresh: 0.0001 # RMSD threshold for stitching segments
 bridge_rmsd_thresh: 0.0001 # RMSD threshold for bridging nodes
 max_nodes_segment: 10 # Max nodes per segment
 max_nodes_bridge: 5 # Max nodes per bridge
 kink_max_nodes: 3 # Max nodes for kink optimizations
 max_seq_kink: 2 # Max sequential kinks
 refine_mode: null # Refinement strategy: peak, minima, or null (auto)
```

---

### `stopt`

chain-of-states 経路最適化（`path-opt`, `path-search`）向けの StringOptimizer 設定。

```yaml
stopt:
 type: string # Optimizer type label
 thresh: gau_loose # StringOptimizer convergence preset
 stop_in_when_full: 300 # Early stop threshold when the string is full
 align: false # Alignment toggle (forced to False in path-opt/path-search; external Kabsch alignment is used instead)
 scale_step: global # Step scaling mode
 max_cycles: 300 # Maximum StringOptimizer iterations
 dump: false # Dump trajectory/restart data
 dump_restart: false # Dump restart checkpoints
 reparam_thresh: 0.0 # Reparameterization threshold
 coord_diff_thresh: 0.0 # Coordinate-difference threshold
 out_dir: ./result_path_opt/ # Output directory
 print_every: 10 # Logging stride
```

---

## TS 最適化セクション

TS 最適化は `--opt-mode` で**2 つのアルゴリズム**を切り替えます:
- `--opt-mode dimer`（または `grad`）→ `hessian_dimer` セクション
- `--opt-mode rsirfo`（または `hess`、デフォルト）→ `rsirfo` セクション

共通オプティマイザ設定（`thresh`, `max_cycles`, `dump`）は上記 `opt` セクションから読み込まれます。

### `hessian_dimer`

Hessian Guided Dimer TS 最適化（tsopt --opt-mode grad）。

```yaml
hessian_dimer:
 thresh_loose: gau_loose # Loose convergence preset
 thresh: baker # Main convergence preset
 update_interval_hessian: 500 # Hessian rebuild cadence
 neg_freq_thresh_cm: 5.0 # Imaginary-frequency detection threshold (cm⁻¹)
 flatten_amp_ang: 0.1 # Flattening amplitude (Å)
 flatten_max_iter: 50 # Flattening iteration cap (下記注記を参照)
 flatten_sep_cutoff: 0.0 # Minimum distance between representative atoms
 flatten_k: 10 # Representative atoms sampled per mode
 flatten_loop_bofill: false # Bofill update for flatten displacements
 mem: 100000 # Memory limit for solver
 device: auto # Device selection for eigensolver
 root: 0 # Targeted TS root index
 dimer:
 length: 0.0189 # Dimer separation (Bohr)
 rotation_max_cycles: 15 # Max rotation iterations
 rotation_method: fourier # Rotation optimizer method
 rotation_thresh: 0.0001 # Rotation convergence threshold
 rotation_tol: 1 # Rotation tolerance factor
 rotation_max_element: 0.001 # Max rotation matrix element
 rotation_interpolate: true # Interpolate rotation steps
 rotation_disable: false # Disable rotations entirely
 rotation_disable_pos_curv: true # Disable when positive curvature detected
 rotation_remove_trans: true # Remove translational components
 trans_force_f_perp: true # Project forces perpendicular to translation
 bonds: null # Bond list for constraints
 N_hessian: null # Hessian size override
 bias_rotation: false # Bias rotational search
 bias_translation: false # Bias translational search
 bias_gaussian_dot: 0.1 # Gaussian bias dot product
 seed: null # RNG seed for rotations
 write_orientations: true # Write rotation orientations
 forward_hessian: true # Propagate Hessian forward
 lbfgs:
 # Same keys as lbfgs section
 thresh: baker
 max_cycles: 10000
```

```{note}
**`flatten_max_iter` の CLI 優先順位の例外。** YAML の `hessian_dimer.flatten_max_iter: 50`（および対応する `rsirfo.flatten_max_iter`）は、通常の `defaults < YAML < CLI` の順序とは異なり、**コマンドラインで `--flatten` が明示的に渡されない限り CLI によって `0` に上書きされます**。完全な動作表は {ref}`ja-flatten-precedence-caveat` を参照してください。
```

---

### `rsirfo`

RS-I-RFO TS 最適化（tsopt --opt-mode hess）。

```yaml
rsirfo:
 thresh: baker # RS-IRFO convergence preset
 max_cycles: 10000 # Iteration cap
 print_every: 100 # Logging stride
 min_step_norm: 1.0e-08 # Minimum accepted step norm
 assert_min_step: true # Assert when steps stagnate
 roots: [0] # Target root indices
 hessian_ref: null # Reference Hessian
 rx_modes: null # Reaction-mode definitions
 prim_coord: null # Primary coordinates to monitor
 rx_coords: null # Reaction coordinates to monitor
 hessian_update: bofill # Hessian update scheme
 hessian_recalc: 500 # Rebuild exact Hessian every N macro steps (rfo から継承)
 hessian_recalc_reset: true # Reset recalc counter after exact Hessian
 max_micro_cycles: 50 # Micro-iterations per macro cycle
 augment_bonds: false # Augment reaction path based on bond analysis
 min_line_search: true # Enforce minimum line-search step
 max_line_search: true # Enforce maximum line-search step
 assert_neg_eigval: false # Require negative eigenvalue at convergence
 track_mode_by_overlap: false # 前回の Hessian との重なりで追跡対象 TS モードを選ぶ
 # Also inherits rfo-like settings: trust_radius, trust_update, etc.
```

```{note}
**`rsirfo.flatten_max_iter` の CLI 優先順位の例外。** YAML で `rsirfo.flatten_max_iter` を設定しても、コマンドラインで `--flatten` が明示的に渡されない限り、CLI によって `0` に上書きされます。{ref}`ja-flatten-precedence-caveat` を参照してください。
```

---

## IRCセクション

(ja-irc-section)=
### `irc` (section)

IRC積分設定。

```yaml
irc:
 step_length: 0.1 # Integration step length
 max_cycles: 125 # Maximum steps along IRC
 downhill: false # Follow downhill direction only
 forward: true # Propagate in forward direction
 backward: true # Propagate in backward direction
 root: 0 # Normal-mode root index
 hessian_init: calc # Hessian initialization source
 hessian_update: bofill # Hessian update scheme
 hessian_recalc: null # Hessian rebuild cadence
 displ: energy # Displacement construction method
 displ_energy: 0.001 # Energy-based displacement scaling
 displ_length: 0.1 # Length-based displacement fallback
 rms_grad_thresh: 0.001 # RMS gradient convergence threshold
 hard_rms_grad_thresh: null # Hard RMS gradient stop
 energy_thresh: 0.000001 # Energy change threshold
 imag_below: 0.0 # Imaginary frequency cutoff
 force_inflection: true # Enforce inflection detection
 check_bonds: false # Check bonds during propagation
 out_dir: ./result_irc/ # Output directory
 prefix: "" # Filename prefix
 max_pred_steps: 500 # Predictor-corrector max steps
 loose_cycles: 3 # Loose cycles before tightening
 corr_func: mbs # EulerPC コレクタ関数: "mbs"（Modified Bulirsch–Stoer、デフォルト）または "rk4"
```

`corr_func` は、予測子–修正子法ベースの IRC 積分器（EulerPC）が使う修正子ステップを選択します。`"mbs"` は pysisyphus 組み込みの Modified Bulirsch–Stoer 実装（デフォルト）、`"rk4"` は古典的な 4 次 Runge–Kutta 修正子を要求します。既定の積分器がシステム上で数値的に不安定な場合にのみ変更してください。通常は `mbs` のままで問題ありません。

---

## 振動解析セクション

(ja-freq-section)=
### `freq` (section)

振動解析設定。

```yaml
freq:
 amplitude_ang: 0.8 # Displacement amplitude for modes (Å)
 n_frames: 20 # Number of frames per mode animation
 max_write: 10 # Maximum number of modes to write
 sort: value # Sort order: "value" or "abs"
 out_dir: ./result_freq/ # Output directory
```

---

### `thermo`

熱化学設定。

```yaml
thermo:
 temperature: 298.15 # Thermochemistry temperature (K)
 pressure_atm: 1.0 # Thermochemistry pressure (atm)
 dump: false # Write thermoanalysis.yaml
```

---

## DFTセクション

(ja-dft-section)=
### `dft` (section)

DFT計算設定。

```yaml
dft:
 func: wb97m-v # Exchange-correlation functional
 basis: def2-tzvpd # Basis set name
 func_basis: null # Combined "FUNC/BASIS" string (overrides func/basis)
 conv_tol: 1.0e-09 # SCF convergence tolerance (hartree)
 max_cycle: 100 # Maximum SCF iterations
 grid_level: 3 # PySCF grid level
 verbose: 0 # PySCF verbosity (0-9)
 out_dir: ./result_dft/ # Output directory root
```

---

## スキャン関連セクション

スキャン座標は `-s/--scan-lists`（インラインまたは YAML ファイル）で指定します（メイン YAML 設定ではありません）。
構文の詳細は[クイックスタート: スキャン](quickstart-scan.md)を参照してください。

(ja-bias-section)=
### `bias`

スキャン・拘束付き最適化で使用する調和バイアス設定。

```yaml
bias:
 k: 300.0 # Harmonic bias strength (eV·Å⁻²)
```

**サブコマンド間で共有されるばね定数。** 同じ物理的な調和ペナルティ（`k`、単位 eV·Å⁻²）がデフォルト値 `300.0` で 3 箇所に現れます:

| YAML キー | 使用元 | CLI フラグ |
|----------|-------|-----------|
| `bias.k` | `opt`（`--dist-freeze` 原子ペアに対する `--bias-k`）、`scan`, `scan2d`, `scan3d` | `--bias-k` |
| `dmf.dmf_options.k_fix` | `path-opt` / `path-search` で `mep_mode: dmf` を使用する場合 | —（YAML 専用） |

調和拘束の強さを調整したい場合はこれらのいずれかを上書きしてください。値を小さく（例: `20.0`）すると、柔らかい誘導項としてジオメトリが緩和しやすくなります。デフォルトの `300.0` はほぼ剛体的に固定する値です。

---

### `bond`

MLIPベースの結合変化検出。

```yaml
bond:
 device: auto # MLIP device for bond analysis
 bond_factor: 1.2 # Covalent-radius scaling for cutoff
 margin_fraction: 0.05 # Fractional tolerance for comparisons
 delta_fraction: 0.05 # Minimum relative change to flag bond formation/breaking
```

---

## 例: 設定ファイルの全体例

```yaml
# pdb2reaction configuration example

geom:
 coord_type: cart
 freeze_atoms: []

calc:
 backend: uma
 charge: 0
 spin: 1
 model: uma-s-1p1 # uma-s-1p1 | uma-m-1p1
 device: auto
 hessian_calc_mode: Analytical # Recommended when VRAM permits
 solvent: none                 # Set to e.g. "water" for implicit solvent

gs:
 max_nodes: 12
 climb: true
 climb_lanczos: true

stopt:
 thresh: gau_loose
 max_cycles: 300
 dump: false
 out_dir: ./result_all/

opt:
 thresh: gau
 lbfgs:                     # トップレベルの `lbfgs:` としても記述可
 max_cycles: 10000
 rfo:                       # トップレベルの `rfo:` としても記述可
 max_cycles: 10000

bond:
 bond_factor: 1.2
 delta_fraction: 0.05

search:
 max_depth: 10
 max_nodes_segment: 10

freq:
 max_write: 10
 amplitude_ang: 0.8

thermo:
 temperature: 298.15
 pressure_atm: 1.0

dft:
 func: wb97m-v
 basis: def2-tzvpd
 grid_level: 3
```

---

## 参照

- [all](all.md) - 一気通貫ワークフロー
- [opt](opt.md) - 単一構造最適化
- [tsopt](tsopt.md) - 遷移状態最適化
- [path-search](path-search.md) - 再帰的 MEP 探索
- [freq](freq.md) - 振動解析
- [dft](dft.md) - DFT計算
- [uma-pysis](uma-pysis.md) - MLIP バックエンドの詳細
