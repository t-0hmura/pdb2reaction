# `opt`

## 概要

> **要約:** L-BFGS（`--opt-mode light`、デフォルト）または RFO（`--opt-mode heavy`）を使用して単一構造を局所極小に最適化します。PDB 入力の場合、リンク水素の親原子は自動的に凍結されます。

`pdb2reaction opt` は、UMAがエネルギー、勾配、ヘシアンを提供しながら、pysisyphus LBFGS（"light"）またはRFOptimizer（"heavy"）エンジンで単一構造の構造最適化を実行します。入力構造は `.pdb`、`.xyz`、`.trj`、または `geom_loader` でサポートされる任意の形式が可能です。設定は**組み込みデフォルト → CLI 上書き → `--args-yaml` 上書き**の順序で適用され（YAMLが最も優先）、軽量なデフォルトを維持しながら選択的にオプションを上書きできます。オプティマイザープリセットは現在LBFGSベースの**`light`**モードがデフォルトです。

開始構造がPDBまたはGaussianテンプレートの場合、フォーマット対応変換は最適化された構造を `.pdb`（PDB 入力）および `.gjf`（Gaussianテンプレート）コンパニオンにミラーリングします（`--convert-files {True\|False}` で制御、デフォルトで有効）。
PDB固有の利便性:
- `--freeze-links`（デフォルト `True`）でリンク水素の親原子を検出し、`geom.freeze_atoms` にマージします（0始まり）。
- 出力変換では `final_geometry.pdb`（および `--dump True` の場合は `optimization.pdb`）を入力PDBを参照して書き出します。
XYZ/GJF入力では `--ref-pdb` が参照 PDB トポロジーを提供しXYZ座標を保持するため、フォーマット対応のPDB/GJF出力変換が可能です。

Gaussian `.gjf` テンプレートは電荷/スピンの既定値を与え、変換が有効な場合に最適化構造を `.gjf` として自動出力します。

## 使用法
```bash
pdb2reaction opt -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [-m MULT] \
                 [--opt-mode light|heavy] [--freeze-links {True\|False}] \
                 [--dist-freeze '[(i,j,target_A), ...]'] [--one-based {True\|False}] \
                 [--bias-k K_eV_per_A2] [--dump {True\|False}] [--out-dir DIR] \
                 [--max-cycles N] [--thresh PRESET] [--args-yaml FILE] \
                 [--convert-files {True\|False}] [--ref-pdb FILE]
```

## ワークフロー
- **オプティマイザー**: `--opt-mode light`（デフォルト）→ L-BFGS; `--opt-mode heavy` → 信頼領域制御付きRational Function Optimizer
- **拘束**: `--dist-freeze` はPythonリテラルタプル `(i, j, target_A)` を解釈し、3番目の要素を省略すると開始距離を拘束します。`--bias-k` はグローバル調和強度（eV·Å⁻²）を設定します。インデックスはデフォルトで1始まりですが、`--one-based False` で0始まりに切り替えられます。
- **電荷/スピン解決**: CLI `-q/-m` は `.gjf` テンプレートメタデータを上書きし、それは `calc` デフォルトを上書きします。`-q` が省略され `--ligand-charge` が与えられている場合は酵素–基質複合体として扱い、`extract.py` の電荷サマリーで総電荷を導出します。明示的な `-q` は常に優先され、`.gjf` 以外で `--ligand-charge` が無い場合は中断します。多重度は省略時 `1` がデフォルトです。
- **凍結原子**: CLIのリンク検出はYAMLの `geom.freeze_atoms` とマージされ、UMA 計算機の `calc.freeze_atoms` に反映されます。
- **ダンプ & 変換**: `--dump True` は `opt.dump=True` を反映し `optimization.trj` を出力します。変換が有効な場合、PDB 入力では軌跡が `optimization.pdb` にミラーされます。`opt.dump_restart` を有効にするとリスタートYAMLが出力されます。
- **終了コード**: `0` 成功、`2` ゼロステップ（`min_step_norm` 未満）、`3` 最適化失敗、`130` 割り込み、`1` 予期せぬエラー。

## CLI オプション

> **注記:** 表示されているデフォルト値は、オプション未指定時に使用される内部デフォルトです。

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる入力構造 | 必須 |
| `-q, --charge INT` | 総電荷。`.gjf` テンプレートまたは `--ligand-charge`（PDB 入力または `--ref-pdb` 付きXYZ/GJF）が提供しない限り必須。両方指定時は `-q` が優先 | テンプレート/導出が適用されない限り必須 |
| `--ligand-charge TEXT` | `-q` が省略された場合に使用される総電荷または残基名ごとのマッピング。PDB 入力（または `--ref-pdb` 付きXYZ/GJF）でextract方式の電荷導出を有効化 | _None_ |
| `--workers`, `--workers-per-node` | UMA予測器の並列度（workers > 1 で解析ヘシアン無効; `workers_per_node` は並列予測器へ転送） | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度（2S+1）。`.gjf` テンプレートまたは `1` にフォールバック | テンプレート/`1` |
| `--dist-freeze TEXT` | 調和拘束用の `(i,j,target_A)` タプルを記述するPythonリテラルとして解析される文字列 | _None_ |
| `--one-based {True\|False}` | `--dist-freeze` インデックスを1始まり（デフォルト）または0始まりとして解釈 | `True` |
| `--bias-k FLOAT` | すべての `--dist-freeze` タプルに適用される調和バイアス強度（eV·Å⁻²） | `10.0` |
| `--freeze-links {True\|False}` | リンク水素親凍結をトグル（PDB 入力のみ） | `True` |
| `--max-cycles INT` | 最適化反復のハードリミット | `10000` |
| `--opt-mode TEXT` | オプティマイザー選択: `light`（LBFGS）または `heavy`（RFO） | `light` |
| `--dump {True\|False}` | 軌跡ダンプ（`optimization.trj`）を出力 | `False` |
| `--convert-files {True\|False}` | PDB 入力用のXYZ/TRJ → PDBコンパニオンおよびGaussianテンプレート用のXYZ → GJFコンパニオンを有効/無効化 | `True` |
| `--ref-pdb FILE` | 入力がXYZ/GJFの場合に使用する参照 PDB トポロジー | _None_ |
| `--out-dir TEXT` | すべてのファイルの出力ディレクトリ | `./result_opt/` |
| `--thresh TEXT` | 収束プリセットの上書き（`gau_loose`、`gau`、`gau_tight`、`gau_vtight`、`baker`、`never`） | `gau` |
| `--args-yaml FILE` | YAML 上書きを提供（セクション `geom`、`calc`、`opt`、`lbfgs`、`rfo`） | _None_ |

## 出力
```
out_dir/
├─ final_geometry.xyz          # 常に書き込み
├─ final_geometry.pdb          # 入力がPDBで変換が有効な場合のみ
├─ final_geometry.gjf          # Gaussianテンプレートが検出され変換が有効な場合
├─ optimization.trj            # ダンプが有効な場合のみ
├─ optimization.pdb            # 軌跡のPDB変換（PDB入力、変換有効時）
└─ restart*.yml                # opt.dump_restartが設定されている場合のオプションのリスタート
```
コンソールには解決済みの `geom`/`calc`/`opt`/`lbfgs`/`rfo` ブロックとサイクル進行、総実行時間が出力されます。

(yaml-configuration-args-yaml)=
## YAML 設定（`--args-yaml`）
YAML 値はCLIを上書きし、CLIはデフォルトを上書きします。

### `geom`
- `coord_type`（`"cart"`）: デカルト座標 vs `"dlc"` 非局在化内部座標
- `freeze_atoms`（`[]`）: 0始まりの凍結インデックス; CLIリンク検出と自動マージ

### `calc`
- UMA設定（`model`、`task_name`、デバイス選択、近傍半径、ヘシアン形式など）
- `charge`/`spin` はCLI オプションをミラー（`.gjf` がある場合はテンプレート値が既定）

### `opt`
LBFGSとRFOの両方で使用される共有オプティマイザー制御:
- `thresh` プリセット、`max_cycles`、`print_every`、`min_step_norm`、収束トグル（`rms_force` など）、`converge_to_geom_rms_thresh`、`overachieve_factor`、`check_eigval_structure`、`line_search`。
- ダンプ/管理項目（`dump`、`dump_restart`、`prefix`、`out_dir`）。

### `lbfgs`
L-BFGS固有で `opt` を拡張: `keep_last`、`beta`、`gamma_mult`、`max_step`、`control_step`、`double_damp`、およびオプションの正則化パラメータ `mu_reg`/`max_mu_reg_adaptions`

### `rfo`
RFOptimizerフィールドで `opt` を拡張: 信頼領域サイジング（`trust_radius`、`trust_min`、`trust_max`、`trust_update`）、`max_energy_incr`、ヘシアン管理（`hessian_update`、`hessian_init`、`hessian_recalc`、`hessian_recalc_adapt`、`small_eigval_thresh`）、マイクロイテレーション制御（`alpha0`、`max_micro_cycles`、`rfo_overlaps`）、DIISヘルパー（`gdiis`、`gediis`、閾値、`gdiis_test_direction`）、および `adapt_step_func`


### YAML例
```yaml
geom:
  coord_type: cart           # coordinate type: cartesian vs dlc internals
  freeze_atoms: []           # 0-based frozen atoms merged with CLI/link detection
calc:
  charge: 0                  # total charge (CLI/template override)
  spin: 1                    # spin multiplicity 2S+1
  model: uma-s-1p1           # UMA model tag
  task_name: omol            # UMA task name
  device: auto               # UMA device selection
  max_neigh: null            # maximum neighbors for graph construction
  radius: null               # cutoff radius for neighbor search
  r_edges: false             # store radial edges
  out_hess_torch: true       # request torch-form Hessian
  freeze_atoms: null         # calculator-level frozen atoms
  hessian_calc_mode: FiniteDifference   # Hessian mode selection
  return_partial_hessian: false         # full Hessian (avoids shape mismatches)
opt:
  thresh: gau                # convergence preset (Gaussian/Baker-style)
  max_cycles: 10000          # optimizer cycle cap
  print_every: 100           # logging stride
  min_step_norm: 1.0e-08     # minimum norm for step acceptance
  assert_min_step: true      # stop if steps fall below threshold
  rms_force: null            # explicit RMS force target
  rms_force_only: false      # rely only on RMS force convergence
  max_force_only: false      # rely only on max force convergence
  force_only: false          # skip displacement checks
  converge_to_geom_rms_thresh: 0.05   # geom RMS threshold when converging to ref
  overachieve_factor: 0.0    # factor to tighten thresholds
  check_eigval_structure: false   # validate Hessian eigenstructure
  line_search: true          # enable line search
  dump: false                # dump trajectory/restart data
  dump_restart: false        # dump restart checkpoints
  prefix: ""                 # filename prefix
  out_dir: ./result_opt/     # output directory
lbfgs:
  thresh: gau                # LBFGS convergence preset
  max_cycles: 10000          # iteration limit
  print_every: 100           # logging stride
  min_step_norm: 1.0e-08     # minimum accepted step norm
  assert_min_step: true      # assert when steps stagnate
  rms_force: null            # explicit RMS force target
  rms_force_only: false      # rely only on RMS force convergence
  max_force_only: false      # rely only on max force convergence
  force_only: false          # skip displacement checks
  converge_to_geom_rms_thresh: 0.05   # RMS threshold when targeting geometry
  overachieve_factor: 0.0    # tighten thresholds
  check_eigval_structure: false   # validate Hessian eigenstructure
  line_search: true          # enable line search
  dump: false                # dump trajectory/restart data
  dump_restart: false        # dump restart checkpoints
  prefix: ""                 # filename prefix
  out_dir: ./result_opt/     # output directory
  keep_last: 7               # history size for LBFGS buffers
  beta: 1.0                  # initial damping beta
  gamma_mult: false          # multiplicative gamma update toggle
  max_step: 0.3              # maximum step length
  control_step: true         # control step length adaptively
  double_damp: true          # double damping safeguard
  mu_reg: null               # regularization strength
  max_mu_reg_adaptions: 10   # cap on mu adaptations
rfo:
  thresh: gau                # RFOptimizer convergence preset
  max_cycles: 10000          # iteration cap
  print_every: 100           # logging stride
  min_step_norm: 1.0e-08     # minimum accepted step norm
  assert_min_step: true      # assert when steps stagnate
  rms_force: null            # explicit RMS force target
  rms_force_only: false      # rely only on RMS force convergence
  max_force_only: false      # rely only on max force convergence
  force_only: false          # skip displacement checks
  converge_to_geom_rms_thresh: 0.05   # RMS threshold when targeting geometry
  overachieve_factor: 0.0    # tighten thresholds
  check_eigval_structure: false   # validate Hessian eigenstructure
  line_search: true          # enable line search
  dump: false                # dump trajectory/restart data
  dump_restart: false        # dump restart checkpoints
  prefix: ""                 # filename prefix
  out_dir: ./result_opt/     # output directory
  trust_radius: 0.1          # trust-region radius
  trust_update: true         # enable trust-region updates
  trust_min: 0.0             # minimum trust radius
  trust_max: 0.1             # maximum trust radius
  max_energy_incr: null      # allowed energy increase per step
  hessian_update: bfgs       # Hessian update scheme
  hessian_init: calc         # Hessian initialization source
  hessian_recalc: 200        # rebuild Hessian every N steps
  hessian_recalc_adapt: null # adaptive Hessian rebuild factor
  small_eigval_thresh: 1.0e-08   # eigenvalue threshold for stability
  alpha0: 1.0                # initial micro step
  max_micro_cycles: 50       # micro-iteration limit
  rfo_overlaps: false        # enable RFO overlaps
  gediis: false              # enable GEDIIS
  gdiis: true                # enable GDIIS
  gdiis_thresh: 0.0025       # GDIIS acceptance threshold
  gediis_thresh: 0.01        # GEDIIS acceptance threshold
  gdiis_test_direction: true # test descent direction before DIIS
  adapt_step_func: true      # adaptive step scaling toggle
```

---

## 関連項目

- [tsopt](tsopt.md) — 極小ではなく遷移状態（鞍点）を最適化
- [freq](freq.md) — 最適化が極小に達したことを確認する振動解析
- [extract](extract.md) — 最適化前にポケットPDBを生成
- [all](all.md) — 端点を事前最適化するエンドツーエンドワークフロー
- [YAML リファレンス](yaml-reference.md) — `opt`、`lbfgs`、`rfo` の完全な設定オプション
- [用語集](glossary.md) — L-BFGS、RFOの定義