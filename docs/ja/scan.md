# `scan`

## 概要

> **要約:** 調和拘束を使用して結合距離をスキャンし、反応座標を駆動します。`--scan-lists` でターゲット距離を指定します。複数ステージは順次実行され、各ステージは前の結果から開始します。

`scan` は、UMA 計算機と調和拘束を使った**段階的な結合長スキャン**を実行します。各タプル `(i, j, targetÅ)` が距離ターゲットを定義します。各積分ステップで一時ターゲットを更新し、拘束井戸を適用したうえで、構造全体を LBFGS（`--opt-mode` light、デフォルト）または RFOptimizer（`--opt-mode` heavy）で緩和します。バイアス付きの走査後、書き出す構造を整えるために無バイアスの前処理/後処理最適化を任意で追加できます。
`--scan-lists` が一度だけ指定された場合は単一ステージ、複数リテラルを与えると**順次ステージ**として実行され、各ステージは直前の緩和構造から開始します。
XYZ/GJF入力では、`--ref-pdb` が参照 PDB トポロジーを提供し、XYZ座標は保持されるため、PDB/GJF出力変換が可能になります。

## 使用法
```bash
pdb2reaction scan -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [-m MULT] \
                  --scan-lists '[(i,j,targetÅ), ...]' [options]
                  [--convert-files {True\|False}] [--ref-pdb FILE]
```

### 例
```bash
# 単一ステージの最小例
pdb2reaction scan -i input.pdb -q 0 --scan-lists '[("TYR,285,CA","MMT,309,C10",1.35)]'

# 2ステージ、LBFGS緩和、軌跡ダンプ
pdb2reaction scan -i input.pdb -q 0 --scan-lists \
    '[("TYR,285,CA","MMT,309,C10",1.35)]' \
    '[("TYR,285,CA","MMT,309,C10",2.20),("TYR,285,CB","MMT,309,C11",1.80)]' \
    --max-step-size 0.20 --dump True --out-dir ./result_scan/ --opt-mode light \
    --preopt True --endopt True

# 単一の --scan-lists の後に複数リテラルを渡す
pdb2reaction scan -i input.pdb -q 0 --scan-lists \
    '[("TYR,285,CA","MMT,309,C10",1.35)]' \
    '[("TYR,285,CA","MMT,309,C10",2.20),("TYR,285,CB","MMT,309,C11",1.80)]'
```

## ワークフロー
1. `geom_loader` で構造を読み込み、CLIの上書き・埋め込みGaussianテンプレート（存在する場合）・デフォルトから電荷/スピンを解決します。`-q` が省略され `--ligand-charge` が与えられている場合、入力を酵素–基質複合体として扱い、PDB 入力（または `--ref-pdb` 付きXYZ/GJF）では `extract.py` の電荷サマリーから総電荷を導出します。
2. `--preopt True` が指定された場合、バイアスをかける前に無バイアスの前処理最適化を実行します。
3. `--scan-lists` で与えられた各ステージリテラルについて `(i, j)` を解析・正規化（デフォルトは1始まり）。PDB 入力では、各エントリは整数インデックスまたは `'TYR,285,CA'` のような原子セレクタ文字列を指定できます。セレクタは空白/カンマ/スラッシュ/バッククォート/バックスラッシュで区切れ、トークン順序は任意（フォールバックは resname, resseq, atom を想定）。
   各結合について `Δ = target − current` を計算し、`h = --max-step-size` として `N = ceil(max(|Δ|) / h)` に分割します。各結合は `δ = Δ / N` ずつ更新されます。
4. すべてのステップを進めながら一時ターゲットを更新し、調和井戸 `E = Σ ½ k (|ri − rj| − target)²` を適用してUMAで最小化します。最適化サイクルは `--relax-max-cycles` で上限が設定され、YAML 値を上書きします。
5. 各ステージの最終ステップ後に、必要に応じて無バイアス緩和（`--endopt True`）を実行し、共有結合の変化を報告して `result.*` を出力します。
6. すべてのステージで繰り返し、`--dump True` の場合のみ軌跡を保存します。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル | 必須 |
| `-q, --charge INT` | 総電荷（CLI > テンプレート > 0）。`-q` が省略され `--ligand-charge` がある場合は電荷導出が行われ、明示的な `-q` が最優先 | `.gjf` テンプレートまたは `--ligand-charge` がない場合は必須 |
| `--ligand-charge TEXT` | `-q` が省略された場合に使用される総電荷または残基名ごとのマッピング。PDB 入力（または `--ref-pdb` 付きXYZ/GJF）で extract 方式の電荷導出を有効化 | _None_ |
| `--workers`, `--workers-per-node` | UMA予測器の並列度（workers > 1 で解析ヘシアンは無効化; `workers_per_node` は並列予測器へ転送） | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度 2S+1（CLI > テンプレート > 1） | `.gjf` テンプレート値または `1` |
| `--scan-lists TEXT` | `(i,j,targetÅ)` タプルを含むPythonリテラル。各リテラルが1ステージ; 1つのフラグの後に複数リテラルを渡す。`i`/`j` は整数インデックスまたは PDB 原子セレクタ（`'TYR,285,CA'`） | 必須 |
| `--one-based {True\|False}` | 原子インデックスを1始まり/0始まりとして解釈 | `True` |
| `--max-step-size FLOAT` | 1ステップあたりの最大距離変化（Å）。ステップ数を決定 | `0.20` |
| `--bias-k FLOAT` | 調和バイアス強度 `k`（eV·Å⁻²）。`bias.k` を上書き | `100` |
| `--relax-max-cycles INT` | 前処理/各バイアスステップ/後処理の最適化サイクル上限。`opt.max_cycles` を上書き | `10000` |
| `--opt-mode TEXT` | `light` → LBFGS、`heavy` → RFOptimizer | `light` |
| `--freeze-links {True\|False}` | PDB 入力時にリンク水素の親を凍結 | `True` |
| `--dump {True\|False}` | バイアス付き軌跡（`scan.trj`/`scan.pdb`）を出力 | `False` |
| `--convert-files {True\|False}` | PDB/Gaussian入力で XYZ/TRJ → PDB/GJF コンパニオン変換を切り替え（軌跡変換はPDBのみ） | `True` |
| `--ref-pdb FILE` | XYZ/GJF入力時の参照 PDB トポロジー（XYZ座標は保持） | _None_ |
| `--out-dir TEXT` | 出力ディレクトリ | `./result_scan/` |
| `--thresh TEXT` | 収束プリセット上書き（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） | `gau` |
| `--args-yaml FILE` | YAML 上書き（`geom`, `calc`, `opt`, `lbfgs`, `rfo`, `bias`, `bond`） | _None_ |
| `--preopt {True\|False}` | スキャン前に無バイアス最適化を実行 | `True` |
| `--endopt {True\|False}` | 各ステージ後に無バイアス最適化を実行 | `True` |

### 共有YAMLセクション
- `geom`, `calc`, `opt`, `lbfgs`, `rfo`: [YAML リファレンス](yaml-reference.md) と同じキー。`opt.dump` は内部的に `False` に固定されるため、ステージ軌跡は `--dump` で制御します。
- `--relax-max-cycles` は明示的に指定された場合のみ `opt.max_cycles` を上書き。未指定の場合はYAMLの `opt.max_cycles` が使われます（デフォルト `10000`）。

### セクション `bias`
- `k`（`100`）: 調和バイアス強度（eV·Å⁻²）。

(section-bond)=
### セクション `bond`
`path_search` と同一のUMAベース結合変化検出:
- `device`（`"cuda"`）: 結合解析用UMAデバイス。
- `bond_factor`（`1.20`）: 共有結合半径のスケーリング。
- `margin_fraction`（`0.05`）: 比較用の相対許容。
- `delta_fraction`（`0.05`）: 形成/切断を判定する最小相対変化。

## 出力
```
out_dir/ (デフォルト: ./result_scan/)
├─ preopt/                   # --preopt が True の場合
│  ├─ result.xyz
│  ├─ result.pdb             # PDB入力かつ変換有効時
│  └─ result.gjf             # Gaussianテンプレートがあり変換有効時
└─ stage_XX/                 # ステージごとのフォルダ
    ├─ result.xyz
    ├─ result.pdb             # 最終構造のPDBミラー（変換有効時）
    ├─ result.gjf             # テンプレートがある場合のGaussianミラー（変換有効時）
    ├─ scan.trj               # --dump True の場合
    └─ scan.pdb               # PDB入力で変換有効時の軌跡コンパニオン（scan.gjf は生成しない）
```
- `geom`/`calc`/`opt`/`bias`/`bond` と最適化ブロックの解決結果、および各ステージの結合変化レポートがコンソールに出力されます。

## 注意事項
- `--scan-lists` は単一フラグの後に複数リテラルを与えてください。フラグの繰り返しは非対応。ターゲット距離は正の値が必要です。原子インデックスは内部で0始まりに正規化されます。PDB 入力ではセレクタ文字列で空白/カンマ/スラッシュ/バッククォート/バックスラッシュ区切りを許容し、トークン順序は任意です。
- `--freeze-links` はユーザー指定の `freeze_atoms` に、PDBのリンクH親原子を追加してポケットを固定します。
- 電荷/スピンはテンプレートに継承されます。`.gjf` 以外の入力では `-q/--charge` が必須ですが、`--ligand-charge` がある場合は例外（PDB 入力、または `--ref-pdb` 付きXYZ/GJF）。明示的な `-q` は常に優先され、多重度は省略時に `1` がデフォルトです。
- ステージ結果（`result.xyz` と任意のPDB/GJFコンパニオン）は常に書き出され、軌跡は `--dump True` の場合のみ保存されます。PDB 入力では変換が有効な場合に `scan.pdb` が生成されます。

## YAML 設定（`--args-yaml`）
YAMLのルートはマッピングでなければなりません。YAML 値はCLIを上書きします。共有セクションは [YAML リファレンス](yaml-reference.md) の定義を再利用します。

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
  out_dir: ./result_scan/    # output directory
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
  out_dir: ./result_scan/    # output directory
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
  out_dir: ./result_scan/    # output directory
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
bias:
  k: 100                    # harmonic bias strength (eV·Å⁻²)
bond:
  device: cuda               # UMA device for bond analysis
  bond_factor: 1.2           # covalent-radius scaling
  margin_fraction: 0.05      # tolerance margin for comparisons
  delta_fraction: 0.05       # minimum relative change to flag bonds
```

---

## 関連項目

- [all](all.md) — 単一構造入力に `--scan-lists` を使用したエンドツーエンドワークフロー
- [path-search](path_search.md) — スキャン端点を中間体としてMEP探索
- [extract](extract.md) — スキャン前にポケットPDBを生成
- [YAML リファレンス](yaml-reference.md) — `bias` と `bond` の完全な設定オプション
- [用語集](glossary.md) — MEP、セグメントの定義