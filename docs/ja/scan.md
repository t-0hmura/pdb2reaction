# `scan`

## 概要

> **要約:** 調和拘束を用いて結合距離をスキャンし、反応座標を駆動します。`--spec`（YAML/JSON、推奨）でターゲット距離を指定し、`--scan-lists` は 入力として利用できます。

### 要点
- **想定場面:** 単一構造から特定の原子間距離を変化させ、もっともらしい反応経路を探索したい場合に使います（`path-search` / `path-opt` の前処理として使うことが多い）。
- **入力:** 1 つの構造 + `--spec scan.yaml`（推奨）または `--scan-lists` の 1 個以上の リテラル（**1 リテラル = 1 ステージ**）。
- **既定値:** `--opt-mode light`（LBFGS）、`--preopt`、`--endopt`、`--max-step-size 0.20 Å`。
- **主な出力:** ステージごとの `result.xyz`（必要に応じて `.pdb`/`.gjf`）。`--dump` なら結合した軌跡も保存。
- **注意:** 可能な限り `--spec` を使ってください。`--scan-lists` は **Python リテラル**のためクォート/エスケープが必要です。

`pdb2reaction scan` は、UMA 計算機と調和拘束を使った段階的な結合長スキャンを実行します。各ステップで一時ターゲットを更新し、拘束ポテンシャルを適用したうえで、構造全体を LBFGS（`--opt-mode light`）または RFOptimizer（`--opt-mode heavy`）で緩和します。


XYZ/GJF 入力では、`--ref-pdb` で参照 PDB トポロジーを指定すると、XYZ 座標を保持したまま PDB/GJF へのフォーマット対応変換が可能になります。

## 最小例

```bash
pdb2reaction scan -i input.pdb -q 0 -m 1 --spec scan.yaml --out-dir ./result_scan
```

## 出力の見方

- `result_scan/stage_01/result.pdb`（または `result.xyz`）
- `result_scan/stage_02/result.pdb`（または `result.xyz`）
- `--dump` 指定時は `result_scan/stage_*/scan_trj.xyz` と `scan.pdb`

## よくある例

1. YAML spec から実行する。

```bash
pdb2reaction scan -i input.pdb -q 0 -m 1 --spec scan.yaml
```

2. リテラル入力を使う。

```bash
pdb2reaction scan -i input.pdb -q 0 -m 1 --scan-lists '[("TYR,285,CA","MMT,309,C10",1.35)]'
```

3. ステージごとの軌跡を保存して確認する。

```bash
pdb2reaction scan -i input.pdb -q 0 -m 1 --spec scan.yaml --dump --out-dir ./result_scan_dump
```

> **Note:** `--spec` / `--scan-lists` の解釈結果を確認したい場合は `--print-parsed` を追加してください。

## 使用法
```bash
pdb2reaction scan -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [-m MULT] \
 [--spec scan.yaml | --scan-lists '[(i,j,targetÅ),...]'] [options] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE]
```

### 例
```bash
# 単一ステージの最小例
pdb2reaction scan -i input.pdb -q 0 --scan-lists '[("TYR,285,CA","MMT,309,C10",1.35)]'

# 2ステージ、LBFGS緩和、軌跡ダンプ
pdb2reaction scan -i input.pdb -q 0 --scan-lists \
 '[("TYR,285,CA","MMT,309,C10",1.35)]' \
 '[("TYR,285,CA","MMT,309,C10",2.20),("TYR,285,CB","MMT,309,C11",1.80)]' \
 --max-step-size 0.20 --dump --out-dir ./result_scan/ --opt-mode light \
 --preopt --endopt

# 単一の --scan-lists の後に複数リテラルを渡す
pdb2reaction scan -i input.pdb -q 0 --scan-lists \
 '[("TYR,285,CA","MMT,309,C10",1.35)]' \
 '[("TYR,285,CA","MMT,309,C10",2.20),("TYR,285,CB","MMT,309,C11",1.80)]'
```

## `--scan-lists` の書式

`--scan-lists` は **Python リテラル文字列**として CLI に渡します。シェルのクォート処理に注意が必要です。

### 基本構造

各リテラルは、三つ組 `(原子1, 原子2, ターゲット距離Å)` の Python リストです。

```
--scan-lists '[(原子1, 原子2, ターゲット距離Å),...]'
```

- リスト全体を **シングルクォート** `'...'` で囲みます（シェルがカッコや空白を解釈しないようにするため）。
- 各三つ組は `原子1`–`原子2` 間の距離を `ターゲット距離Å` まで駆動します。
- **1 リテラル = 1 ステージ**です。複数ステージを実行するには、**1 つの `--scan-lists` フラグの後に複数リテラル**を並べます（フラグの繰り返しは不可）。

### 原子の指定方法

原子は**整数インデックス**または **PDB セレクタ文字列**で指定します。

| 方法 | 例 | 備考 |
| --- | --- | --- |
| 整数インデックス | `(1, 5, 2.0)` | デフォルトは 1 始まり（`--one-based`） |
| PDB セレクタ | `("TYR,285,CA", "MMT,309,C10", 2.0)` | 残基名、残基番号、原子名の三つ組 |

PDB セレクタのトークンは、カンマ `,`、スペース、スラッシュ `/`、バッククォート `` ` ``、バックスラッシュ `\` のいずれでも区切れます。トークンの順序も任意です。

```bash
# 以下はすべて同じ原子を指定:
"TYR,285,CA"
"TYR 285 CA"
"TYR/285/CA"
"285,TYR,CA" # 順序は自由
```

### クォートの規則

```bash
# 正しい: シングルクォートでリストを囲み、セレクタはダブルクォート
--scan-lists '[("TYR,285,CA","MMT,309,C10",1.35)]'

# 正しい: 整数インデックスなら内側のダブルクォートは不要
--scan-lists '[(1, 5, 2.0)]'

# 非推奨: ダブルクォートで外側を囲むとエスケープが必要
--scan-lists "[(\"TYR,285,CA\",\"MMT,309,C10\",1.35)]"
```

### 複数ステージ

1 つの `--scan-lists` フラグの後に、複数のリテラルを並べます。各リテラルが 1 ステージになります。

```bash
# ステージ 1: 1 つの結合を 1.35 Å に駆動
# ステージ 2: 2 つの結合を同時に駆動
--scan-lists \
 '[("TYR,285,CA","MMT,309,C10",1.35)]' \
 '[("TYR,285,CA","MMT,309,C10",2.20),("TYR,285,CB","MMT,309,C11",1.80)]'
```

ステージは順次実行され、各ステージは前ステージの緩和結果から開始します。**`--scan-lists` フラグは繰り返さず**、すべてのリテラルを 1 つのフラグの後に記述してください。

## ワークフロー
1. `geom_loader` で構造を読み込み、CLI の上書き値・埋め込み Gaussian テンプレート（存在する場合）・デフォルト値から電荷とスピンを解決します。`-q` が省略され `--ligand-charge` が与えられている場合は、入力を酵素--基質複合体として扱います。PDB 入力（または `--ref-pdb` 付き XYZ/GJF）では `extract.py` の電荷サマリーから総電荷を導出します。
2. `--preopt` の場合、バイアスをかける前に無バイアスの前処理最適化を実行し、開始構造を緩和します。
3. `--scan-lists` で与えられた各ステージリテラルについて `(i, j)` を解析・正規化します（デフォルトは 1 始まり）。PDB 入力では、各エントリに整数インデックスまたは `'TYR,285,CA'` のような原子セレクタ文字列を指定できます。セレクタの区切りは空白・カンマ・スラッシュ・バッククォート・バックスラッシュのいずれも可で、トークン順序は任意です（フォールバックは resname, resseq, atom を想定）。
 各結合について変位 `Δ = target − current` を計算し、`h = --max-step-size` として `N = ceil(max(|Δ|) / h)` ステップに分割します。各結合は `δ = Δ / N` ずつ更新されます。
4. すべてのステップを順に進め、一時ターゲットを更新しながら調和ポテンシャル `E = Σ ½ k (|ri − rj| − target)²` を適用し、UMA で最小化します。最適化サイクルの上限は `--relax-max-cycles` で設定します（YAML で `opt.max_cycles` が指定されていない場合）。
5. 各ステージの最終ステップ後、必要に応じて無バイアス緩和（`--endopt`）を実行し、共有結合の変化を報告して `result.*` を出力します。
6. すべてのステージについて繰り返します。軌跡は `--dump` の場合のみ保存されます。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル | 必須 |
| `-q, --charge INT` | 総電荷（CLI > テンプレート）。`-q` を省略して `--ligand-charge` がある場合は電荷が導出され、明示的な `-q` が最優先 | `.gjf` テンプレートまたは `--ligand-charge` がない場合は必須 |
| `--ligand-charge TEXT` | 残基別電荷マッピング（例: `GPP:-3,SAM:1`）。PDB の残基電荷から全系の電荷を自動導出します（手動計算不要）。`-q` 省略時に使用（PDB 入力、または `--ref-pdb` 付き XYZ/GJF） | _None_ |
| `--workers`, `--workers-per-node` | UMA 予測器の並列度（workers > 1 で解析ヘシアンは無効化; `workers_per_node` は並列予測器へ転送） | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度 2S+1。`.gjf` テンプレートがあれば継承し、未指定時は `1` | `.gjf` テンプレート値または `1` |
| `--spec FILE` | `stages` を持つ YAML/JSON スキャン仕様。`one_based` を任意指定可能。 | 推奨 |
| `--scan-lists TEXT` | : `(i,j,targetÅ)` タプルを含む Python リテラル。各リテラルが 1 ステージ; 1 つのフラグの後に複数リテラルを渡す。`i`/`j` は整数インデックスまたは PDB 原子セレクタ（`'TYR,285,CA'`） | `--spec` の代替 |
| `--one-based/--zero-based` | 原子インデックスを 1 始まり/0 始まりとして解釈 | `True` |
| `--print-parsed/--no-print-parsed` | `--spec`/`--scan-lists` 解釈後のステージ情報を表示。 | `False` |
| `--max-step-size FLOAT` | 1 ステップあたりのスキャン結合の最大変化量（Å）。ステップ数を決定 | `0.20` |
| `--bias-k FLOAT` | 調和バイアス強度 `k`（eV·Å⁻²） | `300` |
| `--relax-max-cycles INT` | 前処理・各バイアスステップ・後処理における最適化サイクルの上限。YAML で `opt.max_cycles` が指定されていない場合に使用 | `10000` |
| `--opt-mode TEXT` | `light` → LBFGS、`heavy` → RFOptimizer | `light` |
| `--freeze-links/--no-freeze-links` | PDB 入力時にリンク水素の親原子を凍結 | `True` |
| `--dump/--no-dump` | バイアス付き軌跡（`scan_trj.xyz`/`scan.pdb`）を出力 | `False` |
| `--convert-files/--no-convert-files` | PDB/Gaussian 入力で XYZ/TRJ → PDB/GJF コンパニオン変換を切り替え（軌跡変換は PDB のみ） | `True` |
| `--ref-pdb FILE` | XYZ/GJF 入力時の参照 PDB トポロジー（XYZ 座標は保持） | _None_ |
| `--out-dir TEXT` | 出力ディレクトリ | `./result_scan/` |
| `--thresh TEXT` | 収束プリセットの上書き（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） | `gau` |
| `--config FILE` | ベース YAML 設定ファイル（最初に適用） | _None_ |
| `--preopt/--no-preopt` | スキャン前に無バイアス最適化を実行 | `True` |
| `--endopt/--no-endopt` | 各ステージ後に無バイアス最適化を実行 | `True` |

### 共有 YAML セクション
- `geom`, `calc`, `opt`, `lbfgs`, `rfo`: [YAML リファレンス](yaml_reference.md) と同じキーを使用します。`opt.dump` は YAML で設定可能ですが、ステージ軌跡の出力は `--dump` で制御します。
- `--relax-max-cycles` は**明示的に指定され**、かつ YAML で `opt.max_cycles` が設定されていない場合にのみ適用されます（デフォルト `10000`）。

### セクション `bias`
- `k`（`300`）: 調和バイアス強度（eV·Å⁻²）。

(section-bond)=
### セクション `bond`
`path-search` と共通の UMA ベース結合変化検出:
- `device`（`"cuda"`）: 結合解析用 UMA デバイス。
- `bond_factor`（`1.20`）: 共有結合半径のスケーリング係数。
- `margin_fraction`（`0.05`）: 比較時の相対許容値。
- `delta_fraction`（`0.05`）: 結合の形成・切断を判定する最小相対変化量。

## 出力
```
out_dir/ (デフォルト:./result_scan/)
├─ preopt/ # --preopt が True の場合
│ ├─ result.xyz
│ ├─ result.pdb # PDB 入力かつ変換有効時
│ └─ result.gjf # Gaussian テンプレートがあり変換有効時
└─ stage_XX/ # ステージごとのフォルダ
 ├─ result.xyz
 ├─ result.pdb # 最終構造の PDB ミラー（変換有効時）
 ├─ result.gjf # テンプレートがある場合の Gaussian ミラー（変換有効時）
 ├─ scan_trj.xyz # --dump の場合
 └─ scan.pdb # PDB 入力で変換有効時の軌跡コンパニオン（scan.gjf は生成されない）
```
- `geom`/`calc`/`opt`/`bias`/`bond` および最適化ブロックの解決結果と、各ステージの結合変化レポートがコンソールに出力されます。

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes_common_errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- `--scan-lists` には単一フラグの後に複数リテラルを並べてください。フラグの繰り返しには対応していません。ターゲット距離は正の値である必要があります。原子インデックスは内部で 0 始まりに正規化されます。PDB 入力ではセレクタ文字列を使用でき、空白・カンマ・スラッシュ・バッククォート・バックスラッシュで区切れます。トークン順序は任意です。
- `--freeze-links` は、ユーザー指定の `freeze_atoms` に PDB のリンク水素親原子を追加し、ポケット構造を固定します。
- ステージ結果（`result.xyz` と任意の PDB/GJF コンパニオン）は `--dump` の設定にかかわらず常に書き出されます。軌跡は `--dump` の場合のみ保存され、PDB 入力かつ変換が有効な場合は `scan.pdb` も生成されます。


```yaml
geom:
 coord_type: cart # coordinate type: cartesian vs dlc internals
 freeze_atoms: [] # 0-based frozen atoms merged with CLI/link detection
calc:
 charge: 0 # total charge (CLI/template override)
 spin: 1 # spin multiplicity 2S+1
 model: uma-s-1p1 # UMA model tag
 task_name: omol # UMA task name
 device: auto # UMA device selection
 max_neigh: null # maximum neighbors for graph construction
 radius: null # cutoff radius for neighbor search
 r_edges: false # store radial edges
 out_hess_torch: true # request torch-form Hessian
 freeze_atoms: null # calculator-level frozen atoms
 hessian_calc_mode: FiniteDifference # Hessian mode selection
 return_partial_hessian: false # full Hessian (avoids shape mismatches)
opt:
 thresh: gau # convergence preset (Gaussian/Baker-style)
 max_cycles: 10000 # optimizer cycle cap
 print_every: 100 # logging stride
 min_step_norm: 1.0e-08 # minimum norm for step acceptance
 assert_min_step: true # stop if steps fall below threshold
 rms_force: null # explicit RMS force target
 rms_force_only: false # rely only on RMS force convergence
 max_force_only: false # rely only on max force convergence
 force_only: false # skip displacement checks
 converge_to_geom_rms_thresh: 0.05 # geom RMS threshold when converging to ref
 overachieve_factor: 0.0 # factor to tighten thresholds
 check_eigval_structure: false # validate Hessian eigenstructure
 line_search: true # enable line search
 dump: false # dump trajectory/restart data
 dump_restart: false # dump restart checkpoints
 prefix: "" # filename prefix
 out_dir:./result_scan/ # output directory
lbfgs:
 thresh: gau # LBFGS convergence preset
 max_cycles: 10000 # iteration limit
 print_every: 100 # logging stride
 min_step_norm: 1.0e-08 # minimum accepted step norm
 assert_min_step: true # assert when steps stagnate
 rms_force: null # explicit RMS force target
 rms_force_only: false # rely only on RMS force convergence
 max_force_only: false # rely only on max force convergence
 force_only: false # skip displacement checks
 converge_to_geom_rms_thresh: 0.05 # RMS threshold when targeting geometry
 overachieve_factor: 0.0 # tighten thresholds
 check_eigval_structure: false # validate Hessian eigenstructure
 line_search: true # enable line search
 dump: false # dump trajectory/restart data
 dump_restart: false # dump restart checkpoints
 prefix: "" # filename prefix
 out_dir:./result_scan/ # output directory
 keep_last: 7 # history size for LBFGS buffers
 beta: 1.0 # initial damping beta
 gamma_mult: false # multiplicative gamma update toggle
 max_step: 0.3 # maximum step length
 control_step: true # control step length adaptively
 double_damp: true # double damping safeguard
 mu_reg: null # regularization strength
 max_mu_reg_adaptions: 10 # cap on mu adaptations
rfo:
 thresh: gau # RFOptimizer convergence preset
 max_cycles: 10000 # iteration cap
 print_every: 100 # logging stride
 min_step_norm: 1.0e-08 # minimum accepted step norm
 assert_min_step: true # assert when steps stagnate
 rms_force: null # explicit RMS force target
 rms_force_only: false # rely only on RMS force convergence
 max_force_only: false # rely only on max force convergence
 force_only: false # skip displacement checks
 converge_to_geom_rms_thresh: 0.05 # RMS threshold when targeting geometry
 overachieve_factor: 0.0 # tighten thresholds
 check_eigval_structure: false # validate Hessian eigenstructure
 line_search: true # enable line search
 dump: false # dump trajectory/restart data
 dump_restart: false # dump restart checkpoints
 prefix: "" # filename prefix
 out_dir:./result_scan/ # output directory
 trust_radius: 0.1 # trust-region radius
 trust_update: true # enable trust-region updates
 trust_min: 0.0 # minimum trust radius
 trust_max: 0.1 # maximum trust radius
 max_energy_incr: null # allowed energy increase per step
 hessian_update: bfgs # Hessian update scheme
 hessian_init: calc # Hessian initialization source
 hessian_recalc: 200 # rebuild Hessian every N steps
 hessian_recalc_adapt: null # adaptive Hessian rebuild factor
 small_eigval_thresh: 1.0e-08 # eigenvalue threshold for stability
 alpha0: 1.0 # initial micro step
 max_micro_cycles: 50 # micro-iteration limit
 rfo_overlaps: false # enable RFO overlaps
 gediis: false # enable GEDIIS
 gdiis: true # enable GDIIS
 gdiis_thresh: 0.0025 # GDIIS acceptance threshold
 gediis_thresh: 0.01 # GEDIIS acceptance threshold
 gdiis_test_direction: true # test descent direction before DIIS
 adapt_step_func: true # adaptive step scaling toggle
bias:
 k: 300 # harmonic bias strength (eV·Å⁻²)
bond:
 device: cuda # UMA device for bond analysis
 bond_factor: 1.2 # covalent-radius scaling
 margin_fraction: 0.05 # tolerance margin for comparisons
 delta_fraction: 0.05 # minimum relative change to flag bonds
```

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け

- [all](all.md) — 単一構造入力に `--scan-lists` を使用したend-to-endワークフロー
- [path-search](path_search.md) — スキャン端点を中間体として MEP を探索
- [extract](extract.md) — スキャン前にポケット PDB を生成
- [YAML リファレンス](yaml_reference.md) — `bias` と `bond` の完全な設定オプション
- [用語集](glossary.md) — MEP、セグメントの定義
