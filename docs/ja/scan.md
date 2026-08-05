# `scan`

調和拘束で結合長をスキャンして反応座標を駆動します。`pdb2reaction scan` は、単一構造から特定の原子間距離を駆動し、妥当な反応経路を探索する場面（`path-search` / `path-opt` の前処理として使われることが多い）で使用します。各stepで構造をL-BFGSまたはRFOptimizerで緩和します。mmCIF入力は内部PDBで計算し元IDのCIFを出力します。XYZ/GJF入力では`--ref-pdb`にPDBまたはmmCIF topologyを指定できます。

## 実行例

```bash
# 最小例: YAML spec から実行する
pdb2reaction scan -i input.pdb -q 0 -m 1 -s scan.yaml -o ./result_scan
```

```bash
# リテラル入力を使う
pdb2reaction scan -i input.pdb -q 0 -m 1 -s '[("TYR,285,CA","SAM,309,C10",1.35)]'
```

```bash
# ステージごとの軌跡を保存して確認する
pdb2reaction scan -i input.pdb -q 0 -m 1 -s scan.yaml --dump -o ./result_scan_dump
```

コマンド形式:

```bash
pdb2reaction scan -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m MULT] \
 [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] \
 [-s/--scan-lists scan.yaml | '[(i,j,targetÅ),...]'] [options] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE]
```

> **Note:** `-s/--scan-lists` の解釈結果を確認したい場合は `--print-parsed` を追加してください。

## 処理の流れ

1. `geom_loader` で構造を読み込み、電荷を解決します。電荷の解決順序の詳細は {ref}`CLI 規約: 電荷の指定 <ja-charge-specification>` を参照してください。
2. `--preopt` の場合、バイアスをかける前に無バイアスの前処理最適化を実行し、開始構造を緩和します。
3. `-s/--scan-lists`からstage targetを読み取り、indexを正規化します。3-field selector（例: `'TYR,285,CA'`）は順不同です。残基名や番号が重複するときは、位置固定の`CHAIN:RESNAME:RESSEQ[ICODE]:ATOM`を使います。
    各結合について変位 `Δ = target − current` を計算し、`h = --max-step-size` として `N = ceil(max(|Δ|) / h)` ステップに分割します。各結合は `δ = Δ / N` ずつ更新されます。
4. すべてのステップを順に進め、一時ターゲットを更新しながら調和ポテンシャル `E = Σ ½ k (|ri − rj| − target)²` を適用し、MLIP バックエンドで最適化します。最適化サイクルの上限は YAML `opt.max_cycles` から読み、明示した `--relax-max-cycles` がそれを上書きします。
5. 各ステージの最終ステップ後、必要に応じて無バイアス緩和（`--endopt`）を実行し、共有結合の変化を報告して `result.*` を出力します。
6. すべてのステージについて繰り返します。全ステージ連結の `scan_trj.xyz` は常に書き出されます。PDB/CIF/GJF companion には `--convert-files` と参照テンプレートが必要です。`--dump`（CLI フラグ）を指定すると、ステップごとの最適化軌跡ファイルも書き出されます（YAML の `opt.dump` は実行時スコープで上書きされるため無効です）。

## 出力

```
out_dir/ (デフォルト:./result_scan/)
├─ preopt/ # --preopt が True の場合
│ ├─ result.xyz
│ ├─ result.pdb # 参照トポロジーがあり変換有効時
│ ├─ result.cif # mmCIF/oversized-PDB bridge 入力時。元IDを復元
│ └─ result.gjf # Gaussian テンプレートがあり変換有効時
├─ stage_XX/ # ステージごとのフォルダ
│ ├─ result.xyz
│ ├─ result.pdb # 最終構造に対応する PDB（変換有効時）
│ ├─ result.cif # bridge入力。元IDを復元
│ ├─ result.gjf # テンプレートがある場合に対応する Gaussian（変換有効時）
│ ├─ scan_trj.xyz # 常に生成（連結バイアス軌跡）
│ ├─ scan.pdb # 参照トポロジーがあり変換有効時に生成（scan.gjf は生成されない）
│ └─ scan.cif # mmCIF/oversized-PDB bridge 入力時。元IDを復元した軌跡
├─ scan_trj.xyz # 全ステージ連結のスキャン軌跡
├─ scan.pdb # 全ステージ連結の PDB 軌跡（変換有効時）
└─ scan.cif # bridge入力のCIF軌跡
```

- `geom`/`calc`/`opt`/`bias`/`bond` および最適化ブロックの解決結果と、各ステージの結合変化レポートがコンソールに出力されます。

主な確認先:

- `result_scan/stage_01/result.pdb`（または `result.xyz`）
- `result_scan/stage_02/result.pdb`（または `result.xyz`）
- `result_scan/stage_*/scan_trj.xyz`（常に生成。参照トポロジー + 変換有効時は `scan.pdb`、bridge入力時は `scan.cif` も生成）

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 共通bridgeが受け入れる構造file（PDB / mmCIF / XYZ / GJF / TRJ） | 必須 |
| `-q, --charge INT` | 総電荷（CLI > テンプレート）。`-q` を省略して `--ligand-charge` がある場合は電荷が導出され、明示的な `-q` が最優先 | `.gjf` テンプレートまたは `--ligand-charge` がない場合は必須 |
| `-l, --ligand-charge TEXT` | 単一の整数（例: `-1`）でリガンド総電荷を指定するか、残基別マッピング（例: `GPP:-3,SAM:1`）で PDB/mmCIF 残基電荷から全系の電荷を導出。`-q` 省略時に使用（PDB/mmCIF 入力、または `--ref-pdb` 付き XYZ/GJF） | _None_ |
| `--workers`, `--workers-per-node` | UMA 予測器の並列度（`workers_per_node` は並列予測器へ転送）。`workers > 1` と明示的な解析 Hessian は併用不可。{ref}`ja-workers-analytical-error` を参照 | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度 2S+1。`.gjf` テンプレートがあれば継承し、未指定時は `1` | `.gjf` テンプレート値または `1` |
| `-s, --scan-lists TEXT` | YAML/JSONまたはinline literal。`i`/`j`は整数、3-field selector、または位置固定`CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` | 必須 |
| `--one-based/--zero-based` | 原子インデックスを 1 始まり/0 始まりとして解釈。これらは同一フラグの相互排他エイリアス（`--one-based` → `True`、`--zero-based` → `False`） | `True` |
| `--print-parsed/--no-print-parsed` | `-s/--scan-lists` 解釈後のステージ情報を表示 | `False` |
| `--max-step-size FLOAT` | 1 ステップあたりのスキャン結合の最大変化量（Å）。ステップ数を決定 | `0.20` |
| `--bias-k FLOAT` | 調和バイアス強度 `k`（eV·Å⁻²） | `300` |
| `--relax-max-cycles INT` | 前処理・各バイアスステップ・後処理における最適化サイクルの上限。明示値は YAML `opt.max_cycles` を上書き | `10000` |
| `--opt-mode TEXT` | `grad` → L-BFGS、`hess` → RFOptimizer。同じトークンが `tsopt` では Dimer / RS-P-RFO に対応する点については {ref}`ja-opt-mode-semantics` を参照してください | `grad` |
| `--freeze-links/--no-freeze-links` | PDB/mmCIF トポロジー入力時にキャップ水素の親原子を凍結 | `True` |
| `--freeze-atoms TEXT` | 凍結する原子の 1 始まりインデックスをカンマ区切りで明示的に指定（例: `'1,3,5'`）。`--freeze-links` と併用可、任意の入力形式に適用 | _None_ |
| `--dump/--no-dump` | ステップごとの最適化器軌跡ファイルを書き出します（`opt_cfg["dump"]` に転送）。`scan_trj.xyz` はこのフラグに関係なく常に書き出されます | `False` |
| `--convert-files/--no-convert-files` | XYZ/TRJ → PDB/CIF/GJF companionを切り替え | `True` |
| `--ref-pdb FILE` | XYZ/GJF入力の参照PDBまたはmmCIF topology | _None_ |
| `-o, --out-dir TEXT` | 出力ディレクトリ | `./result_scan/` |
| `--thresh TEXT` | 収束プリセットの上書き（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） | `gau` |
| `--config FILE` | ベース YAML 設定ファイル（最初に適用） | _None_ |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP バックエンド | `uma` |
| `--solvent TEXT` | xTB 暗黙溶媒（例: `water`）。`none` で無効化 | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB 溶媒モデル | `alpb` |
| `--preopt/--no-preopt` | スキャン前に無バイアス最適化を実行。**スコープ依存デフォルト:** 単体では `False`、`pdb2reaction all` 経由では `True` に反転されます（{ref}`all → スキャンオプション <ja-scan-options-single-input-runs>` を参照） | `False` |
| `--endopt/--no-endopt` | 各ステージ後に無バイアス最適化を実行 | `False` |
| `--out-json/--no-out-json` | `out_dir` に機械可読な `result.json` を書き出す。スキーマは [JSON 出力スキーマ](json-output.md) を参照 | `False` |

### 共有 YAML セクション
- `geom`, `calc`, `opt`, `lbfgs`, `rfo`: [YAML リファレンス](yaml-reference.md) と同じキーを使用します。`opt.dump` は実行時スコープで常に上書きされる（YAML では設定できない）ため、ステージ軌跡の出力は `--dump`（CLI）で制御します。
- 明示した `--relax-max-cycles` は YAML `opt.max_cycles` を上書きします。省略時は YAML が優先され、いずれもなければデフォルト `10000` です。

### セクション `bias`
- `k`（`300`）: 調和バイアス強度（eV·Å⁻²）。

(ja-section-bond)=
### セクション `bond`
`path-search` と共通の MLIP ベース結合変化検出:
- `device`（`"auto"`）: 結合解析用 MLIP デバイス。
- `bond_factor`（`1.20`）: 共有結合半径のスケーリング係数。
- `margin_fraction`（`0.05`）: 比較時の相対許容値。
- `delta_fraction`（`0.05`）: 結合の形成・切断を判定する最小相対変化量。

## YAML 設定

scan の output directory は command-owned の `-o/--out-dir` で指定します。YAML の optimizer `out_dir` は無視されます。

```yaml
geom:
 coord_type: cart # coordinate type: cartesian vs dlc internals
 freeze_atoms: [] # 1-based frozen atoms merged with CLI/link detection
calc:
 charge: 0 # total charge (CLI/template override)
 spin: 1 # spin multiplicity 2S+1
 model: uma-s-1p2 # uma-s-1p2 | uma-m-1p1
 task_name: omol # UMA task name
 device: auto # MLIP device selection
 max_neigh: null # maximum neighbors for graph construction
 radius: null # cutoff radius for neighbor search
 r_edges: false # store radial edges
 out_hess_torch: true # request torch-form Hessian
 freeze_atoms: null # calculator-level frozen atoms
 hessian_calc_mode: FiniteDifference # Hessian mode selection
 return_partial_hessian: true  # partial Hessian over active DOFs
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
 energy_plateau: true # enable plateau-based early stop
 energy_plateau_thresh: 1.0e-04 # plateau detection threshold
 energy_plateau_window: 50 # plateau detection window (cycles)
 line_search: true # enable line search
 dump: false # dump trajectory/restart data
 dump_restart: false # dump restart checkpoints
 prefix: "" # filename prefix
lbfgs:
 thresh: gau # L-BFGS convergence preset
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
 energy_plateau: true # enable plateau-based early stop
 energy_plateau_thresh: 1.0e-04 # plateau detection threshold
 energy_plateau_window: 50 # plateau detection window (cycles)
 line_search: true # enable line search
 dump: false # dump trajectory/restart data
 dump_restart: false # dump restart checkpoints
 prefix: "" # filename prefix
 keep_last: 7 # history size for L-BFGS buffers
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
 energy_plateau: true # enable plateau-based early stop
 energy_plateau_thresh: 1.0e-04 # plateau detection threshold
 energy_plateau_window: 50 # plateau detection window (cycles)
 line_search: true # enable line search
 dump: false # dump trajectory/restart data
 dump_restart: false # dump restart checkpoints
 prefix: "" # filename prefix
 trust_radius: 0.10 # trust-region radius
 trust_update: true # enable trust-region updates
 trust_min: 0.0001 # minimum trust radius
 trust_max: 0.10 # maximum trust radius
 max_energy_incr: null # allowed energy increase per step
 hessian_update: bfgs # Hessian update scheme
 hessian_init: calc # Hessian initialization source
 hessian_recalc: 500 # rebuild Hessian every N steps
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
 device: auto # MLIP device for bond analysis
 bond_factor: 1.2 # covalent-radius scaling
 margin_fraction: 0.05 # tolerance margin for comparisons
 delta_fraction: 0.05 # minimum relative change to flag bonds
```

`opt`/`lbfgs`/`rfo`/`bias`/`bond` のさらなる YAML オプションとそのデフォルトは [YAML リファレンス](yaml-reference.md) を参照してください。

## スキャンリスト仕様

YAML/JSON ファイル書式、インライン Python リテラル構文、原子セレクタ、クォート規則については
{ref}`CLI 規約: スキャンリスト仕様 <ja-scan-list-spec>` を参照してください。

### 段階的スキャン vs 協奏的スキャン

1 つのリテラル内の `(i, j, target)` タプルの数とリテラルの数が、座標を*同時に*駆動する（協奏的）か*順次*駆動する（段階的）かを決めます。

| モード | 構文 | 使う場面 |
| --- | --- | --- |
| 協奏的 | 1 つの `-s` に複数の座標タプル | 座標が 1 ステップで一緒に動く。機構をステージに分解する必要がない |
| 段階的 | `-s` を順次ステージごとに 1 リテラルずつ | 機構が事前に分かっており、ステップごとの制御とステージ別出力がほしい |

段階的形式では機構を事前に定義し、指定した各座標変化を別々のステージとして扱います。機構が未知あるいは多段階の場合は、ステージを自分で推測せず [`path-search`](path-search.md) に経路を自動分割させます。（4-tuple `(i, j, low, high)` は双方向の 2 ステージスキャンに展開されます。[双方向スキャン](#双方向スキャン4-tuple)を参照。）

```bash
# 協奏的: 2 つの座標を 1 ステージで一緒に駆動
pdb2reaction scan -i reactant.pdb \
    -q 0 -m 1 \
    -s '[("Ca RES 10","Cb RES 11",1.6),("H RES 11","O GLU 20",1.0)]' -o result_concerted
```

段階的スキャンでは 1 つの `-s/--scan-lists` フラグの後に複数のリテラルを並べます。各リテラルが 1 ステージになります。

```bash
# ステージ 1: 1 つの結合を 1.35 Å に駆動
# ステージ 2: 2 つの結合を同時に駆動
-s \
 '[("TYR,285,CA","SAM,309,C10",1.35)]' \
 '[("TYR,285,CA","SAM,309,C10",2.20),("TYR,285,CB","SAM,309,C11",1.80)]'
```

ステージは順次実行され、各ステージは前ステージの緩和結果から開始します。

(ja-scan-direction-barrier-sign)=
### スキャン方向とバリアの符号

`scan` の結果は sampled/final energy を保存しますが、barrier field を認定・出力しません。`scan`（または経路）が**生成物側**から始まり、利用者がそこから障壁を算出する場合、その差は**逆方向**のバリア `E(TS) − E(product)` です。順方向のバリアを引用するには反応物から計算します。

| 実行内容 | 順方向バリア |
| --- | --- |
| 生成物始点のスキャン | `E(TS) − E(reactant)` — 生成物基準の差では**ない** |

これはフラグではなく読み取り時の解釈です。特に結晶構造の生成物複合体から開始した場合は、バリアを引用する前にスキャンがどちらの端点から始まったかを必ず確認してください。

### 双方向スキャン（4-tuple）

3-tuple `(i, j, target)` の代わりに **4-tuple** `(i, j, start, end)` を指定すると、現在の構造から両方向にスキャンします。CLI は各 4-tuple を自動的に 2 ステージに展開します:

1. **パス 1:** `i`--`j` の距離を現在の値から `start` に向けて駆動。
2. **パス 2:** 初期構造を復元し、`i`--`j` の距離を `end` に向けて駆動。

連結軌跡は `start → 初期構造 → end` の順に組み立てられ、出発構造を通る連続的な経路が得られます。

```bash
# 双方向スキャン: 結合 12--45 を現在の構造から
# 1.35 Å（パス 1）と 2.50 Å（パス 2）に向けて駆動
pdb2reaction scan -i input.pdb -q 0 -s '[(12, 45, 1.35, 2.50)]'
```

これは 2 つの手動ステージの間にジオメトリリセットを行うのと同等ですが、スクリプトを書く必要がありません。同じリテラル内で 3-tuple と 4-tuple を混在させることもできます。

```{note}
**4-tuple 使用時のステージ番号。** 1 つの 4-tuple は出力ツリー内で **2 つ** のステージに展開されます。`start` パスは `stage_NN/` に、`end` パスは `stage_NN+1/` に書き込まれます。したがって最初のリテラルとして 1 個の 4-tuple を渡した場合、1 つの統合された `stage_01/` ではなく `stage_01/` と `stage_02/` が作成されます。3-tuple と 4-tuple を混在させた場合、カウンターは 3-tuple ごとに `+1`、4-tuple ごとに `+2` 進みます。
```

## 注意事項

- スキャンの入力は 1 つの構造 + `-s/--scan-lists scan.yaml`（推奨）または `-s/--scan-lists` の 1 個以上のインラインリテラル（1 リテラル = 1 ステージ）です。YAML/JSON ファイルパスはシェルのクォート問題を避けられ、バージョン管理にも向きます。インライン Python リテラルは単純な単一ステージのスキャンには十分です。
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。
- `-s/--scan-lists` には単一フラグの後に複数リテラルを並べます。ターゲット距離は正の値である必要があります。原子インデックスは内部で 0 始まりに正規化されます。PDB/mmCIF トポロジーではセレクタ文字列を使用できます。3フィールドの従来形は順不同ですが、chainまで指定する4フィールド形は `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` の位置固定です。繰り返し残基では4フィールド形を使用してください。
- `--freeze-links` が有効な場合、キャップ水素の親原子は自動的に凍結されます（{ref}`キャップ水素と凍結原子 <ja-link-hydrogen-and-frozen-atoms>` を参照）。
- ステージ結果の `result.xyz` と全ステージ連結の `scan_trj.xyz` は常に書き出されます。対応する PDB/CIF/GJF companion は `--convert-files` が有効で参照トポロジーを利用できる場合に生成されます。`--dump`（CLI）を指定すると、最適化器によるステップごとのダンプが有効になります（YAML の `opt.dump` は実行時スコープで上書きされるため無効です）。

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け

- [all](all.md) — 単一構造入力に `--scan-lists` を使用した一気通貫ワークフロー
- [scan2d](scan2d.md) — 同じ MLIP backend と YAML 設定で 2 距離（d₁, d₂）グリッドスキャン
- [scan3d](scan3d.md) — 3 距離（d₁, d₂, d₃）グリッドスキャン + 等値面出力
- [path-search](path-search.md) — スキャン端点を中間体として MEP を探索
- [extract](extract.md) — スキャン前に活性部位モデル（バインディングポケット） PDB を生成
- [YAML リファレンス](yaml-reference.md) — `bias` と `bond` の完全な設定オプション
- [用語集](glossary.md) — MEP、セグメントの定義
