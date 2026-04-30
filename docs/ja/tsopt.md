# `tsopt`

## 概要

> **要約:** 遷移状態（TS）*候補*を、RS-I-RFO（Restricted-Step Image Rational Function Optimization）（`--opt-mode hess`、デフォルト）で最適化します。RS-I-RFO で収束しない場合や、完全 Hessian の再計算コストが過大な場合の **代替** TS optimizer として Hessian-Guided Dimer（`--opt-mode grad`）が利用可能です。`tsopt` は最後に自動で Hessian 計算と虚振動数チェックを行います。妥当な TS（一次鞍点）では虚振動数が **ちょうど 1 つ** です。端点の接続性は `irc` で確認してください。

### 要点
- **想定場面:** TS 初期構造（`path-opt` / `path-search` の HEI、または自前のもの）を一次鞍点に最適化し、虚振動数チェックを内蔵で行いたい場合。
- **手法:** `--opt-mode hess`（`rsirfo`）= RS‑I‑RFO（デフォルト、一般的により堅牢）。`--opt-mode grad`（`dimer`）= Hessian-Guided Dimer（RS-I-RFO で収束しない場合や完全 Hessian の再計算コストが過大な場合の代替）。`--flatten`（デフォルト無効）で余分な虚振動数モードの除去を制御します。
- **主な出力:** `final_geometry.{xyz,pdb,gjf}`、`vib/imag_*_trj.xyz`（PDB 入力では `.pdb` も）。`--dump` を指定した場合は最適化軌跡も。
- **デフォルト:** `--opt-mode hess`（RS-I-RFO）、`--thresh baker`、`--hessian-calc-mode FiniteDifference`、`--max-cycles 10000`、`--flatten` 無効、バックエンド `uma`。
- **次ステップ:** `tsopt` は最終 Hessian 計算と虚振動数チェックを内部で実行しますが、結果はなお *候補* です。[irc](irc.md) で端点の接続性を確認してください。完全な振動解析や熱化学補正が必要な場合のみ、別途 [freq](freq.md) を実行します。

> **命名規則の注意:** CLI は `grad|dimer`（= Dimer）および `hess|rsirfo`（= RS-I-RFO、デフォルト）を受け付けます。YAML ではトップレベルの `hessian_dimer:`（Dimer）または `rsirfo:`（RS-I-RFO）ブロックを直接指定してください。

XYZ/GJF 入力では `--ref-pdb` で参照 PDB トポロジーを指定し、XYZ 座標を保持したまま PDB/GJF への変換が可能です。TS 初期構造が必要な場合は、2 端点なら [path-opt](path-opt.md)、2 構造以上なら [path-search](path-search.md) で HEI を取得してから `tsopt`（内部で虚振動数チェック済み）→ `irc` の順で検証してください。

## 最小例

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --out-dir ./result_tsopt
```

## 出力の見方

- `result_tsopt/final_geometry.pdb`（または `final_geometry.xyz`）
- `result_tsopt/vib/imag_*_trj.xyz`
- `result_tsopt/vib/imag_*.pdb`（PDB 入力の場合）

## よくある例

1. VRAM に余裕がある場合に dimer モード + 解析的ヘシアンで実行する。

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
 --opt-mode grad --hessian-calc-mode Analytical --out-dir ./result_tsopt_grad
```

2. 最適化軌跡を保存して確認する。

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --dump --out-dir ./result_tsopt_dump
```

3. rsirfo モードを YAML 上書きと併用する。

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
 --opt-mode hess --config tsopt.yaml --out-dir ./result_tsopt_hess
```

4. rsirfo モードで flatten を有効化して実行する。

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
 --opt-mode hess --flatten --out-dir ./result_tsopt_flatten
```

## 使用法
```bash
pdb2reaction tsopt -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m 2S+1] \
 [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] \
 [--opt-mode grad|hess|dimer|rsirfo] [--flatten/--no-flatten] \
 [--freeze-links/--no-freeze-links] [--max-cycles N] [--thresh PRESET] \
 [--hessian-calc-mode Analytical|FiniteDifference] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE]
```

### 例
```bash
# 推奨ベースライン: 電荷/多重度を指定し rsirfo を使う
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --opt-mode hess --out-dir ./result_tsopt/

# YAML 上書き、有限差分ヘシアン、freeze-links処理を伴う dimer モード
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --freeze-links \
 --opt-mode grad --max-cycles 10000 --no-dump \
 --out-dir ./result_tsopt/ --config ./args.yaml \
 --hessian-calc-mode FiniteDifference

# YAMLのみで駆動される rsirfo モード（RS-I-RFO）
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --opt-mode hess \
 --config ./args.yaml --out-dir ./result_tsopt/
```

## ワークフロー
- **電荷/スピン解決**: 電荷の解決順序の詳細は {ref}`CLI 規約: 電荷の指定 <ja-charge-specification>` を参照してください。
- **構造ロードと freeze-links**: 構造は `pysisyphus.helpers.geom_loader` で読み込まれます。`--freeze-links` が有効な場合、リンク水素の親原子は自動的に凍結されます（{ref}`リンク水素と凍結原子 <ja-link-hydrogen-and-frozen-atoms>` を参照）。
- **MLIP ヘシアン（デフォルト: UMA）**: `--hessian-calc-mode` で解析的ヘシアンと有限差分ヘシアンを切り替えます。いずれも活性（PHVA）部分空間を考慮します。凍結原子が存在する場合、MLIP バックエンドは活性ブロックのみを返すことがあります。ヘシアン評価モードの詳細は {ref}`ja-hessian-evaluation` を参照してください。
- **Dimerモード詳細**:
 - Hessian Guided Dimer段階は、正確ヘシアン（活性サブスペース、TR射影）を周期的に評価してダイマー方向を更新します。`root == 0` のときは最小固有対に `torch.lobpcg` を優先し、失敗時は `torch.linalg.eigh` にフォールバックします。
 - `--flatten` が有効な場合、フラット化ループはΔxとΔgを用い、Bofill（SR1/MS ↔ PSBブレンド; `hessian_dimer.flatten_loop_bofill` で切替）で活性ヘシアンを更新します。各ループは虚振動数モード推定 → 1回フラット化 → ダイマー方向再更新 → dimer+LBFGSマイクロ区間 → （任意で）Bofill更新を実行します。虚振動数モードが1つになったら最終的な正確ヘシアンで振動解析を行います。
 - `root != 0` の場合は初期ダイマー方向のみそのrootを使用し、以降の更新は最も負のモード（`root = 0`）に従います。
- **RS-I-RFOモード**: RS-I-RFOを実行し、任意のヘシアン参照やR+S分割セーフガード、マイクロサイクル制御は `rsirfo` セクションで設定します。`--flatten` が有効で収束後も虚振動数モードが複数残る場合、追加モードをフラット化してRS-I-RFOを再実行し、虚振動数モードが1つになるか上限に達するまで繰り返します。
- **モード出力と変換**: 検出された虚振動数モードはすべて `vib/imag_*_trj.xyz` に書き出されます。PDB 入力で変換が有効な場合は `.pdb` としても出力されます。最適化軌跡と最終構造は `--dump` 時に入力テンプレート経由で PDB に変換されます。Gaussian テンプレートでは最終構造のみ `.gjf` が生成されます。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル | 必須 |
| `-q, --charge INT` | 総電荷。`.gjf` テンプレートまたは `--ligand-charge`（PDB 入力または `--ref-pdb` 付きXYZ/GJF）が提供しない限り必須。両方指定時は `-q` が優先 | テンプレート/導出が適用されない限り必須 |
| `-l, --ligand-charge TEXT` | 残基別電荷マッピング（例: `GPP:-3,SAM:1`）。PDB の残基電荷から全系の電荷を自動導出します（手動計算不要）。`-q` 省略時に使用（PDB 入力、または `--ref-pdb` 付き XYZ/GJF） | _None_ |
| `--workers INT` | MLIP予測器の並列度（workers > 1 で解析ヘシアン無効）。診断上の注意は {ref}`ja-workers-fd-downgrade` を参照。 | `1` |
| `--workers-per-node INT` | ノードあたりのワーカー数。並列予測器に渡されます | `1` |
| `-m, --multiplicity INT` | スピン多重度（2S+1） | `.gjf` テンプレート値または `1` |
| `--freeze-links/--no-freeze-links` | PDBのみ。リンク水素の親を凍結（`geom.freeze_atoms` にマージ） | `True` |
| `--freeze-atoms TEXT` | 凍結する原子の 1 始まりインデックスをカンマ区切りで明示的に指定（例: `'1,3,5'`）。`--freeze-links` と併用可、任意の入力形式に適用。 | _None_ |
| `--max-cycles INT` | `opt.max_cycles` に渡されるマクロサイクル上限 | `10000` |
| `--opt-mode TEXT` | 最適化モード: `grad`（`dimer`）または `hess`（`rsirfo`）。`dimer`/`rsirfo` も指定可。サブコマンド別の対応表（`opt` は L-BFGS/RFO、`tsopt` は Dimer/RS-I-RFO）は {ref}`ja-opt-mode-semantics` を参照。 | `hess` |
| `--dump/--no-dump` | 軌跡をダンプ | `False` |
| `-o, --out-dir TEXT` | 出力ディレクトリ | `./result_tsopt/` |
| `--thresh TEXT` | 収束プリセットの上書き（`gau_loose`、`gau`、`gau_tight`、`gau_vtight`、`baker`、`never`） | `baker` |
| `--flatten/--no-flatten` | 余分な虚振動数モードのフラット化ループを有効化（`False` は `flatten_max_iter=0` を強制）。TS 最適化が収束した後、ヘシアン行列の余分な負の固有値モードを反復的にフラット化し、虚振動数が 1 つだけ残るか反復上限に達するまで繰り返します。dimer（dimer ループ）および RS-I-RFO（収束後）の両方に適用 | `False` |
| `--hessian-calc-mode CHOICE` | MLIPヘシアンモード（`Analytical` または `FiniteDifference`） | `FiniteDifference` |
| `--convert-files/--no-convert-files` | PDB または Gaussian 入力用の XYZ/TRJ → PDB/GJF コンパニオン出力を切り替え | `True` |
| `--ref-pdb FILE` | 入力がXYZ/GJFの場合に使用する参照 PDB トポロジー | _None_ |
| `--config FILE` | 明示 CLI オプションより前に適用するベース YAML 設定ファイル | _None_ |
| `--show-config/--no-show-config` | 解決後の設定レイヤーを表示して実行を継続 | `False` |
| `--out-json/--no-out-json` | `out_dir` に機械可読な `result.json` を書き出す。スキーマは [JSON 出力スキーマ](json-output.md) を参照。 | `False` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP バックエンド | `uma` |
| `--solvent TEXT` | xTB 暗黙溶媒（例: `water`）。`none` で無効化 | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB 溶媒モデル | `alpb` |
| `--dry-run/--no-dry-run` | 実行せずに入力/設定を検証し、実行計画を表示。 | `False` |

(ja-flatten-precedence-caveat)=
### `--flatten` 優先順位の注意

```{note}
**`--flatten` はデフォルトで無効です（優先順位の注意）。** `defaults.py` では `flatten_max_iter: 50` が定義されており（下記インライン YAML 表記でも `flatten_max_iter: 50`）、それにもかかわらず CLI の初期化器はコマンドラインで `--flatten` が **明示的に** 渡されない限り `flatten_max_iter = 0` を強制します。実効値は以下のとおりです:

- CLI `--flatten` **未指定** → `flatten_max_iter = 0`（余剰モード除去ループ無効）。YAML の 50 は**無視**されます。
- CLI `--flatten` 指定 → YAML / `defaults.py` の値が有効（デフォルト `flatten_max_iter = 50`）。YAML で `hessian_dimer.flatten_max_iter` や `rsirfo.flatten_max_iter` を上書き可能です。

TS 候補に複数の虚振動数がある場合は、`--flatten` を追加して余分なモードの除去ループを有効にしてください。
```

## 出力
```
out_dir/ (デフォルト:./result_tsopt/)
├─ final_geometry.xyz # 常に書き込み
├─ final_geometry.pdb # 入力がPDBの場合（変換有効時）
├─ final_geometry.gjf # 入力がGaussianの場合（変換有効時）
├─ optimization_all_trj.xyz # --dumpがTrueのときの dimer モードダンプ
├─ optimization_all.pdb # PDB 入力の dimer モードPDB コンパニオン（変換有効時、--dump）
├─ optimization_trj.xyz # --dumpがTrueのときの rsirfo モード軌跡
├─ optimization.pdb # rsirfo モードPDB コンパニオン（変換有効時、--dump）
├─ vib/
│ ├─ imag_±XXXX.Xcm-1_trj.xyz
│ └─ imag_±XXXX.Xcm-1.pdb
└─.dimer_mode.dat # dimer モード方向シード
```

## 終了コード

終了コードは CLI 規約の {ref}`ja-exit-codes` を参照。

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。
- 虚振動数モード**検出**の閾値はデフォルトで 5.0 cm⁻¹（`hessian_dimer.neg_freq_thresh_cm` で変更可能）。この閾値未満の振動数は虚振動数としてカウントされません。選択した `root` は最適化中にどの振動モードを追跡するかを制御します。5 cm⁻¹ 検出閾値と 100 cm⁻¹ 品質ゲートの違いは用語集 {ref}`ja-imaginary-mode-thresholds` を参照。
- `--opt-mode` はワークフロー選択用です（デフォルト: `rsirfo`）。YAML のモードマッピングを手動で変更するのではなく、目的のアルゴリズムに合ったモードを選択してください。
- PHVAの並進/回転射影は `freq` と同じ実装を使用し、メモリ消費を抑えつつ、活性空間の正しい固有ベクトルを保持します。
- 設定の優先順位は {ref}`CLI 規約: 設定の優先順位 <ja-configuration-precedence>` を参照してください。

共通セクションについては [YAML リファレンス](yaml-reference.md) を参照してください。必要な値だけ変更してください。

```{note}
**リファレンスの重複について。** 以下に並ぶ `geom`, `calc`, `opt`, `hessian_dimer`, `rsirfo` の YAML キーは [YAML リファレンス](yaml-reference.md) の正規定義をミラーしたものです。両者が食い違った場合は [YAML リファレンス](yaml-reference.md) と `pdb2reaction/defaults.py` を正とし、本ページの付録は `tsopt` 固有のデフォルト値（例: `out_dir: ./result_tsopt/`、上述の `--flatten` との相互作用）と参照の便宜のためにインライン展開しているだけです。
```

### 共通設定（両モード共通）

`geom` と `calc` のキーは正規定義から変更ありません。詳細は YAML リファレンスの [`geom`](yaml-reference.md#geom) と [`calc`](yaml-reference.md#calc) を参照してください。

`opt` ブロックは [`opt`](yaml-reference.md#opt) と同じキーを使用し、以下が `tsopt` 固有のデフォルトです:

```yaml
opt:
 thresh: baker # tsopt のデフォルト（`opt` は `gau`）
 out_dir: ./result_tsopt/ # tsopt のデフォルト（`opt` は `./result_opt/`）
```

```{note}
**エネルギープラトーによるフォールバック収束。** RS-I-RFO は共通の
`energy_plateau` 設定を参照します。直近 `energy_plateau_window` ステップ（デフォルト 50）の
エネルギーレンジ（max − min）が `energy_plateau_thresh`（デフォルト `1×10⁻⁴ au ≈ 0.06 kcal/mol`）
を下回ると収束と判定されます。大規模 TS 系では MLIP の力のノイズフロア（~4×10⁻⁴ au）が
`baker` max_force 閾値（3×10⁻⁴ au）を上回ることがあり、エネルギーがすでにフラット化していても
力ベース判定に到達しないことがあります。無効化するには `energy_plateau: false` を指定してください。
```

### Dimer モード（`--opt-mode grad`）

`--opt-mode grad`（Hessian Guided Dimer + LBFGS）で使用します。

`hessian_dimer` ブロック全体（内部 `dimer:` とそのネストされた `lbfgs:` を含む）は [`hessian_dimer`](yaml-reference.md#hessian_dimer) に記載されています。内部 `lbfgs:` は [`lbfgs`](yaml-reference.md#lbfgs) セクションを継承し、以下が `tsopt` 固有の上書きです:

```yaml
hessian_dimer:
 dimer:
   lbfgs:
     out_dir: ./result_tsopt/ # tsopt の上書き（defaults.py の値は ./result_opt/）
```

### RS-I-RFO モード（`--opt-mode hess`、デフォルト）

`--opt-mode hess`（RS-I-RFO、デフォルト）で使用します。

`rsirfo` ブロック全体は [`rsirfo`](yaml-reference.md#rsirfo) に記載されています（trust-region と Hessian-update の各キーは [`rfo`](yaml-reference.md#rfo) からも継承）。`tsopt` 固有の上書きは以下のとおりです:

```yaml
rsirfo:
 trust_max: 0.10 # 最大信頼半径 (bohr)
 out_dir: ./result_tsopt/ # tsopt の上書き（defaults.py の値は ./result_opt/）
 hessian_recalc: 500 # N マクロステップごとに exact Hessian を再計算
```

```{tip}
最適化中に TS モードが別のルートに切り替わる場合（例: 複数の虚振動数が存在する場合）は `rsirfo.track_mode_by_overlap: true` を設定してください。
```

```{tip}
TS 収束が遅い場合や最適化中に TS モードが失われる場合は、`rsirfo` セクションの `hessian_recalc` を小さくしてみてください（例: 50--200）。正確なヘシアン再計算の頻度を上げることで、追加のヘシアン評価コストと引き換えに堅牢性が向上します。
```

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け

- [path-search](path-search.md) — TS 候補（HEI）を特定するMEP 探索
- [irc](irc.md) — 最適化されたTSからの反応経路追跡
- [freq](freq.md) — 完全な振動解析と熱化学補正（虚振動数チェックは `tsopt` が内部で実行済み）
- [all](all.md) — 抽出 → MEP → tsopt → IRC（→ オプションで freq/DFT）を連鎖する一気通貫ワークフロー
- [YAML リファレンス](yaml-reference.md) — `hessian_dimer`（Hessian Guided Dimer）と `rsirfo` の完全な設定オプション
- [用語集](glossary.md) — TS、Dimer、RS-I-RFO、ヘシアンの定義
