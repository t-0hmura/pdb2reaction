# `all`

## 概要

`pdb2reaction all` は、抽出から解析までの一連の処理を **まとめて実行する最上位コマンド** です。典型的なフローは次のとおりです。

ポケット抽出 →（任意）段階的 UMA スキャン → 再帰的 MEP 探索（`path-search`, GSM/DMF）→ 全系へのマージ →（任意）TS 最適化 + IRC（`tsopt`）→（任意）振動解析・熱化学（`freq`）→（任意）DFT 一点計算（`dft`）

```{important}
`--tsopt` の出力は **TS 候補** です。`all` は自動的に IRC と freq による検証まで実行しますが、結果の虚数モードと端点極小は必ず目視で確認してください。
```

主なモードは 3 つあります。

- **end-to-end（複数構造）** — 反応順に並べた 2 構造以上（PDB/GJF/XYZ）と基質定義を与えます。`all` がポケット抽出→GSM/DMF による MEP 探索→全系テンプレートへのマージまで行い、必要に応じてセグメントごとに TSOPT / freq / DFT を実行します。
- **単一構造 + 段階的スキャン** — 1 つの構造に対して `--scan-lists` を 1 つ以上与えます。スキャンで得られた中間体列を MEP の端点として用います。
 - `--scan-lists` を 1 つだけ渡すと 1 ステージです。
 - 複数ステージは、`--scan-lists` を 1 回指定した後に複数値として渡します（フラグの繰り返し指定はできません）。
- **TSOPT のみ（ポケット TS 最適化）** — 1 つの入力構造に対し、`--scan-lists` を省略して `--tsopt` を指定します。`-c/--center` がある場合はポケットを抽出し、その系で TS 最適化 + IRC（必要に応じて freq / DFT）だけを実行します。

## 最小例

```bash
pdb2reaction all -i R.pdb P.pdb -c "SAM,GPP" --ligand-charge "SAM:1,GPP:-3" --out-dir ./result_all
```

## 出力の見方

- `result_all/summary.log`
- `result_all/summary.yaml`
- `result_all/path_search/mep.pdb`（または `result_all/path_search/seg_*/`）

## よくある例

1. TS 最適化・熱化学・DFT まで一括実行する。

```bash
pdb2reaction all -i R.pdb P.pdb -c "SAM,GPP" --ligand-charge "SAM:1,GPP:-3" \
 --tsopt --thermo --dft --out-dir ./result_all
```

2. 単一構造 + 段階的スキャンを実行する。

```bash
pdb2reaction all -i A.pdb -c "308,309" --scan-lists "[(12,45,1.35)]" "[(10,55,2.20)]" \
 --multiplicity 1 --out-dir ./result_scan_all
```

テンプレートがある場合の XYZ/TRJ → PDB/GJF 変換（付随ファイルの生成）は、全ステージ共通の `--convert-files/--no-convert-files`（既定: `True`）で制御します。


## 使用法
```bash
pdb2reaction all -i INPUT1 [INPUT2...] -c SUBSTRATE [options]
```

### 例
```bash
# 明示的なリガンド電荷と後処理を伴う複数構造アンサンブル
pdb2reaction all -i reactant.pdb product.pdb -c 'GPP,MMT' \
 --ligand-charge 'GPP:-3,MMT:-1' --multiplicity 1 --freeze-links \
 --max-nodes 10 --max-cycles 100 --climb --opt-mode grad \
 --out-dir ./result_all --tsopt --thermo --dft

# 単一構造段階的スキャン + GSM/DMF + TSOPT/freq/DFT
pdb2reaction all -i single.pdb -c '308,309' \
 --scan-lists '[("TYR,285,CA","MMT,309,C10",2.20),("TYR,285,CB","MMT,309,C11",1.80)]' \
 --opt-mode hess --tsopt --thermo --dft

# TSOPT のみワークフロー（経路探索なし）
pdb2reaction all -i reactant.pdb -c 'GPP,MMT' \
 --ligand-charge 'GPP:-3,MMT:-1' --tsopt --thermo --dft
```

## ワークフロー
1. **活性部位ポケット抽出**（`-c/--center` が指定された場合）
 - 基質は PDB パス、残基 ID（`123,124` または `A:123,B:456`）、または残基名（`GPP,MMT`）で指定可能
 - 抽出オプション: `--radius`、`--radius-het2het`、`--include-H2O`、`--exclude-backbone`、`--add-linkH`、`--selected-resn`、`--verbose`
 - 入力ごとのポケット PDB は `<out-dir>/pockets/` に保存。複数構造が提供された場合、ポケットは残基選択ごとに統合
 - **最初のポケットの総電荷**がスキャン/MEP/TSOPT に伝播

2. **オプションの段階的スキャン（単一入力のみ）**
 - 各 `--scan-lists` 引数は UMA スキャンステージを記述する `(i,j,target_Å)` タプルの Python ライクなリスト
 - 単一リテラルは 1 ステージスキャンを実行し、複数リテラルは**順次**実行
 - `--scan-out-dir`、`--scan-one-based`、`--scan-max-step-size`、`--scan-bias-k`、`--scan-relax-max-cycles`、`--scan-preopt`、`--scan-endopt` などの上書きフラグが利用可能
 - ステージエンドポイント（`stage_XX/result.pdb`）が、後続 MEP ステップへ渡される順序付き中間体となる

3. **ポケットでの MEP 探索（再帰的 GSM/DMF）**
 - 抽出されたポケット（または抽出がスキップされた場合は元の全構造）を使用してデフォルトで `path-search` を実行（出力は `<out-dir>/path_search/`）
 - `--no-refine-path` で再帰的精密化なしのシングルパス `path-opt` GSM/DMFチェーンに切り替え

4. **ポケットを全系にマージ**
 - 参照 PDB テンプレートが存在する場合、マージされた `mep_w_ref*.pdb` およびセグメントごとの `mep_w_ref_seg_XX.pdb` ファイルが `<out-dir>/path_search/` に出力

5. **オプションのセグメントごとの後処理**
 - `--tsopt`: 各HEIポケットでTS 最適化を実行、EulerPC IRCで追跡し、セグメントエネルギーダイアグラムを出力
 - `--thermo`: (R, TS, P) で `freq` を呼び出し振動/熱化学データとUMA Gibbsダイアグラムを取得
 - `--dft`: (R, TS, P) でDFT 一点計算を起動しDFTダイアグラムを構築。`--thermo` と組み合わせるとDFT//UMA Gibbsダイアグラムも生成
- 共有の上書きには `--opt-mode`、`--opt-mode-post`（TSOPT/IRC後最適化のプリセット上書き）、`--flatten/--no-flatten`、`--hessian-calc-mode`、`--tsopt-max-cycles`、`--tsopt-out-dir`、`--freq-*`、`--dft-*`、`--dft-engine`（GPU優先）などが含まれる
 - VRAMが十分な場合は `--hessian-calc-mode` を `Analytical` に設定することを強く推奨

6. **TSOPT のみモード**（単一入力、`--tsopt`、`--scan-lists` なし）
 - MEP/マージステージをスキップ。ポケット（または抽出がスキップされた場合は完全入力）で `tsopt` を実行し、EulerPC IRCを実行
 - 高エネルギー側の IRC 終端を反応物 (R) として識別し、エネルギーダイアグラム一式とオプションの freq/DFT 出力を生成


### 電荷とスピンの優先順位

**電荷の解決（優先度順）:**

| 優先度 | ソース | 適用条件 |
|--------|--------|----------|
| 1 | `-q/--charge` | CLI で明示指定 |
| 2 | ポケット抽出 | `-c` 指定時（アミノ酸・イオン・`--ligand-charge` を合算） |
| 3 | `--ligand-charge` | 抽出スキップ時のフォールバック |
| 4 | `.gjf` テンプレート | 埋め込み電荷/スピン情報 |
| 5 | デフォルト | なし（未解決ならエラー） |

**スピンの解決:** `--multiplicity`（CLI） → `.gjf` テンプレート → デフォルト (1)

> **ヒント:** 非標準の基質には `--ligand-charge` を必ず指定し、正しい電荷伝播を確保してください。

### 入力要件
- 抽出有効（`-c/--center`）: 残基同定のため入力は **PDB** が必須。
- 抽出なし: **PDB/XYZ/GJF** を使用可能。
- 複数構造実行は 2 つ以上の構造が必要。


## CLI オプション

> **注記:** 表示されているデフォルト値は、オプション未指定時に使用される内部デフォルトです。

### 入出力オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH...` | 反応順序の2つ以上の完全構造（`--scan-lists` または `--tsopt` のみ単一入力可） | 必須 |
| `--out-dir PATH` | トップレベル出力ディレクトリ | `./result_all/` |
| `--convert-files/--no-convert-files` | XYZ/TRJ → PDB/GJFコンパニオンのグローバルトグル | `True` |
| `--dump/--no-dump` | MEP(GSM/DMF)軌跡を出力。`path-search`/`path-opt` には常時転送され、`scan`/`tsopt` には明示指定時のみ転送。`freq` はデフォルトで dump=True なので `--no-dump` で無効化。 | `False` |
| `--config FILE` | 先に適用するベース YAML | _None_ |
| `--show-config/--no-show-config` | 実行前に解決済み設定を表示 | `False` |
| `--dry-run/--no-dry-run` | 実行せず検証と計画表示のみ行う | `False` |

### 電荷・スピンオプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--ligand-charge TEXT` | 残基別電荷マッピング（例: `GPP:-3,SAM:1`、推奨）。PDB の残基電荷から全系の電荷を自動導出します（手動計算不要）。`-q` 省略時に使用（PDB 入力、または `--ref-pdb` 付き XYZ/GJF） | _None_ |
| `-q, --charge INT` | 総電荷を強制上書き（`--ligand-charge` より優先） | _None_ |
| `-m, --multiplicity INT` | 全下流ステップへ転送されるスピン多重度 | `1` |

### ポケット抽出オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-c, --center TEXT` | 基質指定（PDBパス、残基ID、または残基名） | 抽出に必須 |
| `-r, --radius FLOAT` | ポケット包含カットオフ（Å） | `2.6` |
| `--radius-het2het FLOAT` | ヘテロ–ヘテロカットオフ（Å） | `0.0` |
| `--include-H2O, --include-h2o/--no-include-h2o` | 水分子を含める（HOH/WAT/TIP3/SOL） | `True` |
| `--exclude-backbone/--no-exclude-backbone` | 非基質アミノ酸の主鎖原子を除去 | `True` |
| `--add-linkH/--no-add-linkH` | 切断結合にリンク水素を付加 | `True` |
| `--selected-resn TEXT` | 強制包含残基 | `""` |
| `--freeze-links/--no-freeze-links` | ポケットPDBでリンクHの親を凍結 | `True` |
| `--verbose/--no-verbose` | 抽出器のINFOログを有効化 | `True` |

### MEP 探索オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--mep-mode [gsm\|dmf]` | MEP 探索アルゴリズム: GSM（Growing String Method）または DMF（Direct Max Flux） | `gsm` |
| `--max-nodes INT` | MEP内部ノード数 | `10` |
| `--max-cycles INT` | MEP最大最適化サイクル | `300` |
| `--climb/--no-climb` | 最初のセグメントでTSクライミングを有効化 | `True` |
| `--opt-mode [grad\|hess]` | ワークフロープリセット（`grad` → LBFGS/Dimer、`hess` → RFO/RSIRFO）。コマンド個別実行では `opt --opt-mode grad|hess`、`tsopt --opt-mode grad|hess` を推奨 | `hess` |
| `--thresh TEXT` | 収束プリセット（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） | `gau` |
| `--preopt/--no-preopt` | MEP前にポケット端点を事前最適化 | `True` |
| `--refine-path/--no-refine-path` | True の場合は再帰的 `path-search`、False の場合は `path-opt` を連結して再帰的精密化なしで実行 | `True` |

### UMA計算機オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--workers`, `--workers-per-node` | UMA並列度（workers > 1 で解析ヘシアン無効） | `1`, `1` |
| `--hessian-calc-mode [Analytical\|FiniteDifference]` | 共有UMAヘシアンエンジン | `FiniteDifference` |

### 後処理オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--tsopt/--no-tsopt` | セグメントごとのTS最適化 + IRC を実行 | `False` |
| `--thermo/--no-thermo` | R/TS/Pで振動解析を実行 | `False` |
| `--dft/--no-dft` | R/TS/PでDFT一点計算を実行 | `False` |
| `--opt-mode-post [grad\|hess]` | TSOPT/IRC後最適化のプリセット上書き（`grad` → Dimer/LBFGS、`hess` → RSIRFO/RFO） | `hess` |
| `--thresh-post TEXT` | IRC後エンドポイント最適化の収束プリセット（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） | `baker` |
| `--flatten/--no-flatten` | 余分な虚数モードのフラット化 | `False` |
TSOPT の最適化モードは、`--opt-mode-post`（指定時）→ `--opt-mode`（明示指定時のみ）→ TSOPT の既定（`hess` → `rsirfo`）の順で決まります。

### TSOPT 上書き

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--tsopt-max-cycles INT` | `tsopt --max-cycles` 上書き | `10000` |
| `--tsopt-out-dir PATH` | tsopt出力サブディレクトリ | _None_ |

### Freq 上書き

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--freq-out-dir PATH` | freq出力ディレクトリ上書き | _None_ |
| `--freq-max-write INT` | 最大モード出力数 | `10` |
| `--freq-amplitude-ang FLOAT` | モードアニメーション振幅（Å） | `0.8` |
| `--freq-n-frames INT` | モードアニメーションフレーム数 | `20` |
| `--freq-sort [value\|abs]` | モードソート方法 | `value` |
| `--freq-temperature FLOAT` | 熱化学温度（K） | `298.15` |
| `--freq-pressure FLOAT` | 熱化学圧力（atm） | `1.0` |

### DFT 上書き

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--dft-engine [gpu\|cpu\|auto]` | バックエンド（`auto` はGPU優先） | `gpu` |
| `--dft-out-dir PATH` | DFT出力ディレクトリ上書き | _None_ |
| `--dft-func-basis TEXT` | 汎関数/基底関数ペア | `wb97m-v/def2-tzvpd` |
| `--dft-max-cycle INT` | 最大SCFサイクル | `100` |
| `--dft-conv-tol FLOAT` | SCF収束閾値 | `1e-9` |
| `--dft-grid-level INT` | PySCFグリッドレベル | `3` |

### スキャンオプション（単一入力）

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--scan-lists TEXT...` | 段階的スキャン: `(i,j,target_Å)` タプル | _None_ |
| `--scan-out-dir PATH` | scan出力ディレクトリ上書き | _None_ |
| `--scan-one-based/--no-scan-one-based` | 1始まり/0始まりインデックス | `True` |
| `--scan-max-step-size FLOAT` | 最大ステップサイズ（Å） | `0.20` |
| `--scan-bias-k FLOAT` | 調和バイアス強度（eV/Å²） | `300` |
| `--scan-relax-max-cycles INT` | 緩和サイクル上限 | `10000` |
| `--scan-preopt/--no-scan-preopt` | scan事前最適化 | `True` |
| `--scan-endopt/--no-scan-endopt` | scanステージ終端最適化 | `True` |

## 出力
```text
out_dir/ (デフォルト:./result_all/)
├─ summary.log # テキスト形式の結果要約
├─ summary.yaml # YAML 形式の結果要約
├─ pockets/ # 抽出実行時の入力ごとのポケット PDB
├─ scan/ # 段階的ポケットスキャン結果（--scan-lists提供時）
├─ path_search/ # MEP結果: 軌跡、マージPDB、ダイアグラム
├─ path_search/post_seg_XX/ # 後処理出力（TS最適化、IRC、freq、DFT）
└─ tsopt_single/ # TSOPT のみ出力とIRCエンドポイント
```


- コンソールにはポケット電荷の解決結果、YAML内容、スキャン段数、MEP進行状況（GSM/DMF）、各ステージ時間のサマリーが出力されます。

### `summary.log` の読み方
ログは番号付きセクションで構成されます:
- **[1] グローバル MEP 概要** – イメージ/セグメント数、MEP 軌跡プロットのパス、MEP 全体のエネルギーダイアグラム。
- **[2] セグメント別MEPサマリー（UMAパス）** – セグメントごとの障壁（`ΔE‡`）、反応エネルギー（`ΔE`）、結合変化サマリー。
- **[3] セグメント別後処理（TSOPT / Thermo / DFT）** – TS虚数振動数チェック、IRC出力、UMA/熱化学/DFTのエネルギーテーブル。
- **[4] エネルギーダイアグラム（概要）** – MEP/UMA/Gibbs/DFT 系の図表と、任意の横断サマリー表。
- **[5] 出力ディレクトリ構造** – 生成ファイルを注釈付きでまとめたツリー。

### `summary.yaml` の読み方
YAML はプログラムから処理しやすい形式の要約です。代表的なトップレベルキーは以下のとおりです。
- `out_dir`, `n_images`, `n_segments` – 実行メタデータと総数。
- `segments` – `index`, `tag`, `kind`, `barrier_kcal`, `delta_kcal`, `bond_changes` を含むセグメント配列。
- `energy_diagrams`（任意） – `labels`, `energies_kcal`, `energies_au`, `ylabel`, `image` などを含む図表データ。

`summary.yaml` には `summary.log` にある整形テーブルやファイルツリーは含まれません。

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes_common_errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- 形式電荷が推定できない場合は `--ligand-charge`（数値または残基別マッピング）を必ず指定し、scan/MEP/TSOPT/DFTへ正しい総電荷を伝播させてください。
- マージ用の参照 PDB テンプレートは元の入力から自動導出されます。`path-search` の `--ref-full-pdb` はこのラッパーでは意図的に隠されています。
- 収束プリセット: `--thresh` の既定は `gau`、`--thresh-post` の既定は `baker`。
- 抽出半径: `--radius` または `--radius-het2het` に `0` を渡すと、内部で `0.001 Å` にクランプされます。
- ダイアグラムのエネルギーは反応物（最初の状態）基準の kcal/mol で報告されます。
- `-c/--center` を省略すると抽出をスキップして全構造を MEP/tsopt/freq/DFT に渡しますが、単一構造実行には `--scan-lists` か `--tsopt` が引き続き必要です。


`all` は YAML の多層指定をサポートします:

- `--config FILE`: ベース設定。

適用順序:

`defaults < config < CLI < override-yaml`

解決後の YAML が**すべての**呼び出されるサブコマンドに転送されます。各ツールは独自ドキュメントに記載のセクションを読み取ります:

| サブコマンド | YAML セクション |
|------------|-----------------|
| [`path-search`](path_search.md) | `geom`, `calc`, `gs`, `opt`, `sopt`, `bond`, `search` |
| [`scan`](scan.md) | `geom`, `calc`, `opt`, `lbfgs`, `rfo`, `bias`, `bond` |
| [`tsopt`](tsopt.md) | `geom`, `calc`, `opt`, `hessian_dimer`, `rsirfo` |
| [`freq`](freq.md) | `geom`, `calc`, `freq`, `thermo` |
| [`dft`](dft.md) | `dft` |


**最小例:**
```yaml
calc:
 model: uma-s-1p1
 hessian_calc_mode: Analytical # VRAM に余裕がある場合推奨
gs:
 max_nodes: 12
 climb: true
dft:
 grid_level: 6
```

すべての YAML オプションの完全なリファレンスについては、**[YAML 設定リファレンス](yaml_reference.md)** を参照してください。

---

## 関連項目

- [インストール](installation.md) — セットアップと依存関係
- [はじめに](getting_started.md) — 初回実行とワークフロー概要
- [概念とワークフロー](concepts.md) — ポケット、セグメント、ステージの全体像
- [extract](extract.md) — 単独のポケット抽出（`all` が内部で呼び出し）
- [path-search](path_search.md) — 単独のMEP 探索（`all` が内部で呼び出し）
- [tsopt](tsopt.md) — 単独のTS最適化
- [freq](freq.md) — 単独の振動解析
- [dft](dft.md) — 単独のDFT計算
- [典型エラー別レシピ](recipes_common_errors.md) — 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) — よくあるエラーと対処法
- [YAML リファレンス](yaml_reference.md) — 全YAML設定オプション
- [用語集](glossary.md) — MEP、TS、IRC、GSM、DMFの定義
