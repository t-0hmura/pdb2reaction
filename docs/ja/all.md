# `all`

## 概要

`pdb2reaction all` は、抽出から解析までの一連の処理を **まとめて実行する最上位コマンド** です。典型的なフローは次のとおりです。

ポケット抽出 →（任意）段階的 UMA スキャン → 再帰的 MEP 探索（`path_search`, GSM/DMF）→ 全系へのマージ →（任意）TS 最適化 + IRC（`tsopt`）→（任意）振動解析・熱化学（`freq`）→（任意）DFT 一点計算（`dft`）

主なモードは 3 つあります。

- **エンドツーエンド（複数構造）** — 反応順に並べた 2 構造以上（PDB/GJF/XYZ）と基質定義を与えます。`all` がポケット抽出→GSM/DMF による MEP 探索→全系テンプレートへのマージまで行い、必要に応じてセグメントごとに TSOPT / freq / DFT を実行します。
- **単一構造 + 段階的スキャン** — 1 つの構造に対して `--scan-lists` を 1 つ以上与えます。スキャンで得られた中間体列を MEP の端点として用います。
  - `--scan-lists` を 1 つだけ渡すと 1 ステージです。
  - 複数ステージは、`--scan-lists` を 1 回指定した後に複数値として渡します（フラグの繰り返し指定はできません）。
- **TSOPT のみ（ポケット TS 最適化）** — 1 つの入力構造に対し、`--scan-lists` を省略して `--tsopt True` を指定します。`-c/--center` がある場合はポケットを抽出し、その系で TS 最適化 + IRC（必要に応じて freq / DFT）だけを実行します。

テンプレートがある場合の XYZ/TRJ → PDB/GJF 変換（付随ファイルの生成）は、全ステージ共通の `--convert-files {True\|False}`（既定: `True`）で制御します。


## 使用法
```bash
pdb2reaction all -i INPUT1 [INPUT2 ...] -c SUBSTRATE [options]
```

### 例
```bash
# 明示的なリガンド電荷と後処理を伴う複数構造アンサンブル
pdb2reaction all -i reactant.pdb product.pdb -c 'GPP,MMT' \
    --ligand-charge 'GPP:-3,MMT:-1' --mult 1 --freeze-links True \
    --max-nodes 10 --max-cycles 100 --climb True --opt-mode light \
    --out-dir result_all_${date} --tsopt True --thermo True --dft True

# 単一構造段階的スキャン + GSM/DMF + TSOPT/freq/DFT
pdb2reaction all -i single.pdb -c '308,309' \
    --scan-lists '[("TYR,285,CA","MMT,309,C10",2.20),("TYR,285,CB","MMT,309,C11",1.80)]' \
    --opt-mode heavy --tsopt True --thermo True --dft True

# TSOPT のみワークフロー（経路探索なし）
pdb2reaction all -i reactant.pdb -c 'GPP,MMT' \
    --ligand-charge 'GPP:-3,MMT:-1' --tsopt True --thermo True --dft True
```

## ワークフロー
1. **活性部位ポケット抽出**（`-c/--center` が提供された場合）
   - 基質はPDB パス、残基ID（`123,124` または `A:123,B:456`）、または残基名（`GPP,MMT`）で指定可能
   - 抽出オプション: `--radius`、`--radius-het2het`、`--include-H2O`、`--exclude-backbone`、`--add-linkH`、`--selected-resn`、`--verbose`
   - 入力ごとのポケット PDBは `<out-dir>/pockets/` に保存。複数構造が提供された場合、ポケットは残基選択ごとに統合
   - **最初のポケットの総電荷**がスキャン/MEP/TSOPTに伝播

2. **オプションの段階的スキャン（単一入力のみ）**
   - 各 `--scan-lists` 引数はUMA スキャンステージを記述する `(i,j,target_Å)` タプルのPythonライクなリスト
   - 単一リテラルは1ステージスキャンを実行; 複数リテラルは**順次**実行
   - スキャンは電荷/スピン、`--freeze-links`、UMA最適化プリセット（`--opt-mode`）、`--args-yaml`、`--preopt` を継承。`--dump` はこのコマンドで明示指定された場合のみスキャンへ転送され、それ以外は scan 側のデフォルト（`False`）を使用
   - `--scan-out-dir`、`--scan-one-based`、`--scan-max-step-size`、`--scan-bias-k`、`--scan-relax-max-cycles`、`--scan-preopt`、`--scan-endopt` などの上書きフラグが利用可能
   - ステージエンドポイント（`stage_XX/result.pdb`）が後続MEPステップに供給される順序付き中間体となる

3. **ポケットでのMEP 探索（再帰的GSM/DMF）**
   - 抽出されたポケット（または抽出がスキップされた場合は元の全構造）を使用してデフォルトで `path_search` を実行
   - `--refine-path False` で再帰的精密化なしのシングルパス `path-opt` GSM/DMFチェーンに切り替え

4. **ポケットを全系にマージ**
   - 参照 PDB テンプレートが存在する場合、マージされた `mep_w_ref*.pdb` およびセグメントごとの `mep_w_ref_seg_XX.pdb` ファイルが `<out-dir>/path_search/` に出力

5. **オプションのセグメントごとの後処理**
   - `--tsopt True`: 各HEIポケットでTS 最適化を実行、EulerPC IRCで追跡し、セグメントエネルギーダイアグラムを出力
   - `--thermo True`: (R, TS, P) で `freq` を呼び出し振動/熱化学データとUMA Gibbsダイアグラムを取得
   - `--dft True`: (R, TS, P) でDFT 一点計算を起動しDFTダイアグラムを構築。`--thermo True` と組み合わせるとDFT//UMA Gibbsダイアグラムも生成
   - 共有の上書きには `--opt-mode`、`--opt-mode-post`（TSOPT/IRC後最適化のプリセット上書き）、`--flatten-imag-mode`、`--hessian-calc-mode`、`--tsopt-max-cycles`、`--tsopt-out-dir`、`--freq-*`、`--dft-*`、`--dft-engine`（GPU優先）などが含まれる
   - VRAMが十分な場合は `--hessian-calc-mode` を `Analytical` に設定することを強く推奨

6. **TSOPT のみモード**（単一入力、`--tsopt True`、`--scan-lists` なし）
   - MEP/マージステージをスキップ。ポケット（または抽出がスキップされた場合は完全入力）で `tsopt` を実行し、EulerPC IRCを実行
   - 高エネルギー側のIRC終端を反応物 (R) として識別し、同じ種類のエネルギーダイアグラムとオプションの freq/DFT 出力を生成


### 電荷とスピンの優先順位

**電荷の解決（優先度順）:**

| 優先度 | ソース | 適用条件 |
|--------|--------|----------|
| 1 | `-q/--charge` | CLI で明示指定 |
| 2 | ポケット抽出 | `-c` 指定時（アミノ酸・イオン・`--ligand-charge` を合算） |
| 3 | `--ligand-charge`（数値） | 抽出失敗時またはスキップ時のフォールバック |
| 4 | `.gjf` テンプレート | 埋め込み電荷/スピン情報 |
| 5 | デフォルト | 0 |

**スピンの解決:** `--mult`（CLI） → `.gjf` テンプレート → デフォルト (1)

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
| `-i, --input PATH...` | 反応順序の2つ以上の完全構造（`--scan-lists` または `--tsopt True` のみ単一入力可） | 必須 |
| `--out-dir PATH` | トップレベル出力ディレクトリ | `./result_all/` |
| `--convert-files {True\|False}` | XYZ/TRJ → PDB/GJFコンパニオンのグローバルトグル | `True` |
| `--dump BOOLEAN` | MEP(GSM/DMF)軌跡を出力 | `False` |
| `--args-yaml FILE` | 全サブコマンドへそのまま転送されるYAML | _None_ |

### 電荷・スピンオプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--ligand-charge TEXT` | 未知残基の残基別マッピングまたは総電荷（推奨） | _None_ |
| `-q, --charge INT` | 総電荷を強制上書き（`--ligand-charge` より優先） | _None_ |
| `-m, --mult INT` | 全下流ステップへ転送されるスピン多重度 | `1` |

### ポケット抽出オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-c, --center TEXT` | 基質指定（PDBパス、残基ID、または残基名） | 抽出に必須 |
| `-r, --radius FLOAT` | ポケット包含カットオフ（Å） | `2.6` |
| `--radius-het2het FLOAT` | ヘテロ–ヘテロカットオフ（Å） | `0.0` |
| `--include-H2O BOOLEAN` | 水分子を含める（HOH/WAT/TIP3/SOL） | `True` |
| `--exclude-backbone BOOLEAN` | 非基質アミノ酸の主鎖原子を除去 | `True` |
| `--add-linkH BOOLEAN` | 切断結合にリンク水素を付加 | `True` |
| `--selected-resn TEXT` | 強制包含残基 | `""` |
| `--freeze-links BOOLEAN` | ポケットPDBでリンクHの親を凍結 | `True` |
| `--verbose BOOLEAN` | 抽出器のINFOログを有効化 | `True` |

### MEP探索オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--max-nodes INT` | MEP内部ノード数 | `10` |
| `--max-cycles INT` | MEP最大最適化サイクル | `300` |
| `--climb BOOLEAN` | 最初のセグメントでTSクライミングを有効化 | `True` |
| `--opt-mode [light\|heavy]` | 最適化プリセット（light → LBFGS/Dimer、heavy → RFO/RSIRFO） | `light` |
| `--thresh TEXT` | 収束プリセット（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`） | `gau` |
| `--preopt BOOLEAN` | MEP前にポケット端点を事前最適化 | `True` |

### UMA計算機オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--workers`, `--workers-per-node` | UMA並列度（workers > 1 で解析ヘシアン無効） | `1`, `1` |
| `--hessian-calc-mode [Analytical\|FiniteDifference]` | 共有UMAヘシアンエンジン | `FiniteDifference` |

### 後処理オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--tsopt BOOLEAN` | セグメントごとのTS最適化 + IRC を実行 | `False` |
| `--thermo BOOLEAN` | R/TS/Pで振動解析を実行 | `False` |
| `--dft BOOLEAN` | R/TS/PでDFT一点計算を実行 | `False` |
| `--opt-mode-post [light\|heavy]` | TSOPT/IRC後最適化のプリセット | _None_ |
| `--thresh-post TEXT` | IRC後エンドポイント最適化の収束プリセット | `baker` |
| `--flatten-imag-mode {True\|False}` | 余分な虚数モードのフラット化 | `False` |

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
| `--scan-one-based BOOLEAN` | 1始まり/0始まりインデックス | `True` |
| `--scan-max-step-size FLOAT` | 最大ステップサイズ（Å） | `0.20` |
| `--scan-bias-k FLOAT` | 調和バイアス強度（eV/Å²） | `100` |
| `--scan-relax-max-cycles INT` | 緩和サイクル上限 | `10000` |
| `--scan-preopt BOOLEAN` | scan事前最適化 | `True` |
| `--scan-endopt BOOLEAN` | scanステージ終端最適化 | `True` |

## 出力
```text
out_dir/ (デフォルト: ./result_all/)
├─ summary.log               # クイック検査用フォーマット済みサマリー
├─ summary.yaml              # YAML バージョンサマリー
├─ pockets/                  # 抽出実行時の入力ごとのポケット PDB
├─ scan/                     # 段階的ポケットスキャン結果（--scan-lists提供時）
├─ path_search/              # MEP結果: 軌跡、マージPDB、ダイアグラム
├─ path_search/post_seg_XX/  # 後処理出力（TS最適化、IRC、freq、DFT）
└─ tsopt_single/             # TSOPT のみ出力とIRCエンドポイント
```


- コンソールにはポケット電荷の解決結果、YAML内容、スキャン段数、MEP進行状況（GSM/DMF）、各ステージ時間のサマリーが出力されます。

### `summary.log` の読み方
ログは番号付きセクションで構成されます:
- **[1] グローバルMEP概要** – 画像/セグメント数、MEP軌跡プロットのパス、MEP全体のエネルギーダイアグラム。
- **[2] セグメント別MEPサマリー（UMAパス）** – セグメントごとの障壁（`ΔE‡`）、反応エネルギー（`ΔE`）、結合変化サマリー。
- **[3] セグメント別後処理（TSOPT / Thermo / DFT）** – TS虚数振動数チェック、IRC出力、UMA/熱化学/DFTのエネルギーテーブル。
- **[4] エネルギーダイアグラム（概要）** – MEP/UMA/Gibbs/DFT 系の図表と、任意の横断サマリー表。
- **[5] 出力ディレクトリ構造** – 生成ファイルを注釈付きでまとめたツリー。

### `summary.yaml` の読み方
YAMLは機械可読サマリーです。代表的なトップレベルキーは以下です:
- `out_dir`, `n_images`, `n_segments` – 実行メタデータと総数。
- `segments` – `index`, `tag`, `kind`, `barrier_kcal`, `delta_kcal`, `bond_changes` を含むセグメント配列。
- `energy_diagrams`（任意） – `labels`, `energies_kcal`, `energies_au`, `ylabel`, `image` などを含む図表データ。

`summary.yaml` には `summary.log` にある整形テーブルやファイルツリーは含まれません。

## 注意事項
- 形式電荷が推定できない場合は `--ligand-charge`（数値または残基別マッピング）を必ず指定し、scan/MEP/TSOPT/DFTへ正しい総電荷を伝播させてください。
- マージ用の参照 PDB テンプレートは元の入力から自動導出されます。`path_search` の `--ref-full-pdb` はこのラッパーでは意図的に隠されています。
- 収束プリセット: `--thresh` の既定は `gau`、`--thresh-post` の既定は `baker`。
- 抽出半径: `--radius` または `--radius-het2het` に `0` を渡すと、内部で `0.001 Å` にクランプされます。
- ダイアグラムのエネルギーは反応物（最初の状態）基準の kcal/mol で報告されます。
- `-c/--center` を省略すると抽出をスキップして全構造を MEP/tsopt/freq/DFT に渡しますが、単一構造実行には `--scan-lists` か `--tsopt True` が引き続き必要です。
- `--args-yaml` で全計算器を単一設定から制御できます。YAMLはCLIを上書きします。

## YAML 設定（`--args-yaml`）

同じ YAML ファイルが**すべての**呼び出されるサブコマンドにそのまま転送されます。各ツールは独自のドキュメントに記載されているセクションを読み取ります:

| サブコマンド | YAML セクション |
|------------|-----------------|
| [`path_search`](path_search.md) | `geom`, `calc`, `gs`, `opt`, `sopt`, `bond`, `search` |
| [`scan`](scan.md) | `geom`, `calc`, `opt`, `lbfgs`, `rfo`, `bias`, `bond` |
| [`tsopt`](tsopt.md) | `geom`, `calc`, `opt`, `hessian_dimer`, `rsirfo` |
| [`freq`](freq.md) | `geom`, `calc`, `freq`, `thermo` |
| [`dft`](dft.md) | `dft` |

> **注記:** CLI と YAML の両方が指定された場合、YAML が優先されます。

**最小例:**
```yaml
calc:
  model: uma-s-1p1
  hessian_calc_mode: Analytical  # VRAM に余裕がある場合推奨
gs:
  max_nodes: 12
  climb: true
dft:
  grid_level: 6
```

すべての YAML オプションの完全なリファレンスについては、**[YAML 設定リファレンス](yaml-reference.md)** を参照してください。

---

## 関連項目

- [はじめに](getting-started.md) — インストールと初回実行
- [概念とワークフロー](concepts.md) — ポケット、セグメント、ステージの全体像
- [extract](extract.md) — 単独のポケット抽出（`all` が内部で呼び出し）
- [path-search](path_search.md) — 単独のMEP探索（`all` が内部で呼び出し）
- [tsopt](tsopt.md) — 単独のTS最適化
- [freq](freq.md) — 単独の振動解析
- [dft](dft.md) — 単独のDFT計算
- [トラブルシューティング](troubleshooting.md) — よくあるエラーと対処法
- [YAML リファレンス](yaml-reference.md) — 全YAML設定オプション
- [用語集](glossary.md) — MEP、TS、IRC、GSM、DMFの定義