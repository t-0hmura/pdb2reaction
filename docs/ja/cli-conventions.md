# CLI 規約

このページでは、`pdb2reaction` の全コマンドで使用される規約を説明します。これらの規約を理解することで、正しいコマンドを記述し、よくあるエラーを回避できます。

---

## ブール値オプション

ブール値オプションは root CLI で正規化されます。
次の2記法を受け付けます。

```bash
# 推奨
--tsopt --thermo --no-dft

# 互換記法として受理
--tsopt True --thermo yes --dft 0
```

`--flag` 単独で定義されているオプションでも、互換のため `--no-flag` と `--flag False` を受理します。
`extract` と `fix-altloc` を含むすべてのサブコマンドが Click を CLI バックエンドとして使用します。

よく使うブール値オプション：
- `--tsopt`, `--thermo`, `--dft` — 後処理ステージの有効化
- `--freeze-links` — リンク水素の親原子を凍結（デフォルト: `True`）
- `--dump` — 軌跡ファイルの出力
- `--preopt`, `--endopt` — 前処理/後処理最適化の切り替え
- `--climb` — MEP 探索でクライミングイメージを有効化
- `--convert-files` — PDB/GJF コンパニオンファイルの生成

---

## 段階的ヘルプ（`all`）

`pdb2reaction all` は 2 段階ヘルプです:

```bash
pdb2reaction all --help # 主要オプションのみ
pdb2reaction all --help-advanced # 全オプション
```

以下のコマンドも同じ段階的ヘルプに対応しています（`--help` で主要オプション、`--help-advanced` で全オプション）: `scan`, `scan2d`, `scan3d`, `opt`, `path-opt`, `path-search`, `tsopt`, `freq`, `irc`, `dft`, `add-elem-info`, `trj2fig`, `energy-diagram`, `extract`, `fix-altloc`。

---

## 残基セレクタ

残基セレクタ（Residue Selector）は、基質や抽出中心として使用する残基を指定します。

### 残基名による指定
```bash
-c 'SAM,GPP' # SAM または GPP という名前の残基をすべて選択
-c 'LIG' # LIG という名前の残基をすべて選択
```

### 残基IDによる指定
```bash
-c '123,456' # 残基 123 と 456
-c 'A:123,B:456' # チェーン A の残基 123、チェーン B の残基 456
-c '123A' # 挿入コード A を持つ残基 123
-c 'A:123A' # チェーン A、残基 123、挿入コード A
```

### PDB ファイルによる指定
```bash
-c substrate.pdb # 別の PDB から座標を使用して基質を特定
```

```{note}
残基名で選択する場合、同名の残基が複数あれば**すべて**が含まれ、警告がログに出力されます。
```

```{warning}
**`--selected-resn` は残基 ID を受け付け、残基名は受け付けません。** `extract` と `all` の `--selected-resn` フラグは、名前に反して **残基 ID**（コロン区切り整数、オプションでチェーン/挿入コード付き、例 `A:123A`）を受け付け、3 文字残基名は受け付けません。残基名ベースの基質選択には `-c/--center 'GPP,SAM'` を使用してください。正式な説明は [`extract` の CLI オプション表](extract.md) を参照してください。
```

---

(ja-charge-specification)=
## 電荷の指定

PDB 入力では、`--ligand-charge` を使うと**非標準残基（基質・補因子・金属イオンなど）の形式電荷（formal charge）だけ**を指定すれば、標準アミノ酸やイオンの電荷と合算して全系の電荷が自動計算されます。大きな酵素–基質系で総電荷を手動で数える必要がなくなります。

### 残基別マッピング（推奨）
```bash
-l 'SAM:1,GPP:-3' # SAM は +1、GPP は -3
-l 'LIG:-2' # LIG は -2
```

### 総電荷の明示的上書き
```bash
-q 0 # 総電荷を 0 に強制
-q -1 # 総電荷を -1 に強制
```

### 電荷の解決順序
1. `-q/--charge`（明示的な CLI 上書き）— 最優先
2. 活性部位モデル（バインディングポケット）抽出（アミノ酸、イオン、`--ligand-charge` の合計） — `all` などで `-c/--center` を指定し抽出が有効な場合のみ
3. フォールバックとしての `--ligand-charge`（抽出がスキップされた場合）
4. `.gjf` テンプレートのメタデータ
5. デフォルト: なし（未解決の場合は実行を中断します。`-q`、`.gjf` 電荷メタデータ、または PDB 入力の `--ligand-charge` のいずれかで解決してください）

```{note}
ステップ 2（抽出による電荷導出）は `-c/--center` を伴う `all` などのコマンドでのみ作動します。`opt`/`tsopt`/`freq` などの単独サブコマンドや、`-c` を指定しない場合は抽出がスキップされ、解決順は 1 → 3 → 4 → 5 となります。
```

```{note}
`--ligand-charge` による導出は、PDB 入力のみ（`--ref-pdb` を付けた XYZ/GJF 入力を含む）で電荷が**まだ解決されていない**場合に適用されます。未解決の場合は、`.gjf` メタデータへフォールバックする前に `--ligand-charge` 導出を先に試行します。
```

```{tip}
非標準の残基（基質、補因子、特殊なリガンド）には必ず `--ligand-charge` を指定し、正しい電荷伝播を確保してください。
```

---

## スピン多重度

```bash
-m 1 # 一重項 (singlet)（デフォルト）
-m 2 # 二重項 (doublet)
-m 3 # 三重項 (triplet)
```

```{note}
`all` を含む全サブコマンドで `-m/--multiplicity` を統一して使用します。
```

---

## 原子セレクタ

原子セレクタ（Atom Selector）は、スキャンや拘束に使用する特定の原子を指定します。指定方法は以下の通りです。

### 整数インデックス（デフォルトは1始まり）
```bash
--scan-lists '[(1, 5, 2.0)]' # 原子 1 と 5、ターゲット距離 2.0 Å
```

### PDB形式のセレクタ文字列
```bash
--scan-lists '[("TYR,285,CA", "SAM,309,C10", 2.20)]'
```

セレクタのフィールドは以下で区切れます。
- 空白: `'TYR 285 CA'`
- カンマ: `'TYR,285,CA'`
- スラッシュ: `'TYR/285/CA'`
- バッククォート: `` 'TYR`285`CA' ``
- バックスラッシュ: `'TYR\285\CA'`

3つのトークン（残基名、残基番号、原子名）は任意の順序で指定できます。パーサーは順序が標準的でない場合にフォールバックヒューリスティックを使用します。

---

## 入力ファイル要件

### PDB ファイル
- **水素原子**を含む必要があります（`reduce`、`pdb2pqr`、または Open Babel で追加）
- 列 77–78 に**元素記号**が必要（欠けている場合は `pdb2reaction add-elem-info` を使用）
- 複数の PDB は**同じ原子を同じ順序**で持つ必要があります（座標のみ異なる）

### XYZ および GJF ファイル
- 活性部位モデル抽出をスキップする場合に使用可能（`-c/--center` を省略）
- `.gjf` ファイルは埋め込みメタデータから電荷/スピンのデフォルトを提供可能

---

(ja-exit-codes)=
## 終了コード

`pdb2reaction` サブコマンドは概ね共通の終了コード規約に従いますが、実際に送出されるコードはサブコマンドごとに異なります（下表「主な発生元」列を参照）。各サブコマンドが返し得る終了コードの一覧は、当該サブコマンドのページを確認してください。

| コード | 意味 | 主な発生元 |
|--------|------|-----------|
| `0` | 成功 | すべてのサブコマンド |
| `1` | 予期しないエラー（未処理例外） | すべてのサブコマンド |
| `2` | ゼロステップ長（ステップノルム下限未満）**または** 依存インポート失敗 | `opt`, `tsopt`, `path-opt`; `dft`（PySCF/GPU4PySCF 未インストール） |
| `3` | 最適化失敗 **または** SCF 非収束 | `opt`, `tsopt`, `path-opt`; `dft` |
| `4` | 軌跡書き出しエラー | `path-opt` |
| `5` | HEI エクスポートエラー | `path-opt` |
| `130` | キーボード割り込み (SIGINT) | すべてのサブコマンド |

`irc` や `freq` のように `0 / 1 / 130` しか使わないサブコマンドも同じ規約に従います。単に、これらは現時点で最適化固有のエラーを送出していないだけです。

---

(ja-opt-mode-semantics)=
## `--opt-mode`（サブコマンド依存）

```{warning}
同じ `--opt-mode` トークンでも、サブコマンドによって**選択される最適化アルゴリズムが異なり**、デフォルトも**統一されていません**。レシピをコピーする前に必ずサブコマンドごとの表を確認してください。
```

| サブコマンド | `grad` エイリアス | `hess` エイリアス | デフォルト |
|------------|------------------|------------------|-----------|
| `opt` | L-BFGS (`lbfgs`) | RFO (`rfo`) | `grad` (L-BFGS) |
| `tsopt` | Dimer (`dimer`) | RS-I-RFO (`rsirfo`) | `hess` (RS-I-RFO) |
| `path-opt`（端点 preopt） | L-BFGS | RFO | `grad` |
| `path-search`（端点 preopt） | L-BFGS | RFO | `grad` |
| `scan` / `scan2d` / `scan3d`（端点 preopt） | L-BFGS | RFO | `grad` |
| `all`（pre-opt 段階、`--opt-mode`） | L-BFGS | RFO | `grad` |
| `all`（post-opt — TSOPT プリセット、`--opt-mode-post`） | Dimer (`dimer`) | RS-I-RFO (`rsirfo`) | `hess` |
| `all`（post-opt — IRC 後エンドポイント最適化、`--opt-mode-post`） | L-BFGS | RFO | `hess` |

**受け付けるエイリアス**もサブコマンド固有です:

- `opt` は `grad` / `lbfgs` と `hess` / `rfo` を受け付けます。
- `tsopt` は `grad` / `dimer` と `hess` / `rsirfo` を受け付けます。

したがって `tsopt` に対する `--opt-mode grad` は L-BFGS 最小化ではなく **Dimer TS 探索**です。サブコマンド間で曖昧さを避けたい場合は、明示的なアルゴリズム名（`--opt-mode lbfgs`, `--opt-mode rsirfo` など）を指定してください。

---

## CLI ↔ YAML 名称の不一致

一部の CLI フラグは YAML の対応キーと微妙に名前が異なり、`all` でラップされたときにリネームされるものもあります。完全なマッピング表は [YAML リファレンス → 主要な CLI→YAML マッピング](yaml-reference.md#ja-common-cli-to-yaml-mapping) にあります。特に混同されやすい 2 ケースを以下に示します:

(ja-pressure-vs-pressure-atm)=
### `--pressure` (CLI) vs `pressure_atm` (YAML)

- **CLI フラグ:** `--pressure FLOAT`（`freq` サブコマンド; `all` では `--freq-pressure` として公開）。
- **YAML キー:** `thermo.pressure_atm`（単位接尾辞付き）。
- 両方とも値は **atm** 単位で扱われ、内部で Pa に変換されます。

(ja-engine-vs-dft-engine)=
### `--engine`（単体 `dft`）vs `--dft-engine`（`all` 内）

- **単体 `dft`** サブコマンドではバックエンド選択フラグは **`--engine`**（値: `gpu`, `cpu`）です。
- **`pdb2reaction all`** 内では同じオプションが **`--dft-engine`** にリネームされます（`all` ラッパーで他の engine 系フラグと衝突しないようプレフィックスで区別するため）。
- YAML では両方とも同じ `dft` セクション設定に解決されます。[YAML リファレンス → `dft` セクション](yaml-reference.md#ja-dft-section) を参照してください。

等価なコマンド:

```bash
# 単体 dft
pdb2reaction dft -i ts.xyz -q 0 --engine gpu

# all ラッパー内で同じ処理
pdb2reaction all -i r.pdb p.pdb -c SAM --dft --dft-engine gpu
```

---

## YAML 設定

詳細設定は多層 YAML で渡せます：

```bash
pdb2reaction -i r.pdb p.pdb -q -1 --config my_settings.yaml --out-dir result/
```

利用可能なすべてのオプションは [YAML リファレンス](yaml-reference.md) を参照してください。

(ja-configuration-precedence)=
### 設定の優先順位

設定は以下の順序で解決されます（後のものが前のものを上書き）：

```
組み込みデフォルト  <  --config (YAML)  <  CLI オプション
```

- **組み込みデフォルト** — すべてのパラメータのハードコード値（`pdb2reaction/defaults.py` を参照）。
- **`--config`** — デフォルトを上書きする YAML ファイル。サイト共通やプロジェクト共通の設定に便利です。
- **CLI オプション** — コマンドラインで明示的に指定されたフラグ（例: `-b orb`）。*明示的に指定された*値のみが YAML を上書きし、CLI デフォルトのままのオプションは YAML の値を隠しません。

この優先順位は `all`, `opt`, `tsopt`, `freq`, `irc`, `scan`, `scan2d`, `scan3d`, `path-opt`, `path-search`, `dft` に共通です。あわせて {ref}`YAML リファレンス: 設定の優先順位 <ja-yaml-configuration-precedence>` を参照してください。

---

## 出力ディレクトリ

`-o/--out-dir` で結果の保存先を指定します：

```bash
-o ./my_results/ # カスタム出力ディレクトリ
```

デフォルトの出力ディレクトリ：
- `all`: `./result_all/`
- `extract`: カレントディレクトリまたは指定の `-o`
- `opt`: `./result_opt/`
- `tsopt`: `./result_tsopt/`
- `path-opt`: `./result_path_opt/`
- `path-search`: `./result_path_search/`
- `scan`: `./result_scan/`
- `scan2d`: `./result_scan2d/`
- `scan3d`: `./result_scan3d/`
- `freq`: `./result_freq/`
- `irc`: `./result_irc/`
- `dft`: `./result_dft/`

---

## 関連項目

- [インストール](installation.md) — セットアップと依存関係
- [はじめに](getting-started.md) — 初回実行とワークフロー概要
- [典型エラー別レシピ](recipes-common-errors.md) — 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) — よくあるエラーと対処法
- [YAML リファレンス](yaml-reference.md) — 全設定オプション
