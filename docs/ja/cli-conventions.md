# CLI 規約

## ブール値オプション

ブール値オプションは root CLI で正規化されます。
次の 2 記法を受け付けます。

```bash
# 推奨
--tsopt --thermo --no-dft

# 互換記法として受理
--tsopt True --thermo yes --dft 0
```

サブコマンドのソース側登録形式（value-style `type=click.BOOL` / flag-pair）に関係なく両形式を受理します。`bool_compat` が value-style flag に `--no-<flag>` synthetic alias を、flag-pair に value-style alias を生成するため、toggle と value 表記はサブコマンド横断で互換です。

よく使うブール値オプション：
- `--tsopt`, `--thermo`, `--dft` — 後処理ステージの有効化
- `--freeze-links` — リンク水素の親原子を凍結（デフォルト: `True`）
- `--dump` — 軌跡ファイルの出力
- `--preopt`, `--endopt` — 前処理/後処理最適化の切り替え
- `--climb` — MEP 探索でクライミングイメージを有効化
- `--convert-files` — 対応する PDB/GJF ファイルの生成

## 段階的ヘルプ（`all`）

`pdb2reaction all` は 2 段階ヘルプです:

```bash
pdb2reaction all --help # 主要オプションのみ
pdb2reaction all --help-advanced # 全オプション
```

以下のコマンドも同じ段階的ヘルプに対応しています（`--help` で主要オプション、`--help-advanced` で全オプション）: `scan`, `scan2d`, `scan3d`, `opt`, `path-opt`, `path-search`, `tsopt`, `freq`, `irc`, `dft`, `sp`, `add-elem-info`, `trj2fig`, `energy-diagram`, `bond-summary`, `extract`, `fix-altloc`。

(ja-verbosity-levels)=

## ログ詳細度 (verbosity)

`-v/--verbose LEVEL` は 0〜3 の整数 (**デフォルト 2**) で、各コマンドのコンソール出力量を決めます。コマンドごとのオプションなので、サブコマンドと一緒に指定します (例: `pdb2reaction opt -v 1 ...`)。4 段階は全コマンド共通で、各コマンドページはそのコマンド固有の出力だけを説明します。

| レベル | 表示内容 |
|---|---|
| `-v 0` | 無出力。成功は終了コードと出力成果物で確認します。 |
| `-v 1` | マイルストーンのみ: バージョン、入力要約、主要設定、出力先、dry-run / 最終ステータス。banner・`[command]`・`[mode]`・config dump は出ません。 |
| `-v 2` | デフォルト。banner、`[command]`、`[mode]`、ステージ進捗、主要な optimizer サイクル表、終了ステータス、Hessian 1 行要約、thermo / DFT 要約、経過時間を追加します。 |
| `-v 3` | デバッグ: resolved config / dry-run plan、backend DEBUG、raw optimizer・内部座標の詳細、`[HessianTiming]`、`[HessianVRAM]`。 |

意味的な失敗はどのレベルでも失敗です。`-v 3` でのみ現れる `Traceback` も実行失敗を意味します。

## 残基セレクタ

残基セレクタ（Residue Selector）は、基質や抽出中心として使用する残基を指定します。

### 残基名による指定
```bash
-c 'SAM,GPP' # SAM または GPP という名前の残基をすべて選択
-c 'LIG' # LIG という名前の残基をすべて選択
```

### 残基 ID による指定
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

(ja-selected-resn-takes-ids)=
### `--selected-resn` は残基 ID を取る（残基名ではない）

```{warning}
**`--selected-resn` は残基 ID を受け付け、残基名は受け付けません。** `extract` と `all` の `--selected-resn` フラグは、名前に反して **残基 ID**（コロン区切り整数、オプションでチェーン/挿入コード付き、例 `A:123A`）を受け付け、3 文字残基名は受け付けません。残基名トークン（例 `'TYR,GLU'`）を渡すと `ValueError("Invalid residue specifier 'TYR'. Use '123', '123A', 'A:123', or 'A:123A'.")` が発生します。残基名ベースの基質選択には `-c/--center 'GPP,SAM'` を使用してください。正式な説明は [`extract` の CLI オプション表](extract.md) を参照してください。
```

---

(ja-charge-specification)=
## 電荷の指定

PDB 入力では、`--ligand-charge/-l` を使うと**非標準残基（基質・補因子・金属イオンなど）の電荷だけ**を指定すれば、標準アミノ酸やイオンの電荷と合算して全系の電荷が自動計算されます。大きな酵素–基質系で総電荷を手動で数える必要がなくなります。

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
2. 活性部位モデル（バインディングポケット）抽出（アミノ酸、イオン、`--ligand-charge/-l` の合計） — `all` などで `-c/--center` を指定し抽出が有効な場合のみ
3. フォールバックとしての `--ligand-charge/-l`（抽出がスキップされた場合）
4. `.gjf` テンプレートのメタデータ
5. デフォルト: なし（未解決の場合は実行を中断します。`-q`、`.gjf` 電荷メタデータ、または PDB 入力の `--ligand-charge/-l` のいずれかで解決してください）

```{note}
ステップ 2（抽出による電荷導出）は `-c/--center` を伴う `all` などのコマンドでのみ作動します。`opt`/`tsopt`/`freq` などの単独サブコマンドや、`-c` を指定しない場合は抽出がスキップされ、解決順は 1 → 3 → 4 → 5 となります。
```

```{note}
`--ligand-charge/-l` による導出は、PDB 入力のみ（`--ref-pdb` を付けた XYZ/GJF 入力を含む）で電荷が**まだ解決されていない**場合に適用されます。未解決の場合は、`.gjf` メタデータへフォールバックする前に `--ligand-charge/-l` 導出を先に試行します。
```

```{tip}
非標準の残基（基質、補因子、特殊なリガンド）には必ず `--ligand-charge/-l` を指定し、正しい電荷伝播を確保してください。
```

## スピン多重度

```bash
-m 1 # 一重項 (singlet)（デフォルト）
-m 2 # 二重項 (doublet)
-m 3 # 三重項 (triplet)
```

```{note}
`all` を含む全サブコマンドで `-m/--multiplicity` を統一して使用します。
```

## 原子セレクタ

原子セレクタ（Atom Selector）は、スキャンや拘束に使用する特定の原子を指定します。指定方法は以下の通りです。

### 整数インデックス（デフォルトは 1 始まり）
```bash
--scan-lists '[(1, 5, 2.0)]' # 原子 1 と 5、ターゲット距離 2.0 Å
```

### PDB 形式のセレクタ文字列
```bash
--scan-lists '[("TYR,285,CA", "SAM,309,C10", 2.20)]'
```

セレクタのフィールドは以下で区切れます。
- 空白: `'TYR 285 CA'`
- カンマ: `'TYR,285,CA'`
- スラッシュ: `'TYR/285/CA'`
- バッククォート: `` 'TYR`285`CA' ``
- バックスラッシュ: `'TYR\285\CA'`

3 つのトークン（残基名、残基番号、原子名）は任意の順序で指定できます。パーサーは順序が標準的でない場合にフォールバックヒューリスティックを使用します。

---

(ja-scan-list-spec)=
### スキャンリスト仕様

`-s/--scan-lists`（`scan` / `scan2d` / `scan3d` / `all` で使用）は YAML/JSON スペックファイルパス、または 1 個以上のインライン Python リテラルを受け付けます。複雑な設定や複数ステージの実行にはファイルが、単純な 1 ステージのみの場合にはインラインリテラルが適しています。

#### YAML/JSON スペックファイルの書式（推奨）

ルートはマッピング形式で、リスト・オブ・タプルのキーは `scan` では `stages`、`scan2d`/`scan3d` では `pairs` です。

```yaml
one_based: true # 任意。未指定時は CLI の --one-based を使用
stages: # scan 用
  - [[1, 5, 1.35]]
  - [[1, 5, 2.20], [2, 8, 1.80]]
```

```yaml
one_based: true # 任意
pairs: # scan2d（要素は 2 つちょうど） / scan3d（要素は 3 つちょうど）
  - [1, 5, 1.30, 3.10]
  - [2, 8, 1.20, 3.20]
```

- `stages` / `pairs` は必須です。
- `scan` の各ステージは `(i, j, target_Å)` 3 要素タプルのリストです。
- `scan2d`/`scan3d` の各軸は `(i, j, low_Å, high_Å)` の 4 要素タプルです。
- インデックスは整数または PDB セレクタのどちらでも指定できます（インラインリテラルと同じ）。

#### インライン Python リテラルの書式

各リテラルは Python リストです。シェルのクォート処理に注意が必要です。

```
-s '[(原子1, 原子2, ターゲット距離Å), ...]'      # scan: 3 要素タプル
-s '[(原子1, 原子2, 下限Å, 上限Å), ...]'         # scan2d / scan3d: 4 要素タプル
```

- リスト全体を **シングルクォート** `'...'` で囲みます（シェルがカッコや空白を解釈しないようにするため）。
- `scan` では **1 リテラル = 1 ステージ**です。複数ステージを実行するには、**1 つの `-s/--scan-lists` フラグの後に複数リテラル**を並べます。
- `scan2d`/`scan3d` ではリテラルは **1 つだけ** を受け付けます（複数ステージは非対応）。`scan2d` ではちょうど 2 つ、`scan3d` ではちょうど 3 つの 4 要素タプルを含む必要があります。

##### 原子の指定方法

原子は**整数インデックス**または **PDB セレクタ文字列**で指定します。

| 方法 | 例 | 備考 |
| --- | --- | --- |
| 整数インデックス | `(1, 5, 2.0)` | デフォルトは 1 始まり（`--one-based`） |
| PDB セレクタ | `("TYR,285,CA", "SAM,309,C10", 2.0)` | 残基名、残基番号、原子名の 3 要素タプル |

PDB セレクタのトークンは、カンマ `,`、スペース、スラッシュ `/`、バッククォート `` ` ``、バックスラッシュ `\` のいずれでも区切れます。トークンの順序も任意です。

```bash
# 以下はすべて同じ原子を指定:
"TYR,285,CA"
"TYR 285 CA"
"TYR/285/CA"
"285,TYR,CA" # 順序は自由
```

##### クォートの規則

```bash
# 正しい: シングルクォートでリストを囲み、セレクタはダブルクォート
-s '[("TYR,285,CA","SAM,309,C10",1.35)]'

# 正しい: 整数インデックスなら内側のダブルクォートは不要
-s '[(1, 5, 2.0)]'

# 非推奨: ダブルクォートで外側を囲むとエスケープが必要
-s "[(\"TYR,285,CA\",\"SAM,309,C10\",1.35)]"
```

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
| `2` | ゼロステップ長（ステップノルム下限未満）**または** 依存インポート失敗 | `opt`, `tsopt`; `dft`（PySCF 未インストール。`--engine gpu` での GPU4PySCF 欠落は ClickException = code 1） |
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
| `path-search`（HEI±1 / kink ノードの単一構造 optimizer） | L-BFGS | RFO | `grad` |
| `scan` / `scan2d` / `scan3d`（grid relaxation） | L-BFGS | RFO | `grad` |
| `all`（pre-opt 段階、`--opt-mode`） | L-BFGS | RFO | `grad` |
| `all`（post-opt — TSOPT プリセット、`--opt-mode-post`） | Dimer (`dimer`) | RS-I-RFO (`rsirfo`) | `hess` |
| `all`（post-opt — IRC 後エンドポイント最適化、`--opt-mode-post`） | L-BFGS | RFO | `hess` |

**受け付けるエイリアス**もサブコマンド固有です:

- `opt` は `grad` / `lbfgs` と `hess` / `rfo` を受け付けます。
- `tsopt` は `grad` / `dimer` と `hess` / `rsirfo` に加え、`trim`（TRIM/Helgaker）と `rsprfo`（RS-P-RFO/Banerjee）も単独の `--opt-mode` 値として受け付けます。
- `scan` / `scan2d` / `scan3d` / `path-opt` / `path-search` / `all` は `grad` / `hess` のみ受け付けます（アルゴリズム名 alias なし）。`all` の `--opt-mode-post` も `grad` / `hess` のみです。

したがって `tsopt` に対する `--opt-mode grad` は L-BFGS 最小化ではなく **Dimer TS 探索**です。曖昧さを避けたい場合は、各サブコマンドが受け付けるアルゴリズム名を使用してください: `opt` では `--opt-mode lbfgs|rfo`、`tsopt` では `--opt-mode dimer|rsirfo`。（他のサブコマンドは `grad` / `hess` のみ受け付けます。）

## CLI ↔ YAML 名称の不一致

一部の CLI フラグは YAML の対応キーと微妙に名前が異なり、`all` でラップされたときにリネームされるものもあります。完全なマッピング表は {ref}`YAML リファレンスの主要な CLI→YAML マッピング <ja-common-cli-to-yaml-mapping>` にあります。特に混同されやすい 2 ケースを以下に示します:

(ja-pressure-vs-pressure-atm)=
### `--pressure` (CLI) vs `pressure_atm` (YAML)

- **CLI フラグ:** `--pressure FLOAT`（`freq` サブコマンド; `all` では `--freq-pressure` として公開）。
- **YAML キー:** `thermo.pressure_atm`（単位接尾辞付き）。
- 両方とも値は **atm** 単位で扱われ、内部で Pa に変換されます。

(ja-engine-vs-dft-engine)=
### `--engine`（単体 `dft`）vs `--dft-engine`（`all` 内）

- **単体 `dft`** サブコマンドではバックエンド選択フラグは **`--engine`**（値: `gpu`, `cpu`）です。
- **`pdb2reaction all`** 内では同じオプションが **`--dft-engine`** にリネームされます（`all` ラッパーで他の engine 系フラグと衝突しないようプレフィックスで区別するため）。
- YAML では両方とも同じ `dft` セクション設定に解決されます。{ref}`YAML リファレンスの dft セクション <ja-dft-section>` を参照してください。

等価なコマンド:

```bash
# 単体 dft
pdb2reaction dft -i ts.xyz -q 0 --engine gpu

# all ラッパー内で同じ処理
pdb2reaction all -i r.pdb p.pdb -c SAM --dft --dft-engine gpu
```

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

- **組み込みデフォルト** — すべてのパラメータのハードコード値（`pdb2reaction/core/defaults.py` を参照）。
- **`--config`** — デフォルトを上書きする YAML ファイル。サイト共通やプロジェクト共通の設定に便利です。
- **CLI オプション** — コマンドラインで明示的に指定されたフラグ（例: `--backend orb`）。*明示的に指定された*値のみが YAML を上書きし、CLI デフォルトのままのオプションは YAML の値を隠しません。

この優先順位は `all`, `opt`, `tsopt`, `freq`, `irc`, `scan`, `scan2d`, `scan3d`, `path-opt`, `path-search`, `dft` に共通です。あわせて {ref}`YAML リファレンス: 設定の優先順位 <ja-yaml-configuration-precedence>` を参照してください。

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

## 関連項目

- [インストール](installation.md) — セットアップと依存関係
- [はじめに](getting-started.md) — 初回実行とワークフロー概要
- [典型エラー別レシピ](recipes-common-errors.md) — 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) — よくあるエラーと対処法
- [YAML リファレンス](yaml-reference.md) — 全設定オプション
