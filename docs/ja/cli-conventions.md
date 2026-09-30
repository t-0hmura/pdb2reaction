# CLI 規約

## ブール値オプション

オン・オフは `--flag` / `--no-flag` で指定します。

```bash
--tsopt --thermo --no-dft
```

旧来の値指定も、互換性のため受け付けます。

よく使うブール値オプション：
- `--tsopt`, `--thermo`, `--dft` — 後処理ステージの有効化
- `--freeze-links` — キャップ水素の親原子を凍結（デフォルト: `True`）
- `--dump` — 軌跡ファイルの出力
- `--preopt`, `--endopt` — 前処理/後処理最適化の切り替え
- `--climb` — MEP 探索でクライミングイメージを有効化
- `--convert-files` — 入力トポロジー／テンプレートに応じた PDB/CIF/GJF companion の生成

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
| `-v 2` | デフォルト。banner、`[command]`、`[mode]`、ステージ進捗、主要なオプティマイザのサイクル表、終了ステータス、Hessian 1 行要約、thermo / DFT 要約、経過時間を追加します。 |
| `-v 3` | デバッグ: resolved config、backend DEBUG、オプティマイザ・内部座標の詳細、`[HessianTiming]`、`[HessianVRAM]`。 |

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

### chain + 残基名による指定
```bash
-c 'A:SAM' # chain A 内の SAM をすべて選択
-c 'A:SAM:123' # chain A、残基名 SAM、残基番号 123 の1残基
```

### PDB/mmCIF ファイルによる指定
```bash
-c substrate.pdb # 別の PDB から座標を使用して基質を特定
-c substrate.cif # 別の mmCIF も使用可能
```

```{note}
残基名で選択する場合、同名の残基が複数あれば**すべて**が含まれ、警告がログに出力されます。
```

(ja-selected-resn-takes-ids)=
### `--selected-resn` も同じ残基selectorを使う

`extract` と `all` の `--selected-resn` は、残基ID・残基名・chain付き
残基名を受け付けます。`A:123A` は挿入コード付きの1残基、`A:SAM`
はchain A内の全SAM、`A:SAM:123` はそのうち1残基を強制包含します。
chainを付けない `TYR` は全matchを包含し、複数match時に警告します。

---

(ja-charge-specification)=
## 電荷の指定

PDB/mmCIF 入力では、`--ligand-charge/-l` を使うと**非標準残基（基質・補因子・金属イオンなど）の電荷だけ**を指定すれば、標準アミノ酸やイオンの電荷と合算して全系の電荷が自動計算されます。大きな酵素–基質系で総電荷を手動で数える必要がなくなります。

### 残基別マッピング（推奨）
```bash
-l 'SAM:1,GPP:-3' # SAM は +1、GPP は -3
-l 'LIG:-2' # LIG は -2
-l -3 # 単一の整数を渡すとリガンド総電荷を直接指定（マッピングの代替）
```

### 総電荷の明示指定
```bash
-q 0 # 総電荷を 0 に指定
-q -1 # 総電荷を -1 に指定
```

### 抽出を行わない場合の電荷解決順序
1. 明示的な `-q/--charge`
2. 明示的な `--ligand-charge/-l` と PDB/mmCIF 残基情報から導出した総電荷
3. どちらの CLI 電荷指定もない場合の `--config` 内 `calc.charge`
4. `.gjf` テンプレートのメタデータ
5. 未解決の場合は実行を中断

この順序は `opt`/`tsopt`/`freq` などの単独サブコマンドと、
`-c/--center` を省略した `all` に適用されます。

**`all -c/--center`:** 抽出処理は、標準残基、イオン、
`--ligand-charge/-l` を合算して総電荷を導出します。明示的な `-q` は
ここでも最優先で、系の総電荷を設定します。抽出値と異なる場合は警告を
表示します。`-q` を省略して `--ligand-charge/-l` を明示した場合は、
その指定を使った抽出由来の電荷が YAML より優先します。CLI の電荷指定が
どちらもない場合は、設定ファイルの `calc.charge` が自動抽出値を上書きし、
それもなければ抽出値を使用します。

```{tip}
非標準の残基（基質、補因子、特殊なリガンド）には必ず `--ligand-charge/-l` を指定し、電荷が正しく伝播するようにしてください。
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
--scan-lists '[("SAM,320,CS1", "GPP,321,C7", 1.60)]'
--scan-lists '[("A:SAM:320:CS1", "A:GPP:321:C7", 1.60)]' # chain付き
```

セレクタのフィールドは以下で区切れます。
- 空白: `'SAM 320 CS1'`
- カンマ: `'SAM,320,CS1'`
- スラッシュ: `'SAM/320/CS1'`
- バッククォート: `` 'SAM`320`CS1' ``
- バックスラッシュ: `'SAM\320\CS1'`

3つのtoken（残基名、残基番号、原子名）は任意の順序で指定できます。残基名や番号が重複するときは、位置固定の4-field形式`CHAIN:RESNAME:RESSEQ[ICODE]:ATOM`を使います。

---

(ja-scan-list-spec)=
### スキャンリスト仕様

`-s/--scan-lists`（`scan` / `scan2d` / `scan3d` / `all` で使用）は 1 個以上のインライン Python リテラルを受け付けます。スタンドアロンの `scan` / `scan2d` / `scan3d` はこれに加えて YAML/JSON スペックファイルパスも受け付けます。複雑な設定や複数ステージの実行にはファイルが、単純な 1 ステージのみの場合にはインラインリテラルが適しています。

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
- `scan` の各ステージは距離 target `(i, j, target_Å)` のリストです。距離・角度・二面角の range `(i,j,low,high)`・`(i,j,k,low,high)`・`(i,j,k,l,low,high)` も書けます。
- `scan2d`/`scan3d` の各軸は距離`(i,j,low,high)`、角度`(i,j,k,low,high)`、二面角`(i,j,k,l,low,high)`のいずれかです。
- indexは整数、3-field selector、または位置固定`CHAIN:RESNAME:RESSEQ[ICODE]:ATOM`で指定できます。

#### インライン Python リテラルの書式

各リテラルは Python リストです。シェルのクォート処理に注意が必要です。

```
-s '[(原子1, 原子2, ターゲット距離Å), ...]'      # scan: 3 要素タプル
-s '[(原子1, 原子2, 下限Å, 上限Å), ...]'         # 距離 range（scan では双方向。scan2d / scan3d の軸も同じ形）
```

- リスト全体を **シングルクォート** `'...'` で囲みます（シェルがカッコや空白を解釈しないようにするため）。
- `scan` では **1 リテラル = 1 ステージ**です。複数ステージを実行するには、**1 つの `-s/--scan-lists` フラグの後に複数リテラル**を並べます。
- `scan2d`/`scan3d` ではリテラルは **1 つだけ** を受け付けます（複数ステージは非対応）。`scan2d` ではちょうど2軸、`scan3d` ではちょうど3軸を含む必要があります。

| コマンド | 受け付けるスキャン仕様 |
| --- | --- |
| `scan` | インライン3-tupleまたは双方向4-tuple。YAML/JSONにも対応 |
| `all --scan-lists` | インラインの距離・角度・二面角の target tuple のみ（range と YAML/JSON は不可） |
| `scan2d` | 距離・角度・二面角の軸を2本含む1リテラル/ファイル |
| `scan3d` | 距離・角度・二面角の軸を3本含む1リテラル/ファイル |

##### 原子の指定方法

原子は**整数インデックス**または **PDB セレクタ文字列**で指定します。

| 方法 | 例 | 備考 |
| --- | --- | --- |
| 整数インデックス | `(1, 5, 2.0)` | デフォルトは 1 始まり（`--one-based`） |
| PDB セレクタ | `("SAM,320,CS1", "GPP,321,C7", 1.60)` | 残基名、残基番号、原子名の 3 要素タプル |
| chain付きselector | `("A:SAM:320:CS1", "A:GPP:321:C7", 1.60)` | 位置固定`CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` |

PDB セレクタのトークンは、カンマ `,`、スペース、スラッシュ `/`、バッククォート `` ` ``、バックスラッシュ `\` のいずれでも区切れます。トークンの順序も任意です。

```bash
# 以下はすべて同じ原子を指定:
"SAM,320,CS1"
"SAM 320 CS1"
"SAM/320/CS1"
"320,SAM,CS1" # 順序は自由
```

##### クォートの規則

```bash
# 正しい: シングルクォートでリストを囲み、セレクタはダブルクォート
-s '[("SAM,320,CS1","GPP,321,C7",1.60)]'

# 正しい: 整数インデックスなら内側のダブルクォートは不要
-s '[(1, 5, 2.0)]'

# 非推奨: ダブルクォートで外側を囲むとエスケープが必要
-s "[(\"SAM,320,CS1\",\"GPP,321,C7\",1.60)]"
```

## 入力ファイル要件

### PDB ファイル
- **水素原子**を含む必要があります（`reduce`、`pdb2pqr`、または Open Babel で追加）
- 列 77–78 に**元素記号**が必要（欠けている場合は `pdb2reaction add-elem-info` を使用）
- 複数の PDB は**同じ原子を同じ順序**で持つ必要があります（座標のみ異なる）

### mmCIFファイル
- multi-character chain、PDB欄を超えるresidue／atom ID、10,000残基以上では`.cif`／`.mmcif`を使用
- 共通bridgeが内部PDBで計算し、出力CIFへ元IDを復元
- 反応順に並べる複数入力は、形式変換後も同じ原子identityと順序が必要

### XYZ および GJF ファイル
- 活性部位モデル抽出をスキップする場合に使用可能（`-c/--center` を省略）
- `.gjf` ファイルは埋め込みメタデータから電荷/スピンのデフォルトを提供可能

---

(ja-exit-codes)=
## 終了コード

| コード | 意味 |
|---|---|
| `0` | 成功、または利用できる一部の結果あり |
| `1` | 未収束、利用できる結果なし、実行中の例外、出力失敗 |
| `2` | 入力・CLI 引数・設定の指定ミス |
| `130` | ユーザー中断（SIGINT） |

JSON の有無で終了コードは変わりません。IRC の上限停止だけでは失敗とせず、all は TS と端点の最適化を判定します。

(ja-opt-mode-semantics)=

## `--opt-mode`（サブコマンド依存）

```{warning}
同じ `--opt-mode` トークンでも、サブコマンドによって**選択される最適化アルゴリズムが異なり**、デフォルトも**統一されていません**。レシピをコピーする前に必ずサブコマンドごとの表を確認してください。
```

| サブコマンド | `grad` エイリアス | `hess` エイリアス | デフォルト |
|------------|------------------|------------------|-----------|
| `opt` | L-BFGS (`lbfgs`) | RFO (`rfo`) | `grad` (L-BFGS) |
| `tsopt` | Dimer (`dimer`) | RS-P-RFO (`rsprfo`) | `hess` (RS-P-RFO) |
| `path-opt`（端点 preopt） | L-BFGS | RFO | `grad` |
| `path-search`（HEI±1 / kink ノードの単一構造 optimizer） | L-BFGS | RFO | `grad` |
| `scan` / `scan2d` / `scan3d`（grid relaxation） | L-BFGS | RFO | `grad` |
| `all`（pre-opt 段階、`--opt-mode`） | L-BFGS | RFO | `grad` |
| `all`（post-opt — TSOPT プリセット、`--opt-mode-post`） | Dimer (`dimer`) | RS-P-RFO (`rsprfo`) | `hess` |
| `all`（post-opt — IRC 後エンドポイント最適化、`--opt-mode-post`） | L-BFGS | RFO | `hess` |

**受け付けるエイリアス**もサブコマンド固有です:

- `opt` は `grad` / `lbfgs` と `hess` / `rfo` を受け付けます。
- `tsopt` は `grad` / `dimer` と `hess` / `rsprfo` に加え、`rsirfo`（RS-I-RFO）、`trim`（TRIM/Helgaker）も単独の `--opt-mode` 値として受け付けます。
- `scan` / `scan2d` / `scan3d` / `path-opt` / `path-search` / `all` は `grad` / `hess` のみ受け付けます（アルゴリズム名 alias なし）。`all` の `--opt-mode-post` も `grad` / `hess` のみです。

したがって `tsopt` に対する `--opt-mode grad` は L-BFGS 最小化ではなく **Dimer TS 探索**です。曖昧さを避けたい場合は、各サブコマンドが受け付けるアルゴリズム名を使用してください: `opt` では `--opt-mode lbfgs|rfo`、`tsopt` では `--opt-mode dimer|rsirfo`。（他のサブコマンドは `grad` / `hess` のみ受け付けます。）

## CLI ↔ YAML 名称の不一致

一部の CLI フラグは YAML の対応キーと微妙に名前が異なり、`all` でラップされたときにリネームされるものもあります。完全なマッピング表は {ref}`YAML リファレンスの主要な CLI→YAML マッピング <ja-common-cli-to-yaml-mapping>` にあります。特に混同されやすい 2 ケースを以下に示します:

(ja-pressure-vs-pressure-atm)=
### `--pressure` (CLI) vs `pressure_atm` (YAML)

- **CLI フラグ:** `--pressure FLOAT`（`freq` サブコマンド; `all` では `--freq-pressure` として提供されます）。
- **YAML キー:** `thermo.pressure_atm`（単位接尾辞付き）。
- 両方とも値は **atm** 単位で扱われ、内部で Pa に変換されます。

(ja-engine-vs-dft-engine)=
### DFT engine

単体`dft`と`all`はどちらも`--dft-engine gpu|cpu`を使います。対応するYAML
keyは`dft.engine`です。{ref}`YAML リファレンスの dft セクション <ja-dft-section>`
も参照してください。

```bash
# 単体 dft
pdb2reaction dft -i ts.xyz -q 0 --dft-engine gpu

# all ラッパー内で同じ処理
pdb2reaction all -i r.pdb p.pdb -c SAM --tsopt --dft --dft-engine gpu
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
- **既知のデフォルト例外**: `flatten_max_iter` は YAML 適用前に 0 で初期化されます。toggle 未指定なら明示した YAML 値を保持し、`--flatten` は YAML／組み込みの正値を有効化、`--no-flatten` は 0 を強制します。{ref}`ja-flatten-precedence-caveat` を参照してください。

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
