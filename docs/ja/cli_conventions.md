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
`extract` と `fix-altloc` は parser wrapper（argparse バックエンド）ですが、root CLI で同じ bool 正規化が適用されます。

よく使うブール値オプション：
- `--tsopt`, `--thermo`, `--dft` — 後処理ステージの有効化
- `--freeze-links` — リンク水素の親原子を凍結（デフォルト: `True`）
- `--dump` — 軌跡ファイルの出力
- `--preopt`, `--endopt` — 前処理/後処理最適化の切り替え
- `--climb` — MEP 探索でクライミングイメージを有効化
- `--convert-files` — PDB/GJF コンパニオンファイルの生成

---

## Progressive Help (`all`)

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

---

## 電荷の指定

PDB 入力では、`--ligand-charge` を使うと**非標準残基（基質・補因子・金属イオンなど）の形式電荷（formal charge）だけ**を指定すれば、標準アミノ酸やイオンの電荷と合算して全系の電荷が自動計算されます。大きな酵素–基質系で総電荷を手動で数える必要がなくなります。

### 残基別マッピング（推奨）
```bash
--ligand-charge 'SAM:1,GPP:-3' # SAM は +1、GPP は -3
--ligand-charge 'LIG:-2' # LIG は -2
```

### 総電荷の明示的上書き
```bash
-q 0 # 総電荷を 0 に強制
-q -1 # 総電荷を -1 に強制
```

### 電荷の解決順序
1. `-q/--charge`（明示的な CLI 上書き）— 最優先
2. ポケット抽出（アミノ酸、イオン、`--ligand-charge` の合計）
3. フォールバックとしての `--ligand-charge`（抽出がスキップされた場合）
4. `.gjf` テンプレートのメタデータ
5. デフォルト: なし（未解決の場合は実行を中断します。`-q`、`.gjf` 電荷メタデータ、または PDB 入力の `--ligand-charge` のいずれかで解決してください）

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
- ポケット抽出をスキップする場合に使用可能（`-c/--center` を省略）
- `.gjf` ファイルは埋め込みメタデータから電荷/スピンのデフォルトを提供可能

---

## YAML 設定

詳細設定は多層 YAML で渡せます：

```bash
pdb2reaction -i r.pdb p.pdb -q -1 --config my_settings.yaml --out-dir result/
```

利用可能なすべてのオプションは [YAML リファレンス](yaml_reference.md) を参照してください。

### 設定の優先順位

設定は以下の順序で解決されます（後のものが前のものを上書き）：

```
組み込みデフォルト  <  --config (YAML)  <  CLI オプション  <  --override-yaml
```

- **組み込みデフォルト** — すべてのパラメータのハードコード値。
- **`--config`** — デフォルトを上書きする YAML ファイル。サイト共通やプロジェクト共通の設定に便利です。
- **CLI オプション** — コマンドラインで明示的に指定されたフラグ（例: `--backend orb`）。*明示的に指定された*値のみが YAML を上書きし、CLI デフォルトのままのオプションは YAML の値を隠しません。
- **`--override-yaml`** — 最後に適用される YAML レイヤー。実行ごとの微調整に使用します。

この優先順位は `all`, `opt`, `tsopt`, `freq`, `irc`, `scan`, `scan2d`, `scan3d`, `path-opt`, `path-search`, `dft` に共通です。

---

## 出力ディレクトリ

`--out-dir` で結果の保存先を指定します：

```bash
--out-dir ./my_results/ # カスタム出力ディレクトリ
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
- `dft`: `./result_dft/`

---

## 関連項目

- [インストール](installation.md) — セットアップと依存関係
- [はじめに](getting_started.md) — 初回実行とワークフロー概要
- [典型エラー別レシピ](recipes_common_errors.md) — 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) — よくあるエラーと対処法
- [YAML リファレンス](yaml_reference.md) — 全設定オプション
