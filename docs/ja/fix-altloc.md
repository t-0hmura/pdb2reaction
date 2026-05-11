# `fix-altloc`

## 概要

> **要約:** PDB ファイルから代替位置（altLoc）指示子を除去します。各原子について占有率が最も高いコンフォーマーを残し、列 17 をスペースに置換します。

### 要点
- **想定場面:** PDB に altLoc 指示子が残っていて、下流のパイプライン段階（`extract` / `opt` / `tsopt` など）が入力を受け付けない場合。クラスター抽出や構造最適化の前に実行します。
- **手法:** 各原子について占有率最高のものを優先選択（同点はファイル内出現順）し、保持側の列 17 のみをスペース 1 文字に置換。座標・B-factor などの他カラムは書き換えません。
- **主な出力:** クリーン化済 PDB を `<input>_clean.pdb`（ファイル入力）または `<input>_clean/`（ディレクトリ入力）に出力。`-o` で上書き先指定可。`--inplace` 指定で元ファイルを上書きします。
- **デフォルト値:** `--recursive False`、`--inplace False`、`--overwrite False`、`--force False`。altLoc 文字が無いファイルは既定でスキップされ、`--force` で強制処理に切り替わります。
- **次のステップ:** クリーン化された PDB を [`extract`](extract.md) や [`all`](all.md) に渡してください。下流段階は単一コンフォーマー入力を前提とします。

`fix-altloc` は、PDB ファイルから代替位置（altLoc）指示子を削除します。
各原子について占有率（occupancy）に基づいて最良のコンフォーマーを選択し、
重複を除去します。

### アルゴリズム
1. PDB の altLoc カラム（列 17、1-based）を空白（スペース 1 文字）に置換。
 - 1 文字の置換のみで、フォーマットのシフトや再構成は行いません。
2. 代替位置（A/B/... やカスタムラベル H/L など）により同じ原子が複数回
    出現する場合、「最良」のものを保持:
 - 占有率（occupancy）が最も高いものを優先
 - 同点の場合（または占有率が欠損している場合）、ファイル内で最初に出現したものを保持

### 処理対象レコード
- `ATOM` / `HETATM`: altLoc の選択とブランク化
- `ANISOU`: 対応する ATOM/HETATM レコード（同じシリアル番号）が保持される場合のみ保持

### 自動スキップ動作
デフォルトでは、ファイルに**altLoc 文字が含まれていない**場合（列 17 がすべて空白）、
ファイルは**スキップ**され、出力は書き込まれません。altLoc の有無に関わらず
処理を行うには `--force` を使用してください。

## 使用法
```bash
pdb2reaction fix-altloc -i INPUT.pdb [-o OUTPUT.pdb] [OPTIONS]
```

## 例
```bash
# 単一ファイルを処理（出力: INPUT_clean.pdb）
pdb2reaction fix-altloc -i 1abc.pdb

# 出力ファイルを指定
pdb2reaction fix-altloc -i 1abc.pdb -o 1abc_fixed.pdb

# ディレクトリを再帰的に処理
pdb2reaction fix-altloc -i ./structures -o ./cleaned --recursive

# 入力ファイルをその場で上書き（.bakバックアップを作成）
pdb2reaction fix-altloc -i ./structures --inplace --recursive

# altLocが検出されなくても強制的に処理
pdb2reaction fix-altloc -i 1abc.pdb -o 1abc_fixed.pdb --force
```

## ワークフロー
1. 入力ファイルに非空白の altLoc 文字（列 17）が含まれているかチェック。
 - altLoc が見つからず `--force` が設定されていない場合、ファイルをスキップ。
2. 各 ATOM/HETATM レコードについて、altLoc フィールドを無視した同一性キーを構築:
 - レコード名、原子名、残基名、チェーン ID、残基番号、挿入コード、segID
3. 同じ同一性キーを持つ原子の中から、以下の基準で選択:
 - 最高の占有率（列 55–60）
 - 同点の場合、ファイル内で最初に出現したもの
4. 出力を書き込み:
 - 選択された原子のみを保持
 - altLoc カラム（17）を空白（スペース 1 文字）に置換
 - ANISOU レコードは保持された原子に一致するもののみフィルタリング

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 入力 PDB ファイルまたはディレクトリ | 必須 |
| `-o, --out PATH` | 出力ファイル（入力がファイルの場合）またはディレクトリ（入力がディレクトリの場合） | 入力がファイル: `<input>_clean.pdb`、入力がディレクトリ: `<input>_clean/` |
| `--recursive/--no-recursive` | 入力がディレクトリの場合、`*.pdb` ファイルを再帰的に処理 | `False` |
| `--inplace/--no-inplace` | 入力ファイルをその場で上書き（`.bak` バックアップを作成） | `False` |
| `--overwrite/--no-overwrite` | 既存の出力ファイルの上書きを許可 | `False` |
| `--force/--no-force` | altLoc が検出されなくてもファイルを処理 | `False` |

## 出力
- 代替位置が削除された PDB ファイル:
 - 入力がファイル: デフォルトは `<input>_clean.pdb`（`-o/--out` が省略された場合）
 - 入力がディレクトリ: デフォルトは `<input>_clean/`（サブパスを保持）
 - `-o/--out` 指定時: `OUTPUT.pdb`
 - `--inplace` 設定時: 元のファイルを上書き（バックアップは `<input>.pdb.bak` として保存）

## `all` ワークフローとの統合
`pdb2reaction all` ワークフローを実行する際、`fix-altloc` は **`add-elem-info` の後**
（元素フィールドが欠損していた場合）、**活性部位モデル抽出の前**に自動的に実行されます。
これにより:
1. まず元素記号が補完される
2. 代替位置が単一のコンフォーマーに解決される
3. 抽出ステップがクリーンで曖昧さのない座標を受け取る

altLoc 文字が検出された場合のみファイルが処理され、そうでない場合は元のファイルが
そのまま渡されます。

## altLoc 状態間で原子数が異なる場合の処理

異なる altLoc 状態で異なる原子が含まれている場合（例：altLoc A には N, CA, CB, CG、
altLoc B には N, CA, CB, CD がある場合）、`fix-altloc` は以下のように正しく処理します：

- **重複原子**（複数の altLoc で同じ残基＋原子名、例：N, CA, CB）：
  占有率に基づいて最良のものを選択（最高値を優先、同点の場合はファイル内で最初のもの）
- **ユニーク原子**（1 つの altLoc にのみ存在、例：A の CG、B の CD）：
  **すべてのユニーク原子が出力に保持されます**

これにより、出力構造にはすべての altLoc 状態の原子が含まれ、
真の重複のみが単一のコンフォーマーに解決されます。

**例:**
```
入力:
 ATOM 1 N AALA A 1... 0.50 # altLoc A
 ATOM 2 CA AALA A 1... 0.50 # altLoc A
 ATOM 3 CG AALA A 1... 0.50 # altLoc A のみ
 ATOM 4 N BALA A 1... 0.40 # altLoc B
 ATOM 5 CA BALA A 1... 0.40 # altLoc B
 ATOM 6 CD BALA A 1... 0.40 # altLoc B のみ

出力:
 ATOM 1 N ALA A 1... 0.50 # A から（占有率が高い）
 ATOM 2 CA ALA A 1... 0.50 # A から（占有率が高い）
 ATOM 3 CG ALA A 1... 0.50 # 保持（Aのみ）
 ATOM 6 CD ALA A 1... 0.40 # 保持（Bのみ）
```

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- 原子シリアル番号は**再番号付けされません**（重複削除後にギャップが残る場合があります）。
- `CONECT` およびその他の接続/注釈レコードは**更新されません**。
- 変更されるのは列 17（altLoc）のみ。座標、占有率、B 因子、電荷、挿入コード、
  レコード順序は保持されます（重複削除を除く）。
- MODEL/ENDMDL ブロックは独立して処理されます。

## API 使用法
プログラム的な使用のために、モジュールは以下をエクスポートします:
```python
from pdb2reaction.fix_altloc import has_altloc, fix_altloc_file

# ファイルにaltLocがあるかチェック
if has_altloc(Path("input.pdb")):
 # altLocを修正
 was_processed = fix_altloc_file("input.pdb", "output.pdb", overwrite=True)
```
