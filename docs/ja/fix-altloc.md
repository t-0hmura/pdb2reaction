# `fix-altloc`

`fix-altloc` は、PDB ファイルから代替位置（altLoc）を解決します。
原子ごとではなく、残基ごとに一つの非空白 altLoc ラベルを選びます。
選択基準は、その残基内のラベル付き原子の平均 occupancy が最も高いこと、
同値の場合はファイル内で最初に現れることです。空白（共通）原子は残し、
非選択ラベルの原子は削除し、残したレコードの列 17 を空白にします。
構造workflowは同じ残基単位規則を共通input bridgeで自動適用します。この
commandはcleaned PDB自体が必要な場合、directory batch処理、または計算前の
conformer確認に使用します。

## 実行例

```bash
# コマンド形式
pdb2reaction fix-altloc -i INPUT.pdb [-o OUTPUT.pdb] [OPTIONS]

# 単一ファイルを処理（出力: INPUT_clean.pdb）
pdb2reaction fix-altloc -i 1abc.pdb

# 出力ファイルを指定
pdb2reaction fix-altloc -i 1abc.pdb -o 1abc_fixed.pdb

# ディレクトリを再帰的に処理
pdb2reaction fix-altloc -i ./structures -o ./cleaned --recursive

# 入力ファイルをその場で上書き（.bakバックアップを作成）
pdb2reaction fix-altloc -i ./structures --inplace --recursive
```

## 処理の流れ

1. 入力ファイルに非空白の altLoc 文字（列 17）が含まれているかチェック。
 - altLoc が見つからず `--force` が設定されていない場合、ファイルをスキップ。
2. ラベル付き ATOM/HETATM レコードを残基ごと（残基名、チェーン ID、
   残基番号、挿入コード、segID）に分ける。
3. 各残基で一つの非空白ラベルを選択:
 - ラベル内の解析可能な occupancy（列 55–60）の平均が最大
 - 同値または occupancy が無い場合、ファイル内で最初に現れるラベル
4. 空白（共通）原子と選択したラベルの原子を残し、残る重複を occupancy と
   出現順で解決する。
5. 出力を書き込み:
 - 空白（共通）原子と選択 conformer のみを保持
 - altLoc カラム（17）を空白（スペース 1 文字）に置換
 - ANISOU レコードは保持された原子に一致するもののみフィルタリング

### 処理対象レコード
- `ATOM` / `HETATM`: altLoc の選択とブランク化
- `ANISOU`: 対応する ATOM/HETATM レコード（同じシリアル番号）が保持される場合のみ保持

### altLoc 状態間で原子数が異なる場合の処理

異なる altLoc 状態で原子集合が異なる場合（例：A に N, CA, CB, CG、B に
N, CA, CB, CD）も、選択した残基ラベルの原子のみを残します。非選択ラベルに
しか無い原子は削除されます。これにより、A/B を混合した実在しないハイブリッド
残基の生成を防ぎます。

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
```

## 出力
- 代替位置が削除された PDB ファイル:
 - 入力がファイル: デフォルトは `<input>_clean.pdb`（`-o/--out` が省略された場合）
 - 入力がディレクトリ: デフォルトは `<input>_clean/`（サブパスを保持）
 - `-o/--out` 指定時: `OUTPUT.pdb`
 - `--inplace` 設定時: 元のファイルを上書き（バックアップは `<input>.pdb.bak` として保存）

## Python API
プログラムから利用する場合、モジュールは以下をエクスポートします:
```python
from pdb2reaction.io.pdb_fix import has_altloc, fix_altloc_file

# ファイルにaltLocがあるかチェック
if has_altloc(Path("input.pdb")):
 # altLocを修正
 was_processed = fix_altloc_file("input.pdb", "output.pdb", overwrite=True)
```

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 入力 PDB ファイルまたはディレクトリ | 必須 |
| `-o, --out PATH` | 出力ファイル（入力がファイルの場合）またはディレクトリ（入力がディレクトリの場合） | 入力がファイル: `<input>_clean.pdb`、入力がディレクトリ: `<input>_clean/` |
| `--recursive/--no-recursive` | 入力がディレクトリの場合、`*.pdb` ファイルを再帰的に処理 | `False` |
| `--inplace/--no-inplace` | 入力ファイルをその場で上書き（`.bak` バックアップを作成） | `False` |
| `--overwrite/--no-overwrite` | 既存の出力ファイルの上書きを許可 | `False` |
| `--force/--no-force` | altLoc が検出されなくてもファイルを処理 | `False` |

すべてのフラグの一覧は、生成される[コマンドリファレンス](../reference/commands/index.md)にあります。

## 注記
- デフォルトでは、ファイルに**altLoc 文字が含まれていない**場合（列 17 がすべて空白）、ファイルは**スキップ**され、出力は書き込まれません。altLoc の有無に関わらず処理を行うには `--force` を使用してください。
- 原子シリアル番号は**再番号付けされません**（重複削除後にギャップが残る場合があります）。
- `CONECT` およびその他の接続/注釈レコードは**更新されません**。
- 生存したレコードの座標、occupancy、B 因子、電荷、挿入コード、相対順は保持します。
  非選択 altLoc レコードは削除し、生存レコードの列 17 は空白にします。
- MODEL/ENDMDL ブロックは独立して処理されます。
- 残基単位の occupancy 規則もヒューリスティックです。活性部位の接触や機構に基づいて
  conformer を選ぶ必要がある場合は、構造エディタで明示選択し、目視確認してください。
- `all` とstandalone geometry commandは共通bridgeでaltLocを残基単位に解決します。`fix-altloc` commandはcleaned PDBを別fileとして保持したい場合に使います。

## 関連項目
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。
- [`extract`](extract.md) — 共通bridgeでaltLocを自動解決する活性部位model抽出。
- [`all`](all.md) — 共通altLoc/大規模構造bridgeを使うend-to-end workflow。
