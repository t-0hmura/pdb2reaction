# `add-elem-info`

## 概要
`add-elem-info` は、PDBファイルの ATOM/HETATM レコードにある元素記号カラム（77–78）を修復します。

### 出力の挙動
- `-o/--out` が**省略**され、`--overwrite` が **`True` でない**場合、出力は `<input>_add_elem.pdb` に書き込まれます（末尾の `.pdb` を `_add_elem.pdb` に置換）。
- `--overwrite True` **かつ** `-o/--out` が**省略**されている場合、**入力ファイルをその場で上書き**します。`-o/--out` が指定された場合、`--overwrite` は無視されます。

## 使用法
```bash
pdb2reaction add-elem-info -i INPUT.pdb [-o OUTPUT.pdb] [--overwrite {True\|False}]
```

## 例
```bash
# 元素カラムを補完して "<input>_add_elem.pdb" に出力
pdb2reaction add-elem-info -i 1abc.pdb

# 出力ファイルを指定
pdb2reaction add-elem-info -i 1abc.pdb -o 1abc_fixed.pdb

# 入力ファイルをその場で上書き
pdb2reaction add-elem-info -i 1abc.pdb --overwrite True
```

## ワークフロー
1. `Bio.PDB.PDBParser` で入力を解析し、`extract.py` で使用する残基定義（`AMINO_ACIDS`、`WATER_RES`、`ION`）に合わせる。
2. 各原子について、原子名・残基名・HETATMフラグを組み合わせて元素を推定:
   - `ION` 辞書にある単原子イオン残基は対応元素を使用
   - タンパク質/核酸/水はH/DやSeの特例を扱い、C/N/O/P/Sは先頭文字で判定（炭素側鎖ラベルはC）
   - それ以外のリガンドは原子名の接頭辞で判定し、ハロゲン認識や重水素→水素の正規化でフォールバック
3. `PDBIO` で書き出し:
   - デフォルト出力: `<input>_add_elem.pdb`（`-o/--out` 省略かつ `--overwrite` が `True` でない場合）
   - `-o/--out`: 指定パスへ書き込み（この場合 `--overwrite` は無視）
   - `--overwrite True`（`-o/--out` なし）: 入力パスを上書き
4. 割り当て/再割り当て数、元素別合計、未解決原子のリスト（最大50件）を要約表示。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 入力PDBファイル | 必須 |
| `-o, --out PATH` | 出力パス。指定した場合 `--overwrite` は無視される | _None_ → `<input>_add_elem.pdb` |
| `--overwrite {True\|False}` | `-o/--out` が省略された場合に入力を上書き | `False` |

## 出力
- 元素記号が補完/修正されたPDBファイル:
  - デフォルト: `<input>_add_elem.pdb`（`-o/--out` 省略かつ `--overwrite` が `True` でない場合）
  - `-o/--out` 指定時: `OUTPUT.pdb`（`--overwrite` の値に関わらず）
  - `--overwrite True` を `-o/--out` なしで指定: `INPUT.pdb` をその場で上書き
- コンソールに、処理/割り当て原子数、元素別カウント、未解決原子（最大50件）を出力

## 注意事項
- 変更されるのは列 77–78 のみ。座標、占有率、B因子、電荷、altloc、挿入コード、レコード順序は保持されます。
- すべてのモデル/チェーン/残基に渡る ATOM/HETATM レコードを処理します。
- 重水素は水素に正規化され、セレン（`SE*`）やハロゲンは自動認識されます。