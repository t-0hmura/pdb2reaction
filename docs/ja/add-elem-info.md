# `add-elem-info`

PDB ファイルの ATOM/HETATM レコードの元素記号カラム（77–78）を修復します。原子ごとの元素は原子名と残基コンテキストから推定し、`Bio.PDB.PDBParser` で再パースして原子単位で判定します。書き換えるのはカラム 77–78 のみです。PDB の元素カラム（77–78）が欠落・誤記で、`extract` / `opt` / `tsopt` などの下流サブコマンドが入力を受け付けないときに使用します。`all` はプリフライトで `add-elem-info` を自動呼び出しするため、手動実行は単独サブコマンドの前処理用です。

## 実行例

```bash
# 元素カラムを補完して "<input>_add_elem.pdb" に出力
pdb2reaction add-elem-info -i 1abc.pdb

# 出力ファイルを指定
pdb2reaction add-elem-info -i 1abc.pdb -o 1abc_fixed.pdb

# 入力ファイルをその場で上書き
pdb2reaction add-elem-info -i 1abc.pdb --overwrite
```

## 処理の流れ
1. `Bio.PDB.PDBParser` で入力を解析し、`extract.py` で使用される残基定義（`AMINO_ACIDS`、`WATER_RES`、`ION`）と同じ定義を使用して分類する。
2. 各原子について、原子名・残基名・ HETATM フラグを組み合わせて元素を推定:
 - `ION` 辞書に登録された単原子イオン残基は対応する元素を使用
 - タンパク質/核酸/水は H/D や Se の特例を扱い、C/N/O/P/S は先頭文字で判定（炭素側鎖ラベルは C）
 - その他のリガンドは原子名の接頭辞で判定し、ハロゲン認識や重水素→水素の正規化にフォールバック
3. `PDBIO` で書き出し:
 - デフォルト出力: `<input>_add_elem.pdb`（`-o/--out` 省略かつ `--overwrite` が `True` でない場合）
 - `-o/--out`: 指定パスへ書き込み（この場合 `--overwrite` は無視）
 - `--overwrite`（`-o/--out` なし）: 入力パスを上書き
4. 割り当て/再割り当て数、元素別合計、未解決原子のリスト（最大 50 件）を要約表示。

## 出力
- 元素記号が補完/修正された PDB ファイル:
 - デフォルト: `<input>_add_elem.pdb`（`-o/--out` 省略かつ `--overwrite` が `True` でない場合）
 - `-o/--out` 指定時: `OUTPUT.pdb`（`--overwrite` の値に関わらず）
 - `--overwrite` を `-o/--out` なしで指定: `INPUT.pdb` をその場で上書き
- コンソールに、処理/割り当て原子数、元素別カウント、未解決原子（最大 50 件）を出力

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 入力 PDB ファイル | 必須 |
| `-o, --out PATH` | 出力パス。指定した場合 `--overwrite` は無視される | _None_ → `<input>_add_elem.pdb` |
| `--overwrite/--no-overwrite` | `-o/--out` が省略された場合に入力を上書き | `False` |

完全なフラグ一覧は生成された [command reference](../reference/commands/index.md) を参照してください。

## 注記
- 変更されるのは列 77–78 のみ。座標、占有率、B 因子、電荷、altloc、挿入コード、レコード順序は保持されます。
- すべてのモデル/チェーン/残基にわたる ATOM/HETATM レコードを処理します。
- 重水素は水素に正規化され、セレン（`SE*`）やハロゲンは自動認識されます。
- すでに有効な元素記号を持つ PDB に対して再実行しても no-op（原子はそのまま通過）です。`pdb2reaction all` を実行すると、`add-elem-info` は **元素記号が欠落している PDB 入力に対してのみ、プリフライトステップとして自動的に呼び出されます**（その後 `fix-altloc`、続いて活性部位モデルの抽出が行われます）。元素記号がすでに補完されている入力はそのまま通過します。そのため `pdb2reaction all` の前に `add-elem-info` を手動で実行する**必要はありません**。手動実行が必要となるのは、元素カラムが欠落した PDB に対して `extract` や `opt` などの個別サブコマンドを直接使用する場合だけです。詳細は [all](all.md) を参照してください。

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
- [extract](extract.md) -- 元素カラム修正後の活性部位モデル抽出
- [all](all.md) -- 一気通貫ワークフローの入口
