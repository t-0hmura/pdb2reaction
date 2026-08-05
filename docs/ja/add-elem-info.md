# `add-elem-info`

PDB の ATOM/HETATM record の元素記号列（77–78）を修復します。固定列の原子名と残基情報から元素を推定します。`all` の preflight が修復するのは空欄の元素列だけで、非空欄の値は化学的に検証せず保持します。誤った非空欄記号を直す場合、または空欄のある入力を単独 command へ渡す前に、この utility を明示実行してください。

## 実行例

```bash
# 元素列を補完して "<input>_add_elem.pdb" に出力
pdb2reaction add-elem-info -i 1abc.pdb

# 出力ファイルを指定
pdb2reaction add-elem-info -i 1abc.pdb -o 1abc_fixed.pdb

# 入力ファイルをその場で上書き
pdb2reaction add-elem-info -i 1abc.pdb --overwrite
```

## 処理の流れ
1. PDB の生レコードを読み、`extract.py` で使用される残基定義（`AMINO_ACIDS`、`WATER_RES`、`ION`）と同じ定義に従って分類する。
2. 各原子について、原子名と残基名を組み合わせて元素を推定:
 - `ION` 辞書に登録された単原子イオン残基は対応する元素を使用
 - タンパク質/核酸/水は H/D や Se の特例を扱い、C/N/O/P/S は先頭文字で判定（炭素側鎖ラベルは C）
 - その他のリガンドは原子名の接頭辞で判定し、ハロゲン認識や重水素→水素の正規化にフォールバック
3. ATOM/HETATM の列 77–78 だけを置換し、その他の列とレコードを変更せずに書き出し:
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
- 修復対象の ATOM/HETATM レコードの列 77–78 を除き、各入力行はそのまま保持されます。HEADER/REMARK/CONECT/ANISOU と従来形式の電荷列（79–80）も保持されます。
- すべてのモデル/チェーン/残基にわたる ATOM/HETATM レコードを処理します。
- 重水素は水素に正規化され、セレン（`SE*`）やハロゲンは自動認識されます。
- inference では既存の非空欄元素記号を保持しますが、writer が record formatting を正規化する場合があります。`all` は元素記号が空欄の PDB input にだけ preflight を適用します。詳細は [all](all.md) を参照してください。

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
- [extract](extract.md) -- 元素列修正後の活性部位モデル抽出
- [all](all.md) -- 一気通貫ワークフローの入口
