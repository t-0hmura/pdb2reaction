# `add-elem-info`（PDB の元素欄の修復）

`add-elem-info` サブコマンドは、PDB ファイルの ATOM・HETATM レコードの**元素記号**（77–78 列）を**補い、または直します**。

---

## 主な用途

* **元素欄の無い PDB**: モデリングツールや分子動力学（MD）の出力で 77–78 列が空欄の構造
* **誤った元素記号**: `--overwrite-elem` で全原子の欄を推定し直す
* **ほかのコマンドの入力の準備**: `extract` などを単独で使う前に `add-elem-info` を実行してください

---

## 基本的な実行例

### 1. `<input>_add_elem.pdb` に書き出す

元素欄を補い、入力と同じ場所に書き出します。

```bash
pdb2reaction add-elem-info -i 1abc.pdb
```

端末に `[OK] Wrote: 1abc_add_elem.pdb` と、`total atoms`・`assigned/updated`・`kept existing` の件数が出ます。`[WARN]` の行が無ければ、すべての原子に元素が入っています。

### 2. 出力ファイルを指定する

```bash
pdb2reaction add-elem-info -i 1abc.pdb -o 1abc_fixed.pdb
```

### 3. 入力ファイルを上書きする

入力ファイルそのものを置き換えます。

```bash
pdb2reaction add-elem-info -i 1abc.pdb --overwrite
```

---

## 処理の仕組みと計算仕様

1. **レコードの読み込み**:
すべての行を読み、MODEL ブロックを追いながら、ATOM・HETATM レコードだけを対象にします。
2. **有効な記号の保持**:
元素欄にすでに有効な記号（水の仮想サイトでは `EP`）があれば、そのまま残します。空欄や認識できない記号の欄を修復し、`--overwrite-elem` ではすべての欄を推定し直します。
3. **元素の推定**:
4 文字の原子名（13–16 列）と残基名から元素を決めます。イオンの残基はそのイオンの元素、アミノ酸・核酸・水は通常の命名に従います。そのほかのリガンドは原子名が始まる列で判定し、14 列から始まる `NA` は N、13 列からなら Na です。14 列から始まる `CL1` などの LEaP のハロゲンや `HG11` などの水素の名前も扱います。
4. **書き出しと要約**:
修復したレコードの 77–78 列だけを変えて書き出し、端末に要約を表示します。

---

## 主な出力ファイル

* **修復した PDB**: デフォルトは `<input>_add_elem.pdb`。`-o` を指定するとそのパス、`-o` なしで `--overwrite` を指定すると入力ファイルそのものです。
* **端末の要約**: `total atoms`、`assigned/updated`（元素欄が変わった原子の数）、`kept existing`、`assignment breakdown`（元素ごとの件数）。割り当てられなかった原子は変更せず、`[WARN] Could not confidently assign N atoms; left unchanged.` の後に最大 50 件を表示します。列挙された原子は、77–78 列に元素記号を右詰めで手で書き込んでください。`--overwrite-elem` も同じ規則で推定するので、これらは割り当てられません。

---

## 主な CLI オプション

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス | （必須） | 入力 PDB ファイル |
| `-o, --output` | パス | `None` | 出力 PDB ファイル。省略時は `<input>_add_elem.pdb` |
| `--overwrite/--no-overwrite` | フラグ | `False` | `-o` を省略したときに入力ファイルを上書き |
| `--overwrite-elem/--no-overwrite-elem` | フラグ | `False` | 有効な記号が入っている元素欄も推定し直す |

全オプションは [自動生成のオプションの一覧（英語のみ）](../reference/commands/add_elem_info.md) を参照してください。

---

## 使用上の注意点

* **変わる部分**: 修復する ATOM・HETATM レコードの 77–78 列だけです。HEADER・REMARK・CONECT・ANISOU と電荷の欄（79–80 列）を含め、ほかの行はそのまま書き出します。
* **特殊な名前**: 重水素は H に、セレン（`SE*`）は Se になり、ハロゲンも自動で認識します。
* **2 つの別のフラグ**: `--overwrite-elem` はどの元素欄を推定し直すかを、`--overwrite` はファイルの書き出し先だけを決めます。
* **入力と同じパスへの書き出し**: `-o` が入力ファイルを指すとき（シンボリックリンク経由を含む）は `--overwrite` が必要で、無いとエラーで止まります。
* **終了コード**: {ref}`終了コード <ja-exit-codes>`を参照してください。

---

## 関連ドキュメント

* [extract](extract.md) — 修復した PDB から活性部位モデルを抽出
* [all](all.md) — 空欄の元素欄を自動で修復する全工程のワークフロー
* [トラブルシューティング](troubleshooting.md) — 実行に失敗したときの対処
