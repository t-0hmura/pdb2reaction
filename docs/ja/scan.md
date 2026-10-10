# `scan`（拘束付き座標スキャン）

`scan` サブコマンドは、1 つの構造の中の距離・角度・二面角を調和拘束で少しずつ動かし、各点でそれ以外の自由度を緩和して、1 つの構造から反応経路の候補を作ります。1 つのリテラルや YAML の 1 つのステージに書いた座標は 1 つの**ステージ**として一緒に動きます。リテラルを複数並べるとステージが順に実行され、各ステージは前のステージの緩和後の構造から始まります。

---

## 主な用途

* **1 つの構造からの経路づくり**: 反応物の反応する結合を動かして、中間体や生成物に近い構造を作り、[`path-search`](path-search.md) に渡す
* **反応の順序の検討**: 結合形成とプロトン移動を 1 つのステージで動かす場合と、別のステージに分ける場合とで、エネルギーの変化を比べる
* **`all` のスキャン段の単独実行**: [`all`](all.md) が `-s` で行うスキャンを、刻み幅や拘束を変えて単独で実行し直す

計算バックエンドにはデフォルトの **UMA**（Meta）のほか、`-b/--backend` オプションで **ORB**、**MACE**、**AIMNet2**、DFT（`dft`）も選択可能です。独立した 2 つまたは 3 つの座標でエネルギーの格子を作るには、[`scan2d`](scan2d.md) または [`scan3d`](scan3d.md) を使います。

---

## 基本的な実行例

例の `input.pdb` は、同梱の酵素構造から [extract](extract.md) で切り出したクラスターモデルで、電荷は `-l 'SAM:1,GPP:-3'` から求めます。

```bash
pdb2reaction extract -i examples/1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' -o input.pdb
```

この PDB は chain の欄が空なので、原子は残基名・残基番号・原子名の 3 項目を任意の順序で、カンマか空白で区切って `'SAM,320,CS1'` や `'CS1 SAM 320'` のように書きます。

### 1. YAML スペックファイルからの実行

ステージをファイルに書き、`--out-json` を付けて結果の要約も出力します。

```yaml
# scan.yaml
stages:
  - [["SAM,320,CS1", "GPP,321,C7", 1.60]]
  - [["GPP,321,H11", "GLU,186,OE2", 0.90]]
```

```bash
pdb2reaction scan -i input.pdb -l 'SAM:1,GPP:-3' -m 1 -s scan.yaml --out-json -o ./result_scan
```

端末には各ステージの共有結合の変化が出て、最後に `====== Scan summary ======` が出ます。`result_scan/result.json` には `scientific_status` が入り、`success` でないときの理由は `scientific_status_reasons` に出ます。

### 2. インラインリテラルでの指定

単純な 1 ステージのスキャンは、コマンドラインに直接書けます。

```bash
pdb2reaction scan -i input.pdb -l 'SAM:1,GPP:-3' -m 1 -s '[("SAM,320,CS1","GPP,321,C7",1.60)]'
```

### 3. 2 つの座標を 1 つのステージで動かす

同じリテラルの中の座標は一緒に（協奏的に）動きます。

```bash
pdb2reaction scan -i input.pdb -l 'SAM:1,GPP:-3' -m 1 \
    -s '[("CS1 SAM 320","GPP 321 C7",1.60),("GPP 321 H11","GLU 186 OE2",0.90)]' -o ./result_concerted
```

### 4. 2 つのステージを順に実行する

1 つの `-s` の後にリテラルを複数並べます。ステージ 2 はステージ 1 の緩和後の構造から始まります。

```bash
pdb2reaction scan -i input.pdb -l 'SAM:1,GPP:-3' -m 1 \
    -s '[("SAM,320,CS1","GPP,321,C7",1.60)]' '[("GPP,321,H11","GLU,186,OE2",0.90)]' -o ./result_staged
```

### 5. 双方向スキャン

[4-tuple](#双方向スキャン4-tuple) を使うと、1 つの距離を入力構造から両方向にスキャンします。

```bash
pdb2reaction scan -i input.pdb -l 'SAM:1,GPP:-3' -m 1 -s '[("SAM,320,CS1","GPP,321,C7",1.60,3.00)]'
```

### 6. 軌跡の保存

`--dump` を付けると、各ステップの最適化の軌跡も保存します。

```bash
pdb2reaction scan -i input.pdb -l 'SAM:1,GPP:-3' -m 1 -s scan.yaml --dump -o ./result_scan_dump
```

---

## 処理の仕組みと計算仕様

1. **構造の読み込み**:
{ref}`電荷 <ja-charge-specification>` は `-q` または `-l` から決まります。`--preopt` を付けると、まず拘束なしで構造を最適化します。収束しなかった場合は入力構造を使います。
2. **ステージのステップ分割**:
座標ごとに変化量 Δ = 目標値 − 現在値 を求め、ステージを N = ceil(max(|Δ| / h)) ステップに分けます。h は距離では `--max-step-size`（Å）、角度では `--max-angle-step-size`、二面角では `--max-dihedral-step-size`（度）です。各座標は 1 ステップに Δ / N ずつ動くので、ステージ内のすべての座標が同時に目標値に着きます。
3. **拘束付きの緩和**:
各ステップで、調和拘束 E = ½ k (q − q_target)² がスキャンする座標 q をそのステップの目標値に保ち（k は `--restraint-k`）、残りの構造を `--opt-mode grad`（デフォルト）では L-BFGS、`hess` では RFO で緩和します。各ステップのエネルギーは、拘束を外して計算した値を記録します。
4. **ステージの終わり**:
`--endopt` を付けると、ステージの最後の構造を拘束なしでもう一度最適化します。そのあと、ステージの最初と最後の構造を比べて共有結合の変化を調べ、ステージの結果を書き出します。
5. **次のステージ**:
次のステージはこの結果から始まります。最後のステージが終わると、全ステージの軌跡を 1 つのファイルにつなぎます。

### 双方向スキャン（4-tuple）

目標値 `(i, j, target)` の代わりに範囲 `(i, j, low, high)` を指定すると、入力構造から両方向にスキャンします。範囲は 2 つのステージに展開されます。

1. **1 回目**: `i`–`j` の距離を現在の値から `low` に向けて動かす。
2. **2 回目**: 入力構造に戻し、`i`–`j` の距離を `high` に向けて動かす。

つないだ軌跡は `low → 入力構造 → high` の順になり、出発構造を通る連続した経路になります。角度の範囲 `(i, j, k, low, high)` と二面角の範囲 `(i, j, k, l, low, high)` も同じようにスキャンします。

(ja-section-bond)=
### 結合変化の検出

両原子の共有結合半径の和に `bond_factor`（デフォルト `1.20`）を掛けた値を T とします。2 原子の距離が T − `margin_fraction` × T（デフォルト `0.05`）以下なら、結合しているとみなします。結合の形成・切断として報告するのは、距離が `delta_fraction` × T（デフォルト `0.05`）以上変わった組だけです。`path-search` も同じ基準を使います。キーは YAML の [`bond`](yaml-reference.md#bond) の節にあります。

---

(ja-scan-checking-result)=
## 結果の判定

| 確認する場所 | 見るもの |
| --- | --- |
| 端末（各ステージ） | `[stage k] Covalent-bond changes (start vs final): Yes` とできた結合・切れた結合の一覧、または `No` と `(no covalent changes detected)` |
| 端末（実行の最後） | `====== Scan summary ======`：各ステージの目標値・ステップ数・結合変化 |
| `result.json`（`--out-json`） | `scientific_status`：すべてのステージの全ステップが収束し（`--endopt` を付けたときはその最適化も収束し）、エネルギーが有限なら `success`、一部のステージだけなら `partial`、1 つも無ければ `failed` |
| `result.json`（`--out-json`） | `stages[].converged`、`stages[].bond_changes.changed`、`stages[].final_energy_hartree`、各ステップのエネルギー `stages[].energies_hartree` |

`partial` の終了コードは 0、`failed` は 1 です。収束しなかったステージについては {ref}`max_cycles とプラトー停止 <ja-troubleshooting-max-cycles>` を参照してください。収束して狙った結合変化が起きたスキャンは経路の候補になり、エネルギーが最も高いステップは [`tsopt`](tsopt.md) に渡す TS 候補になります。このステップは `scan_trj.xyz` から {ref}`取り出せます <ja-trajectory-one-frame>`。

---

(ja-scan-direction-barrier-sign)=
## 障壁の向き

`scan` はエネルギーを記録しますが、障壁は出力しません。スキャンから障壁を読む場合、順方向の障壁は常に反応物から計算します。

| 実行内容 | 開始構造との差 | 順方向障壁 |
| --- | --- | --- |
| 反応物から始めたスキャン | `E(TS) − E(reactant)` | 開始構造との差と同じ |
| 生成物から始めたスキャン | `E(TS) − E(product)`。**逆方向**の障壁 | `E(TS) − E(reactant)`。開始構造との差では**ない**。E(reactant) は最適化した反応物のエネルギー（例: [`opt`](opt.md) で最適化した IRC の端点） |

これを切り替えるオプションはありません。障壁を引用する前に、スキャンがどちらの端点から始まったかを確認してください。結晶構造の生成物複合体から始めた場合は特に注意してください。

---

## 主な出力ファイル

`--out-dir` に次のファイルを書きます。

```text
result_scan/
├─ preopt/
│  └─ result.xyz                    # 事前最適化した構造（--preopt 指定時）
├─ stage_01/                        # ステージごとのディレクトリ（stage_NN）
│  ├─ result.xyz                    # ステージの final geometry
│  ├─ scan_trj.xyz                  # ステージ内の各ステップの構造とエネルギー
│  └─ scan_s0001_optimization_trj.xyz  # 各ステップの最適化の軌跡（--dump 指定時）
├─ scan_trj.xyz                     # 全ステージをつないだ軌跡
└─ result.json                      # 結果の要約（--out-json 指定時）。summary.json も同じ内容
```

PDB・mmCIF 入力では、各 `result.xyz` を `result.pdb`、各 `scan_trj.xyz` を `scan.pdb` として同じディレクトリにも書き、Gaussian 入力では構造を `result.gjf` でも書きます。{ref}`mmCIF の入力 <ja-mmcif-input>` と、PDB の欄に入りきらない大きな PDB の入力では、元の識別子を保った `.cif` も書きます。

* **ステージの結果**: `stage_NN/result.*` はステージ NN の終わりの構造です。[`path-search`](path-search.md) には、`all` と同じく、開始構造に続けて `stage_NN/result.*` をステージの順に渡します。`--preopt` を付けたときの開始構造は `preopt/result.*` です。
* **エネルギーの変化**: `scan_trj.xyz` の各フレームのコメント行には、拘束を外したエネルギー（Hartree）が入っています。[`trj2fig`](trj2fig.md) で図にできます。

---

## 主な CLI オプション

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス | （必須） | 入力構造ファイル（`.pdb`, `.cif`, `.mmcif`, `.xyz` 等） |
| `-q, --charge` | 整数 | `None` | 系全体の総電荷。`-l` を使う場合と `.gjf` 入力のほかは必須 |
| `-m, --multiplicity` | 整数 | `1` | スピン多重度（2S+1） |
| `-l, --ligand-charge` | 文字列 | `None` | リガンドの総電荷（例: `-1`）または残基名ごとの電荷（例: `'GPP:-3,SAM:1'`）。`-q` を省いたときに使用（PDB/mmCIF 入力または `--ref-pdb`） |
| `-s, --scan-lists` | 文字列 | （必須） | YAML/JSON スペックファイル、または 1 つ以上のインラインリテラル（1 つが 1 ステージ）。距離の目標値 `(i,j,target)`、または距離 `(i,j,low,high)`・角度 `(i,j,k,low,high)`・二面角 `(i,j,k,l,low,high)` の範囲 |
| `-o, --out-dir` | パス | `./result_scan/` | 出力先ディレクトリ |
| `--one-based/--zero-based` | フラグ | `--one-based` | `-s` の原子インデックスを 1 始まり / 0 始まりとして読む |
| `--max-step-size` | 浮動小数点数 | `0.2` | 1 ステップあたりの距離の最大変化量（Å） |
| `--max-angle-step-size` | 浮動小数点数 | `5.0` | 1 ステップあたりの角度の最大変化量（度） |
| `--max-dihedral-step-size` | 浮動小数点数 | `10.0` | 1 ステップあたりの二面角の最大変化量（度） |
| `--restraint-k` | 浮動小数点数 | `300.0` | 拘束の強さ k（距離は eV/Å²、角度は eV/rad²）。別名 `--bias-k` |
| `--preopt/--no-preopt` | フラグ | `False` | スキャンの前に入力構造を拘束なしで最適化 |
| `--endopt/--no-endopt` | フラグ | `False` | 各ステージの結果を拘束なしで最適化 |
| `--dump/--no-dump` | フラグ | `False` | 各ステップの最適化の軌跡を出力 |
| `--opt-mode` | `grad` / `hess` | `grad` | 緩和の方法：L-BFGS / RFO（`tsopt` では同じ語が別の最適化法を指す。{ref}`コマンドごとの --opt-mode <ja-opt-mode-semantics>` を参照） |
| `--freeze-links/--no-freeze-links` | フラグ | `True` | クラスター境界のキャップ水素の親原子を自動固定 |
| `--out-json/--no-out-json` | フラグ | `False` | 結果の要約を `result.json` に出力（[JSON 出力の一覧](json-output.md)） |

全オプションは [自動生成のオプションの一覧（英語のみ）](../reference/commands/scan.md) を参照してください。

---

## 使用上の注意点

* **`--preopt` は呼び出し方で変わる**: `scan` を単独で実行したときは、`--preopt` を付けたときだけ事前最適化をします。`all` の中では `all --preopt`（デフォルトで有効）に従い、`all --scan-preopt/--no-scan-preopt` で上書きできます。
* **`all -s` のタプル**: {ref}`スキャンリスト仕様 <ja-scan-list-spec>` の `scan` の形で書きます。
* **インラインでは目標値と範囲を混ぜない**: 1 つのインラインリテラルの中でも、1 回の実行のリテラルどうしでも、目標値 `(i,j,target)` と範囲のどちらか一方だけを使います。両方を組み合わせるときは、YAML/JSON スペックの `stages:` に並べてください。
* **範囲を使うときのステージ番号**: 範囲 1 つは `low` 向きと `high` 向きの 2 つのステージになります。インラインでは、1 つのリテラルの範囲がすべてこの 2 つのステージを共有します。YAML の `stages:` では、目標値だけのステージは 1 つのままで、範囲を含むステージは、目標値 1 つにつき 1 つ、範囲 1 つにつき 2 つのステージに分かれます。
* **目標の距離は正の値**にしてください。また、1 つのステージに同じ座標を 2 回書くことはできません。
* **計算せずに指定を確かめる**: `--dry-run` は入力・電荷とスピン・`-s` を読み、ステージの数を表示して、最適化をせずに終了します。
* **サイクル数の上限**: `--relax-max-cycles`（デフォルト `100000`）が各緩和のサイクル数を制限します。指定すると YAML の `opt.max_cycles` より優先されます。
* **YAML での拘束の強さ**: `--config` では、`--restraint-k` を指定しないときに {ref}`bias.k <ja-bias-section>` が使われます。

---

## 関連ドキュメント

* {ref}`スキャンリスト仕様 <ja-scan-list-spec>` — YAML/JSON スペックファイル、インラインリテラル、原子の指定
* [scan2d](scan2d.md) — 2 つの座標のエネルギーマップ
* [scan3d](scan3d.md) — 3 つの座標のエネルギー格子
* [path-search](path-search.md) — スキャンの結果からの最小エネルギー経路（MEP）探索
* [all](all.md) — 1 つの構造と `-s` からのスキャンを含む一気通貫ワークフロー
* [トラブルシューティング](troubleshooting.md) — 異常終了時の原因切り分けと対処法
