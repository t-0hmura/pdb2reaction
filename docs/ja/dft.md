# `dft`（DFT 一点計算）

`dft` サブコマンドは、1 つの構造に対して GPU4PySCF（GPU）または PySCF（CPU）で **DFT（密度汎関数理論）一点計算**を行います。**エネルギー**と、Mulliken・meta-Löwdin・IAO（内在的原子軌道）のポピュレーション解析による**原子電荷**を出力します。デフォルトの手法は ωB97M-V/def2-svp です。DFT 用の追加パッケージが必要です。[詳細なインストール手順](installation.md#詳細なインストール手順)の手順 7 を参照してください。

`all --dft` は MLIP（機械学習原子間ポテンシャル）で求めた反応物（R）・遷移状態（TS）・生成物（P）に DFT 一点計算を行います。`-b dft` はコマンドのすべての計算を DFT で行います。違いは [MLIP の TS を DFT で確かめる](dft-backend.md) を参照してください。

## 主な用途

* **MLIP 構造での DFT エネルギー**: MLIP で最適化した R・TS・P の一点計算
* **電荷分布の把握**: 原子ごとの電荷と、開殻系のスピン密度
* **陰溶媒中のエネルギー**: `--solvent` による PySCF の PCM（分極連続体モデル）または SMD（密度に基づく溶媒和モデル）

---

## 基本的な実行例

### 1. GPU での一点計算

中性の一重項について、GPU でエネルギーと電荷を計算します。

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 --out-dir ./result_dft
```

端末に `E_total (Hartree): …` と `E_total (kcal/mol): …` が出て、`result_dft/result.yaml` に `energy.converged: true` があれば成功です。

### 2. SCF を厳しくし、基底を大きくする

SCF（自己無撞着場）の収束を厳しくし、基底を大きくします。

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 \
  --func-basis 'wb97m-v/def2-tzvpd' --scf-tol 1e-10 --scf-max-cycles 200 \
  --out-dir ./result_dft_tight
```

### 3. CPU だけで計算する

GPU の無いマシンでは、CPU の PySCF で計算できます。

```bash
pdb2reaction dft -i input.pdb -q 0 -m 1 --dft-engine cpu --out-dir ./result_dft_cpu
```

### 4. リガンドの電荷から総電荷を求める

`-q` を省略して `-l` でリガンドの形式電荷を与えると、`dft` は PDB 中のアミノ酸残基とイオンの電荷を足して総電荷を求め、その内訳を端末に表示します。

```bash
pdb2reaction dft -i input.pdb -l 'SAM:1,GPP:-3' -m 1 --out-dir ./result_dft_ligand
```

---

## 処理の仕組みと計算仕様

1. **構造の読み込み**:
PDB・mmCIF・XYZ・GJF を読み込み、PySCF に渡す座標を `input_geometry.xyz` に保存します。XYZ・GJF 入力では、`-l` に必要な PDB/mmCIF のトポロジーを `--ref-pdb` で与えます。`dft` は PDB・mmCIF・GJF のファイルを書き出しません。
2. **SCF**:
`--func-basis` で汎関数と基底を、`--dft-engine` で GPU4PySCF（`gpu`、デフォルト）か PySCF（`cpu`）を選びます。閉殻は RKS、開殻は UKS で計算します。デフォルトで有効な低メモリモードでは、密度フィッティングを使わずに J と K を直接組み立て、GPU の閉殻では GPU4PySCF の低メモリ版 RKS を使います。`--no-dft-low-memory` では密度フィッティングを使います。
3. **電荷と結果ファイル**:
SCF の後に Mulliken・meta-Löwdin・IAO の電荷とスピン密度を求め、Hartree と kcal/mol のエネルギーとともに `result.yaml` に書き出します。失敗した解析の列は `null` になり、警告が出ます。

---

## 主な出力ファイル

`--out-dir` に以下のファイルを書き出します。

```text
result_dft/
├─ input_geometry.xyz   # PySCF に渡した構造
├─ result.yaml          # エネルギー、収束、エンジン、原子ごとの電荷とスピン密度
├─ result.json          # 機械可読な要約（--out-json 指定時）
└─ summary.json         # result.json の写し。result.json を読む（--out-json 指定時）
```

* **`energy`**（`result.yaml`）: `hartree`・`kcal_per_mol`・`converged`・`used_gpu`・`used_lowmem`・`engine`。`engine` は `gpu4pyscf(rks_lowmem)`・`gpu4pyscf`・`pyscf(cpu)` のいずれかです。
* **`charges [index, element, mulliken, lowdin, iao]`**: 1 原子 1 行の表で、`index` は 0 始まりです。端末にも同じ表が出ます。
* **`spin_densities [index, element, mulliken, lowdin, iao]`**: 同じ形の表です。閉殻でも書き出し（値はすべて 0）、端末には開殻のときだけ表示します。
* **`result.json`**: 電荷とスピン密度を `mulliken`・`lowdin`・`iao` の配列で持ち、電荷・多重度・汎関数・基底・SCF の設定も記録します。[JSON 出力リファレンス](json-output.md#dft) を参照してください。

---

## 主な CLI オプション

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス | （必須） | 入力構造ファイル（`.pdb`, `.cif`, `.xyz`, `.gjf` 等） |
| `-q, --charge` | 整数 | `None` | 系全体の総電荷。`-l`、YAML の `calc.charge`、`.gjf` 入力のどれも無ければ必須 |
| `-m, --multiplicity` | 整数 | `1` | スピン多重度（2S+1）。`.gjf` 入力ではファイルの値を使用 |
| `-l, --ligand-charge` | 文字列 | `None` | 残基ごとの形式電荷（例: `'SAM:1,GPP:-3'`）またはリガンドの総電荷。PDB/mmCIF 入力か `--ref-pdb` が必要 |
| `--func-basis` | 文字列 | `wb97m-v/def2-svp` | 汎関数と基底（`汎関数/基底` の形） |
| `--scf-tol` | 浮動小数点数 | `1e-9` | SCF の収束閾値（Hartree） |
| `--scf-max-cycles` | 整数 | `100` | SCF の最大反復回数 |
| `--dft-grid-level` | 整数 | `3` | 数値積分グリッドのレベル（PySCF の `grids.level`） |
| `--dft-engine` | `gpu` / `cpu` | `gpu` | GPU4PySCF か CPU の PySCF |
| `--dft-low-memory/--no-dft-low-memory` | フラグ | `True` | J と K を直接組み立てる。`--no-dft-low-memory` で密度フィッティングを使用 |
| `--solvent` | 文字列 | `none` | PySCF の PCM/SMD に渡す溶媒名（例: `water`）。`none` は気相 |
| `--solvent-model` | `pcm` / `smd` | `smd` | 陰溶媒モデル |
| `--dft-nprocs` | 整数 | auto | PySCF の CPU スレッド数（スケジューラとホストから自動検出） |
| `--dft-memory` | 文字列 | auto | PySCF のホスト RAM の上限（例: `64GB`）。GPU メモリではない |
| `-o, --out-dir` | パス | `./result_dft/` | 出力先ディレクトリ |

全オプションの一覧は [自動生成 CLI リファレンス](../reference/commands/dft.md) を参照してください。

> **補足:** YAML（`--config`）では、{ref}`dft <ja-dft-section>` の節で同じ設定を指定できます。`dft.pyscf` は PySCF のオブジェクトに属性を名前で渡し、SCF が収束しにくいときは `pyscf: {mf: {level_shift: 0.2}}` のように使えます。電荷と多重度は、ほかのコマンドと同じく `calc.charge`・`calc.spin` に書きます。`calc.spin` は多重度 2S+1 で、PySCF の 2S ではありません。優先されるのは、コマンドラインで指定した `-q`・`-l`・`-m`、YAML、`.gjf` のヘッダーの順です。

---

(ja-notes)=
## 使用上の注意点

* **基底のコスト**: `def2-tzvpd` は `def2-svp` よりはるかに重い計算です。原子数や GPU メモリの決まった上限は無く、コストは基底関数の数・元素・汎関数・グリッド・GPU で決まります。まず代表構造を 1 つ計算し、メモリの最大使用量を確かめてください。足りないときは、基底を小さくするか、メモリの大きい GPU を使ってください。
* **GPU**: GPU4PySCF が動かないとき、`dft` は `--dft-engine cpu` を勧めるエラーで止まり、自動では CPU に切り替えません。新しい世代の GPU では、メモリ不足や未対応カーネルのエラーがメモリ量ではなく GPU4PySCF と CuPy の版から来ることがあるので、まず版とトレースバックを確かめてください。
* **CPU**: `--dft-engine cpu` で実用になる系の大きさは手法とマシンで変わるため、代表構造の一点計算で時間を測ってください。
* **一時ファイル**: PySCF は一時ファイルを `$PYSCF_TMPDIR` に書きます。`/tmp` の小さい計算ノードでは、実行前に空き容量の十分なディスクへ向けてください。
* **x86 以外のマシン**: GPU4PySCF のビルド済みホイールが対応しないことがあります。その場合は GPU4PySCF を[ソース](https://github.com/pyscf/gpu4pyscf)からビルドしてください。
* **補助基底**: `--no-dft-low-memory` では、選んだ基底に対する PySCF のデフォルトの補助基底を使います。自分で指定する必要はありません。
* **IAO 解析**は難しい系で失敗することがあります。
* **溶媒**: この `--solvent` は PySCF の PCM・SMD で、MLIP バックエンドの xTB 溶媒補正とは別物です。`--solvent-model` は小文字の `pcm` か `smd` だけを受け付けます。PCM で `dft` が知らない溶媒名を使うときは、YAML の `dft.pyscf.with_solvent.eps` に誘電率を指定してください。
* **SCF が収束しないとき**: `dft` は `WARNING: SCF did not converge to the requested tolerance.` を表示し、`converged: false` として `result.yaml` を書いた上で、終了コード 1 で終わります。低メモリモードでは、メモリに余裕があれば `--no-dft-low-memory` の別名 `--no-lowmem` で再実行するよう提案します。
* **多重度**: 1 未満は受け付けません。
* **前回の結果**: 実行の最初に、出力ディレクトリに残っている `result.yaml`・`result.json`・`summary.json` を削除します。
* **終了コード**: {ref}`終了コード <ja-exit-codes>`を参照してください。

---

## 関連ドキュメント

* [MLIP の TS を DFT で確かめる](dft-backend.md) — ワークフローでの `-b dft` と `--dft`、DFT の設定、GPU メモリ
* [sp](sp.md) — `-b dft` を含む任意のバックエンドでの一点エネルギーと力
* [all](all.md) — 全工程のワークフロー。`--dft` で R・TS・P に DFT 一点計算を追加
* [MLIP バックエンド](backends.md) — バックエンドの選び方
* [トラブルシューティング](troubleshooting.md) — 実行に失敗したときの対処
