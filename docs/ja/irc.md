# `irc`（固有反応座標）

`irc` サブコマンドは、最適化した遷移状態（TS）から、EulerPC（Euler 予測子–修正子法）で固有反応座標（IRC）を両方向へたどります。各分岐の軌跡と、2 つの端点の候補を書き出します。この端点を [`opt`](opt.md) で最適化すると、TS がどの反応物（R）と生成物（P）をつなぐかが分かります。

---

## 主な用途

* **TS の確認**: [`tsopt`](tsopt.md) で得た TS が、意図した R と P をつなぐかを確かめる
* **R と P の取得**: 端点を [`opt`](opt.md) で最適化し、この TS の R と P の構造を得る
* **`all` の IRC 段のやり直し**: [`all`](all.md) の IRC を、設定を変えて単独でたどり直す

デフォルトの計算バックエンドは、Meta が公開した学習済みの[機械学習原子間ポテンシャル（MLIP）](backends.md)の **UMA** です。`-b/--backend` で **ORB**、**MACE**、**AIMNet2**、**DFT** も選べます。

---

## 基本的な実行例

### 1. 両方向の IRC

TS から両方向へたどり、`--out-json` で結果の要約も書き出します。

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --out-json --out-dir ./result_irc
```

端点の候補は `finished_first.xyz` と `finished_last.xyz` です。

### 2. 順方向だけ、大きいステップ

順方向の分岐だけを、最大ステップ 0.2 bohr でたどります。

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --no-backward --step-size 0.2 --out-dir ./result_irc_forward
```

### 3. 解析 Hessian

最初の Hessian を、有限差分ではなく解析的に計算します。

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --hessian-calc-mode Analytical --out-dir ./result_irc_analytical
```

### 4. 小さいステップでの再試行

分岐が数フレームで止まるときは、最大ステップを 0.05 bohr にして再試行します。

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --step-size 0.05 --out-dir ./result_irc_small_step
```

### 5. サイクルの上限までたどる

`--never-stop` を付けると、勾配とエネルギーによる停止の条件を無視し、各分岐を `--max-cycles` までたどります。

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --step-size 0.05 --never-stop \
    --max-cycles 250 --out-dir ./result_irc_continue
```

---

## 処理の仕組みと計算仕様

1. **出発の方向**: TS で Hessian を計算するか `--read-hess` のファイルから読み、剛体運動を [`freq`](freq.md#固定境界での剛体モード) と同じように除いてから、`--root` 番目の固有ベクトル（`0` が最小の固有値）を反応モードとします。そのモードが虚振動でなければ、エラーで止まります。
2. **EulerPC による積分**: 各分岐（順方向、次に逆方向）は TS から始まります。各ステップでは、質量加重の最急降下方向に沿って Euler 予測子で進みます。予測子の勾配は、Bofill 式で更新する現在の Hessian を使った 2 次の Taylor 展開で見積もります。続いて、DWI（距離加重補間）面の上で、改良型の Bulirsch–Stoer 法による修正子をかけます。分岐は、TS の近くを出た後に RMS 勾配が 1 × 10⁻³ hartree/bohr を下回ったとき、エネルギーが上がったとき、1 ステップのエネルギー変化が 1 × 10⁻⁶ hartree 以下になったとき、または `--max-cycles` に達したときに止まります。
3. **経路の書き出し**: 各分岐、TS を通る経路全体、その経路の両端の構造を書き出します。PDB/mmCIF の入力では、軌跡を PDB にも変換します。

---

## IRC の成否の判定

IRC が収束しなくても、端点を `opt` で最適化して狙った R と P に着けば、その結果は使えます。

| 確かめること | 見る場所 |
| --- | --- |
| 出発点が TS か | 端末の `Transition vector is mode 0 with wavenumber … cm⁻¹.` の行の波数が負 |
| 各分岐の止まり方 | `result.json` の `forward_integration_converged` / `backward_integration_converged`。RMS 勾配が閾値を下回ったときは `true`、エネルギーで止まったときやサイクルの上限では `false` |
| 経路に沿って変わる結合 | `result.json` の `bond_changes`（`finished_first` から `finished_last` への `formed` と `broken`） |
| どちらの端が R でどちらが P か | `finished_first.xyz` と `finished_last.xyz` を [`opt`](opt.md) で最適化し、意図した R と P と比べる。first / last の順では決まらない |

`result.json` の `scientific_status` の `success`（終了コード 0）は、各分岐の止まり方によらず、積分がエラーなく終わったことを示します。`irc` は端点を判定しないので、端点が狙った R と P かは自分で確かめてください。

端点が意図した R と P でないときは、{ref}`TS が取れないとき <ja-ts-search-fails>` を参照してください。

---

## 主な出力ファイル

実行が終わると、`--out-dir` に次のファイルができます。

```text
result_irc/
├─ finished_irc_trj.xyz    # IRC 経路全体：順方向の端 → TS → 逆方向の端
├─ finished_irc.pdb        # 同じ経路の PDB（PDB/mmCIF 入力）
├─ finished_first.xyz      # 最初のフレーム：順方向の端（--no-forward では TS）
├─ finished_last.xyz       # 最後のフレーム：逆方向の端（--no-backward では TS）
├─ {forward,backward}_{first,last}.xyz  # 各分岐の両端
├─ forward_irc_trj.xyz     # 順方向の分岐（実行したとき）
├─ forward_irc.pdb         # 同じ分岐の PDB（PDB/mmCIF 入力）
├─ backward_irc_trj.xyz    # 逆方向の分岐（実行したとき）
├─ backward_irc.pdb        # 同じ分岐の PDB（PDB/mmCIF 入力）
└─ result.json             # 結果の要約（--out-json）
```

{ref}`mmCIF の入力 <ja-mmcif-input>`と、PDB の欄に入りきらない大きな PDB の入力では、元の識別子を保った `.cif` も書きます。

* **端点の候補**: `finished_first.xyz` と `finished_last.xyz` を [`opt`](opt.md) で最適化します。各分岐の両端のうち、`forward_first.xyz` と `backward_last.xyz` は TS から遠い端、`forward_last.xyz` と `backward_first.xyz` は TS からの最初のステップです。
* **経路**: `finished_irc_trj.xyz` か `finished_irc.pdb` を PyMOL や VMD で開くと、反応の動きを見られます。
* **要約**: `--out-json` を付けると、[`result.json`](json-output.md) に各分岐のフレーム数、各分岐の止まり方、`bond_changes`、両端と TS のエネルギーが記録されます。
* **ファイル名の接頭辞**: YAML で `irc.prefix: trial` とすると、軌跡と構造のファイルの名前が `trial_` で始まります。`result.json` の `files` にも接頭辞つきの名前が記録されます。
* **端末**: 各分岐のステップの表と実行時間が出ます。{ref}`-v 3 <ja-verbosity-levels>` では、実際に使った `geom`・`calc`・`irc` の設定も出ます。

---

## 主な CLI オプション

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス | （必須） | TS の構造（`.pdb`, `.cif`, `.mmcif`, `.xyz`, `.gjf`）。軌跡は 1 フレームを `.xyz` に切り出してから指定（{ref}`軌跡から 1 フレームを取り出す <ja-trajectory-one-frame>` を参照） |
| `-q, --charge` | 整数 | `None` | 系全体の総電荷。`-l` を使う場合と `.gjf` 入力のほかは必須 |
| `-l, --ligand-charge` | 文字列 | `None` | リガンドの総電荷（例: `-1`）または残基名ごとの電荷（例: `'GPP:-3,SAM:1'`）。`-q` を省いたときに使用（PDB/mmCIF 入力または `--ref-pdb`） |
| `-m, --multiplicity` | 整数 | `1` | スピン多重度（2S+1） |
| `-b, --backend` | 文字列 | `uma` | バックエンド（`uma`, `orb`, `mace`, `aimnet2`, `dft`） |
| `--max-cycles` | 整数 | `125` | 分岐ごとの IRC ステップの上限 |
| `--step-size` | 実数 | `0.10` | 最大ステップ長（bohr、質量加重しない Cartesian 座標） |
| `--never-stop/--no-never-stop` | フラグ | `False` | 勾配とエネルギーによる停止の条件を無視し、`--max-cycles` までたどる |
| `--forward/--no-forward` | フラグ | `True` | 順方向の分岐を実行 |
| `--backward/--no-backward` | フラグ | `True` | 逆方向の分岐を実行 |
| `--root` | 整数 | `0` | 反応モードとする Hessian の固有ベクトル。固有値の昇順に 0 から数える |
| `--hessian-calc-mode` | `FiniteDifference` / `Analytical` | `FiniteDifference` | 最初の Hessian の計算方法 |
| `--read-hess` | パス | `None` | Hessian を計算せず、`.npy` ファイル（`freq` や `tsopt --dump-hess` で書いたものなど）から読んで始める |
| `--freeze-links/--no-freeze-links` | フラグ | `True` | キャップ水素の親原子を固定（PDB/mmCIF 入力または `--ref-pdb`） |
| `--out-json/--no-out-json` | フラグ | `False` | 結果の要約を `result.json` に出力（[JSON 出力の一覧](json-output.md)） |
| `-o, --out-dir` | パス | `./result_irc/` | 出力先ディレクトリ |

全オプションは [自動生成のオプションの一覧（英語のみ）](../reference/commands/irc.md) を参照してください。

> **補足:** YAML（`--config`）の `irc` ブロックのキーは、停止の閾値も含めて YAML 設定の一覧の {ref}`irc <ja-irc-section>` にすべて載っています。

---

## 使用上の注意点

* **すぐ止まる分岐**: 分岐が 3 フレーム以下で終わると、端末に `[irc] IRC stopped after only a few frames in …` の警告が出ます。ステップが大きすぎると EulerPC が不安定になることがあるので、先に例 4 を試してください。
* **`--never-stop` でも止まる場合**: `--never-stop` を付けても、数値的な失敗や外部からの中断では止まります。軌跡を確かめて端点を最適化し、その先の経路が役に立つときだけ `--max-cycles` を増やしてください。
* **`--root` は 0 から数える**: TS 最適化が成功すると、反応モードの虚振動が 1 つ出るので、n_imag = 1 の TS では `--root 0`（ただ 1 つの負の固有値）のままにしてください。`1`、`2` などは、反応モードより固有値の小さい（より負の）疑似モードがあると分かっているときだけ使います。
* **変えられない設定**: YAML の `geom.coord_type` と `calc.return_partial_hessian` にかかわらず、`irc` は Cartesian 座標と、動ける原子だけの Hessian を使います。
* **`--read-hess` のファイル**: [`freq`](freq.md) と同じ `.npy` ファイルです。`irc.hessian_init: calc`（デフォルト）が必要です。ファイルを使ったときは、`result.json` の `rigid_projection.hessian_source` が `"file"` になります。
* **解析 Hessian と `--uma-workers`**: UMA では、`--hessian-calc-mode Analytical` は 1 より大きい `--uma-workers` と併用できず、エラーで止まります。解析 Hessian には `--uma-workers 1` を使ってください。速度とメモリ量はバックエンド・モデル・系の大きさによって変わるので、先に対象の系で試してください。
* **端点の最適化**: `finished_first.xyz` と `finished_last.xyz` は `.xyz` だけで書かれるので、TS の PDB を `--ref-pdb` で渡し、キャップ水素の親原子を{ref}`固定 <ja-freeze-atoms-and-restraints>`したまま最適化してください。

  ```bash
  pdb2reaction opt -i result_irc/finished_first.xyz --ref-pdb ts.pdb -q 0 -m 1 --out-dir ./result_opt_first
  pdb2reaction opt -i result_irc/finished_last.xyz --ref-pdb ts.pdb -q 0 -m 1 --out-dir ./result_opt_last
  ```

* **固定原子**: `--freeze-links` に加えて、`--freeze-atoms`（1 始まり）でほかの原子も固定できます。除いた剛体運動と最初の Hessian は、`result.json` の `rigid_projection` に記録されます。
* **大きな系**: `--hess-device cpu` を付けると、最初の Hessian と IRC の Hessian の演算を CPU で行い、GPU のメモリに収めます。
* **分岐は少なくとも 1 つ**: `--no-forward` と `--no-backward` を両方付けると、エラーで止まります。

---

## 関連ドキュメント

* [tsopt](tsopt.md) — IRC の前に TS を最適化する
* [opt](opt.md) — IRC の端点を R と P へ最適化する
* [freq](freq.md) — 振動解析と熱化学補正
* [all](all.md) — `tsopt` の後に IRC を実行し、端点まで最適化する一連のワークフロー
* [トラブルシューティング](troubleshooting.md) — 実行が失敗したときの切り分け
* [YAML 設定の一覧](yaml-reference.md) — `irc` のすべての設定
* [用語集](glossary.md) — IRC などの用語
* {ref}`終了コード <ja-exit-codes>` — 終了ステータスの意味
