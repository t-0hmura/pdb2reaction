# MLIP の TS を DFT で確かめる

MLIP で妥当な経路が見つかったら、その TS をそのまま DFT での TS 構造最適化にもっていくことにも `pdb2reaction` は対応しています。TS 最適化 → IRC → 端点の最適化 → 振動数計算のワークフローを、GPU4PySCF を用いることで GPU で高速化された DFT 計算により実行可能です。

主役は MLIP による経路探索で、DFT は MLIP で見つけた TS 候補を確かめるための追加の機能です。

---

## 主な用途

- **TS を DFT で詰める**：TS 最適化 → IRC → 端点の最適化 → 振動数計算を DFT で実行します（`-b dft`）。
- **MLIP の計算に DFT のエネルギーを足す**：MLIP で得た R・TS・P に DFT の一点計算を行います（`--dft`）。

## 処理の流れ

1. **MLIP で探す**：経路を作り、条件を変えて試し、いちばん有望な TS 候補を選びます。
2. **DFT で詰める**：その TS を入力にして、`-b dft` 付きの TS-only モードを実行します。コマンドは[基本的な実行例](#基本的な実行例)の例 2 にあります。
3. **確かめる**：MLIP のときと同じく、TS 最適化が成功すると、反応モードの虚振動が 1 つ出ます。求めた段がすべて収束すると `====== Pipeline summary ======` の下に `Scientific status: success` と出ます。続けて、モードと IRC の端点を[結果の確認](quickstart-tsopt.md#結果の確認)のとおりに確かめてください。

## 基本的な実行例

### 1. 小さなモデルで MLIP の探索を行う

DFT で扱える大きさのモデルで MLIP の探索を行います。出てくる TS の `result_mlip/segments/seg_01/ts.pdb` を例 2 の入力にします。`seg_01` は最初の反応セグメントです。セグメントが複数あるときは、詰めたい段のものを選んでください。

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -r 0 --selected-resn '44,63,186' --tsopt -o ./result_mlip
```

`-r 0` で切り出しの半径を 0 Å にすると、距離で近くの残基を足すのを止め、`-c` と `--selected-resn` の残基からモデルを組みます。[`examples/`](https://github.com/t-0hmura/pdb2reaction/tree/main/examples) の同梱例の 44・63・186 番は、SAM のメチル炭素（CS1）にいちばん近い 3 残基です。同梱の PDB は chain の欄が空です。chain が空の PDB では、残基を名前か番号で指定してください。自分の系では、反応に関わる残基を選んでください。

### 2. TS を DFT で詰める

例 1 の TS を 1 つだけ入力にすると、TS-only モードになります。DFT 用の[追加パッケージ](#使用上の注意点)が要ります。

```bash
pdb2reaction all -i result_mlip/segments/seg_01/ts.pdb -l 'SAM:1,GPP:-3' \
    --tsopt --thermo -b dft -o ./result_dft
```

## モデルは 300 原子くらいまでに

DFT で最適化できるのは、多くても 300 原子くらいまでです。300 と比べるのはキャップ水素を含む原子数で、[モデルを確かめる](model-setup.md#モデルを確かめる)の端末の行から N + M として求めます。MLIP で探す前にこの大きさのモデルを作っておけば、出てきた TS をそのまま DFT に渡せます。

モデルの削り方は [モデルを削る](model-setup.md#モデルを削る)、手で削ったモデルの境界の原子の固定は {ref}`原子の固定と距離の拘束 <ja-freeze-atoms-and-restraints>` を見てください。

## `-b dft` と `--dft` の違い

| オプション | DFT で計算するもの | 使う場面 |
|---|---|---|
| `-b dft` | 計算のすべて（MEP 探索、TS 最適化、IRC、端点の最適化、振動数） | TS 候補を DFT で詰めて確かめる |
| `--dft` | MLIP で得た R・TS・P の一点計算（`all` だけ） | MLIP の構造で DFT のエネルギーを得る |

`-b dft` は 11 のコマンドで使えます：`all`、`opt`、`tsopt`、`irc`、`freq`、`scan`、`scan2d`、`scan3d`、`path-opt`、`path-search`、`sp`。

## 主な出力ファイル

`-b dft` の出力の並びは、同じモードの MLIP の計算と同じです（TS-only モードは[期待される出力](quickstart-tsopt.md#期待される出力)）。`--dft` を付けると、次のファイルが増えます。

| ファイル | 内容 |
|---|---|
| `segments/seg_NN/dft/{R,TS,P}/result.yaml` | 各状態の DFT 一点計算の結果 |
| `segments/seg_NN/energy_diagram_DFT.png` | MLIP の構造での DFT のエネルギー図 |
| `segments/seg_NN/energy_diagram_G_DFT_plus_MLIP.png` | DFT のエネルギーに MLIP の熱補正を足した図（`--thermo` のとき） |
| `energy_diagram_DFT_all.png`、`energy_diagram_G_DFT_plus_MLIP_all.png` | 全セグメントをまとめた同じ図（出力ディレクトリの直下） |

## 主な CLI オプション

| オプション | 説明 | 既定値 |
|---|---|---|
| `-b, --backend dft` | DFT を計算に使います（GPU4PySCF。`--dft-engine cpu` で CPU の PySCF）。 | `uma` |
| `--func-basis TEXT` | 汎関数と基底を `FUNCTIONAL/BASIS` の形で指定します。`-b dft` と `--dft` の両方に効きます。 | `wb97m-v/def2-svp` |
| `--solvent TEXT`、`--solvent-model [pcm\|smd]` | `-b dft` の連続溶媒です。溶媒名を指定すると SMD になり、`--solvent-model pcm` で PCM に変わります。 | `none`（気相） |
| `--dft/--no-dft` | R・TS・P に DFT の一点計算を足します（`all` だけ）。 | `--no-dft` |
| `--dft-solvent TEXT`、`--dft-solvent-model [pcm\|smd]` | `--dft` の一点計算だけの連続溶媒です。`--dft` なしで指定するとエラーで止まります。 | `none`、`smd` |

ほかの DFT のオプションは [`all` のオプションの一覧（英語のみ）](../reference/commands/all.md)にあります。

> **補足:** YAML では、`-b dft` の設定は `calc.dft` に書き、`--dft` の一点計算は `dft` コマンドと同じく最上位の `dft` セクションを読みます。どちらのブロックでも、`pyscf` で PySCF のオブジェクトに名前ごとに属性を渡せます。たとえば収束しにくい SCF には、`calc.dft.pyscf` か `dft.pyscf` の下に `mf: {level_shift: 0.2}` と書きます。キーの一覧は [YAML 設定の一覧](yaml-reference.md#calc)と {ref}`dft セクション <ja-dft-section>` にあります。

## 使用上の注意点

- **DFT 用の追加パッケージ**：PyTorch の wheel が `cu130`・`cu132` なら `pip install "pdb2reaction[dft]"`、`cu126` なら `pip install "pdb2reaction[dft-cuda12]"` で入れてください。GPU が無いときは `--dft-engine cpu` を付けてください。
- **電荷**：残基を除くと全体の電荷が変わります。DFT を流す前に、例 1 の端末の出力の `Total active site model charge` を確かめてください。
- **キャップ水素**：`ts.pdb` にはモデルのキャップ水素（`LKH`/`HL`）が残っているので、DFT の計算でも MLIP のときと同じく、その親原子が固定されます。
- **組み合わせ**：`-b dft` と `--dft` は一緒に使えず、実行の始めにエラーで止まります。`-b dft` の計算の後に DFT の一点計算を足すときは、別のジョブで `pdb2reaction sp -b dft` か `pdb2reaction dft` を実行してください。`--dft` と `--thermo` には `--tsopt` が必要です。
- **図のファイル名とキーの名前**：`-b dft` でも、図のファイル名は `energy_diagram_MLIP.png` と `energy_diagram_G_MLIP.png`（`--thermo` のとき）、`summary.json` のブロックの名前は `mlip` と `gibbs_mlip` のままです。中身は DFT の値で、図の題には DFT と出ます。
- **メモリとスレッド数**：`-b dft` と `--dft` は、[`pdb2reaction dft`](dft.md#主な-cli-オプション) と同じ低メモリモード・`--dft-nprocs`・`--dft-memory` を使います。GPU のメモリが足りないときはモデルを削り、`--dft` でメモリが足りないときは `--dft` を外して `pdb2reaction dft` を別に実行してください。
- **SCF が収束しないとき**：初期推測をやり直しても SCF が収束しないときは、`PySCF SCF did not converge with either the reused density or a fresh guess.` で計算が止まります。メモリに余裕があれば `--no-dft-low-memory` の密度フィッティングか、[YAML のレベルシフト](#主な-cli-オプション)で収束しやすくなることがあります。
- **SCF のチェックポイント**：ファイルがとても大きくなることがあるので、既定では保存しません。`--save-scf-checkpoint` で保存し、`--scf-checkpoint PATH` で使うファイルを選べます。パスを省くと、`all` 以外のコマンドは `<out-dir>/_work/dft_scf/state.chk` に、`all` は状態ごとに別のファイルに保存します。チェックポイントは、手法・原子の順番・座標が今の構造と一致するときだけ使われます。

## 関連ドキュメント

- [`dft`](dft.md)：DFT 一点計算とポピュレーション解析
- [`sp`](sp.md)：任意のバックエンドでの一点のエネルギーと力
- [クイックスタート：TS-only モード](quickstart-tsopt.md)：`all --tsopt` で TS 候補を確かめる
- [クラスターモデルの組み方](model-setup.md)：活性部位モデルを組む・削る・広げる
- [インストール](installation.md#詳細なインストール手順)：手順 7 で DFT 用の追加パッケージを入れる
- [MLIP バックエンド](backends.md)：バックエンドの選び方
- [トラブルシューティング](troubleshooting.md)：計算が失敗したとき
