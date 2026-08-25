# クイックスタート: スキャンを起点とする `pdb2reaction all`

## 目的

単一構造に対して `-s/--scan-lists` で 1 つ以上の距離座標を定義し、`pdb2reaction all` を実行します。拘束スキャンに続いて MEP を最適化し、必要に応じて遷移状態（TS）最適化と固有反応座標（IRC）解析へ進みます。

## 事前に必要なもの

- 入力構造: PDB/mmCIF、XYZ、または GJF
- 対象状態に対応した電荷（`-q/--charge` または `-l/--ligand-charge`）・多重度（`-m`）

## スキャンワークフローの選択

| 目的 | コマンド | 座標の扱い |
| --- | --- | --- |
| 拘束付き端点構造と軌跡のみを生成 | `pdb2reaction scan` | 複数距離の協奏scan・多段階scanに対応 |
| スキャン定義経路から MEP と任意の TS/IRC 解析へ継続 | `pdb2reaction all -s ...` | スキャン端点を `path-opt`、または `--refine-path` 時の `path-search` へ渡す |
| 2 次元または 3 次元の energy landscape を探索し PES を描画 | `pdb2reaction scan2d` / `scan3d` | 独立距離グリッドの直積 |

`scan` とスキャン起点の `all` では、1 リテラルが 1 ステージです。
同一リテラル内の複数タプルは協奏的に駆動し、複数リテラルは多段階scanとして順次実行されます。

## `-s/--scan-lists` インラインリテラル

`-s/--scan-lists` はコマンドライン上で Python リテラル文字列を直接受け取ります。原子セレクタの構文（残基/原子トークン、区切り文字、順序）と外側/内側のクォートのルールについては、{ref}`CLI 規約: スキャンリスト仕様 <ja-scan-list-spec>` を参照してください。

### 基本構文

各リテラルは 3 要素タプル `(atom1, atom2, target_distance_Å)` のリストです。3 番目の要素は必ず **ångström** 単位の目標距離で、ちょうど 3 要素が必要です。1 リテラル = 1 ステージ。

```bash
# 単一ステージ、整数原子インデックス（デフォルトで1-based）
pdb2reaction all -i input.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' -m 1 \
 -s '[(1, 5, 1.35)]' -o ./result_scan

# 単一ステージ、PDBセレクタ文字列
pdb2reaction all -i input.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' -m 1 \
 -s '[("TYR,285,CA", "SAM,309,C10", 1.35)]' -o ./result_scan
```

全系 PDB/mmCIF から活性部位モデルを抽出する場合は `-c/--center` を指定します。抽出済みクラスターモデル、または XYZ/GJF をそのまま解析する場合は省略します。`-m/--multiplicity` のデフォルトは `1`（一重項）ですが、ここでは明示しています。

### 複数ステージ

複数のリテラルを渡すと、各リテラルが順番に 1 ステージとして実行されます:

```bash
# ステージ1: 1つの結合を 1.35 Å に駆動
# ステージ2: 2つの結合を同時に駆動
pdb2reaction all -i input.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' -m 1 -s \
  '[("TYR,285,CA","SAM,309,C10",1.35)]' \
  '[("TYR,285,CA","SAM,309,C10",2.20),("TYR,285,CB","SAM,309,C11",1.80)]' \
  -o ./result_scan
```

ステージは順番に実行され、各々は前ステージの緩和結果から開始します。

## 期待される出力

```text
result_scan/
├── summary.log
├── summary.json
├── mep.pdb                        # MEP 全体パス（ルートへ配置）
├── energy_diagram_MEP.png         # MEP エネルギー図（ルートへ配置）
└── _work/                         # パイプラインの作業領域（削除可）
    ├── scan/
    │   ├── preopt/                # 事前最適化構造
    │   ├── stage_01/              # スキャンステージ 1
    │   │   ├── result.{xyz,pdb}   # 最終拘束端点（--scan-endopt 時のみ無拘束最適化後）
    │   │   ├── scan_trj.xyz       # スキャン軌跡
    │   │   └── scan.pdb
    │   └── stage_02/              # スキャンステージ 2（マルチステージ時）
    └── path_opt/                  # MEP 探索（--refine-path 時は path_search/）
        └── hei_seg_01.{xyz,pdb}   # MEP の最高エネルギー像
```

この最小コマンドは MEP ステージ終了時に停止するため、`segments/` は作成しません。
`--tsopt` を追加し、検証に成功するとセグメント別の正規 R/TS/P と IRC 出力が作成され、
`--thermo` も追加すると `freq/` が作成されます。

### 出力の検証

1. `_work/scan/stage_01/scan_trj.xyz` — 結合距離の変化を PyMOL で確認
2. `mep.pdb` と `_work/path_opt/hei_seg_01.pdb` — 最適化後の MEP と最高エネルギー像を確認
3. `summary.log` — 障壁高さと結合変化

入力、抽出、電荷、スキャン原子対応を含む `all` 全体を事前検証します:

```bash
pdb2reaction all -i input.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
  -s '[(1, 5, 1.35)]' --dry-run
```

standalone `scan --print-parsed` はスキャン仕様だけを検証し、`all` 固有の抽出や原子index remappingは検証しません。

## 補足

- `-s/--scan-lists` は `all` ではインライン Python リテラルのみを受け取ります。単独の `scan` サブコマンドはこれに加えて YAML/JSON スペックファイルパスも受け取れます（[scan](scan.md) を参照）。
- スキャンエンジンは共通ですが、parent option名とpreopt既定値が異なります:

  | コマンド | step / bias | 緩和上限 | preopt / endopt |
  | --- | --- | --- | --- |
  | `pdb2reaction all` | `--scan-max-step-size 0.20 Å`; `--scan-bias-k 300 eV/Å²` | `--scan-relax-max-cycles 100000` | scan preoptはparent `--preopt`（デフォルトON）を継承; endoptはOFF（`--no-scan-endopt`） |
  | `pdb2reaction scan` | `--max-step-size 0.20 Å`; `--bias-k 300 eV/Å²` | `--relax-max-cycles 100000` | `--no-preopt`; `--no-endopt` |

  step 幅は対応する CLI flag で上書きします。バイアス強度は CLI flag または YAML の `bias.k` で上書きできます（[scan](scan.md) / [yaml-reference](yaml-reference.md#bias) 参照）。
- 各スキャンステージは、最終緩和構造に対する結合変化チェック（`has_bond_change`）で終了します。ステージごとの結果は scan ログに記録され、`--out-json` 指定時には scan 出力ディレクトリに書き出される集約 `result.json`（その `stages` 配列内）にも記録されます。
- 再帰的 MEP refinement（`path-search`）は scan の端点を**無条件**に入力として取り込みます。実行されるかどうかのゲートは、scan stage の bond-change フラグ（`has_bond_change`）ではなく `--refine-path` です。
- scanの正常終了が示すのは拘束付き端点の生成です。無拘束極小、素反応、遷移状態の成立は、無拘束最適化、鞍点次数、IRCで別途検証します。
- `pdb2reaction all --help-advanced` で全オプション（スキャン制御を含む）を確認できます。
- 単独の `scan` サブコマンド（MEP 精密化なし）については [scan](scan.md) を参照してください。

## 次のステップ

- 全オプション: [all](all.md)
- TS 候補の検証: [クイックスタート: TS のみモード](quickstart-tsopt-freq.md)
