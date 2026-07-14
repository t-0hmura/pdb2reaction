# 出力ディレクトリのレイアウト

このページでは、各 `pdb2reaction` サブコマンドが出力ディレクトリに書き出すファイルと、エージェントや下流スクリプトが依拠すべき規約をまとめます。

## ファイル名の規約

| ファイル名 | 書き出し元 | 用途 |
|---|---|---|
| `summary.json` | `all`、`path-search`（常時） | aggregate workflow の標準 JSON envelope（[JSON 出力リファレンス](json-output.md)）。最初にこれを読んでください。 |
| `summary.json` | 成功したstage runでは `--out-json` 指定時（default `--no-out-json`）。捕捉runtime errorはflagなしでもbest-effort error envelopeを書く場合あり | stage JSON envelope。純粋utility（`fix-altloc`、`add-elem-info`、`bond-summary`など）はこのstage envelopeを使いません。 |
| `result.json` | stage `summary.json` と同条件（`opt`、`tsopt`、`freq`、`irc`、`sp`、scan系、`path-opt`、`dft`、`extract`） | 別名file。`summary.json` と同一payloadです。 |
| `summary.log` | `path-search`、`all` | 人間可読な実行ログ（セグメント／ステージごとに 1 行）。 |
| `final_geometry.xyz` | `opt`、`tsopt` | 最適化された構造（XYZ、完全精度）。 |
| `mep.pdb` / `mep.cif` / `mep_trj.xyz` | `path-search` | 反応経路のフレーム。変換が有効な mmCIF／oversized-PDB topology では `.cif` companion も追加。 |
| `final_geometries_trj.xyz` / `hei.xyz` | `path-opt` | スタンドアロンの path-opt 軌跡と最高エネルギーイメージ（変換が有効な場合は `.pdb` / `.cif` / `.gjf` companion も生成）。 |
| `mep_plot.png` | `path-search` | MEP のエネルギープロファイル（PNG）。（`all` では代わりに整形済みの `energy_diagram_MEP.png` をルートに配置します。） |
| `finished_irc_trj.xyz` / `forward_irc_trj.xyz` / `backward_irc_trj.xyz` | `irc` | IRC 軌跡（参照 topology があれば `.pdb`、bridge topology では `.cif` companion も生成）。 |
| `frequencies_cm-1.txt` | `freq` | 振動モードの一覧。 |
| `*.cif` / `*.gjf` | 各種（`--convert-files` 指定時） | 入力 template に応じた、元IDを保持する mmCIF または Gaussian companion。 |

## デフォルトの `--out-dir`

| サブコマンド | デフォルトの `--out-dir` |
|---|---|
| `all` | `./result_all/` |
| `opt` | `./result_opt/` |
| `tsopt` | `./result_tsopt/` |
| `freq` | `./result_freq/` |
| `irc` | `./result_irc/` |
| `dft` | `./result_dft/` |
| `scan` / `scan2d` / `scan3d` | `./result_scan*/` |
| `path-opt` / `path-search` | `./result_path_*/` |
| `sp` | `./result_sp/` |
| `extract` | `./`（`model.pdb` を書き出し。入力が複数の場合は `model_<input>.pdb`） |

`--out-dir <path>`（または `-o`）で上書きします。

## スタンドアロン と `all` の比較

単独で実行したサブコマンドは**フラット**な結果ディレクトリを書き出します。同じ書き出し処理でも、`all` によってオーケストレーションされると構造化されたツリーにネストされます。この 2 つのレイアウトは設計上異なります。

- **スタンドアロンのサブコマンド** → 上記のファイルを含むフラットな `result_<subcmd>/`。`segments/` も `_work/` もありません。これらは `all` が 1 回の実行で複数の書き出し処理を協調させるときにのみ現れます。
- **`all` の内部では、リーフの書き出し処理はそのままネストされます。** `segments/seg_NN/<subcmd>/` にあるセグメント別のリーフ出力は、スタンドアロンの `result_<subcmd>/` と構造的に同一です — `all` は同じ書き出し処理に別の出力ディレクトリを渡すだけです。
- **`path-search` / `path-opt` はエンジンの例外です。** スタンドアロンで実行すると、それぞれの出力が成果物となります: `path-search` → `result_path_search/`（`summary.log`、`mep.pdb`、bridge入力時の`mep.cif`、`mep_trj.xyz`、`mep_plot.png`、`energy_diagram_MEP.png`）、`path-opt` → `result_path_opt/`（`final_geometries_trj.xyz`、`hei.xyz`）。`all` の内部では、その生のエンジン出力は `_work/path_opt/`（`--refine-path` 指定時は `_work/path_search/`）下のスクラッチとして扱われ、マージされた成果物（`mep.pdb`、bridge入力時の`mep.cif`、`mep_trj.xyz`、`mep_w_ref.pdb` / `.cif`、`energy_diagram_MEP.png`）のみがパイプラインのルートに配置されます。
したがって `all` のツリーには 3 つのゾーンがあります。

```text
result_all/
├─ summary.log · summary.json                 # ルートに書き出し
├─ mep.{pdb,cif} · mep_w_ref.{pdb,cif} · mep_trj.xyz # CIF は bridge 入力時
├─ energy_diagram_MEP.png · energy_diagram_*.png
├─ segments/
│  └─ seg_NN/                                  # 反応セグメント別の成果物（2桁番号）
│     ├─ reactant.{pdb,cif,xyz,gjf} · ts.* · product.* # 正準の R/TS/P
│     └─ ts/ · irc/ · freq/{R,TS,P}/ · dft/         # ステージ別の作業ファイル（--tsopt / --thermo / --dft）
└─ _work/                                      # パイプラインのスクラッチ（rm -rf 可）
   ├─ models/ · scan/ · add_elem_info/ · fix_altloc/
   └─ path_opt/                                # MEP エンジンの生出力（--refine-path 時は path_search/）
```

TSOPT のみのモードでは MEP ステージがないため、`_work/path_opt/` は存在せず、成果物は `segments/seg_01/` 下に置かれます。モードごとの完全な内訳は [all](all.md) を参照してください。

## エージェント向けレシピ

```python
# Read whichever subcommand's output, single filename across the board.
import json
from pathlib import Path

summary = json.loads((Path(out_dir) / "summary.json").read_text())

if summary["status"] == "error":
    chain = summary.get("error_class_chain", [])
    if "OptimizationError" in chain:
        # retry with looser convergence threshold
        ...
    else:
        raise RuntimeError(summary["error"])
```

成功したstage runで `summary.json` / `result.json` を書くのは `--out-json`
指定時だけです。捕捉runtime errorはflagなしでもbest-effort error envelopeを書く場合が
ありますが、validation exitやoutput setup前の失敗では何も書かれない場合があります。
書かれたenvelopeはschema version + status（error pathではclass chain）を保持します。
