# 出力ディレクトリのレイアウト

このページでは、各 `pdb2reaction` サブコマンドが出力ディレクトリに書き出すファイルと、エージェントや下流スクリプトが依拠すべき規約をまとめます。

## ファイル名の規約

| ファイル名 | 書き出し元 | 用途 |
|---|---|---|
| `summary.json` | `all`、`path-search`、および `write_result_json` を実行するすべてのステージ別サブコマンド | 正式な JSON エンベロープ（[JSON 出力リファレンス](json-output.md) を参照）。最初にこれを読んでください。`write_result_json` を呼ばない純粋なユーティリティ系サブコマンド（例: `fix-altloc`、`add-elem-info`、`bond-summary`）はこれを出力しません。 |
| `result.json` | `write_result_json` を呼ぶステージ別サブコマンド（`opt`、`tsopt`、`freq`、`irc`、`sp`、`scan` / `scan2d` / `scan3d`、`path-opt`、`dft`、`extract`） | 別名ファイル — `summary.json` と内容（ペイロード）が同一です。単一ファイル名の規約に従う場合は `summary.json` を読んでください。`result.json` は同じ内容を持ち、`summary.json` のみを利用するなら削除して構いません。 |
| `summary.log` | `path-search`、`all` | 人間可読な実行ログ（セグメント／ステージごとに 1 行）。 |
| `final_geometry.xyz` | `opt`、`tsopt` | 最適化された構造（XYZ、完全精度）。 |
| `mep.pdb` / `mep_trj.xyz` | `path-search` | 反応経路のフレーム（PDB / XYZ）。 |
| `final_geometries_trj.xyz` / `hei.xyz` | `path-opt` | スタンドアロンの path-opt 軌跡（経路全体）と最高エネルギーイメージ（変換が有効な場合は `.pdb` / `.gjf` コンパニオンも）。 |
| `mep_plot.png` | `path-search` | MEP のエネルギープロファイル（PNG）。（`all` では代わりにスタイル付きの `energy_diagram_MEP.png` をルートに昇格します。） |
| `finished_irc_trj.xyz` / `forward_irc_trj.xyz` / `backward_irc_trj.xyz` | `irc` | IRC 軌跡（経路全体およびブランチ別。参照 PDB が利用可能な場合は `.pdb` コンパニオンも）。 |
| `frequencies_cm-1.txt` | `freq` | 振動モードの一覧。 |
| `*.gjf` | 各種（`--convert-files` 指定時） | Gaussian 形式のコンパニオン構造。 |

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
- **`path-search` / `path-opt` はエンジンの例外です。** スタンドアロンで実行すると、`path-search` 自体が成果物です（独自の `summary.log`、`mep.pdb`、`mep_trj.xyz`、`mep_plot.png`、`energy_diagram_MEP.png` を持つ `result_path_search/`）。`all` の内部では、その生の出力は `_work/path_search/` 下のスクラッチとして扱われ、マージされた成果物（`mep.pdb`、`mep_trj.xyz`、`mep_w_ref.pdb`、`energy_diagram_MEP.png`）のみがパイプラインのルートに昇格されます。
したがって `all` のツリーには 3 つのゾーンがあります。

```text
result_all/
├─ summary.log · summary.json                 # authored at the root
├─ mep.pdb · mep_w_ref.pdb · mep_trj.xyz       # MEP deliverables promoted from the engine
├─ energy_diagram_MEP.png · energy_diagram_*.png
├─ segments/
│  └─ seg_NN/                                  # 2-digit per-reactive-segment deliverables
│     ├─ reactant.{pdb,xyz,gjf} · ts.* · product.*   # canonical R/TS/P
│     └─ ts/ · irc/ · freq/{R,TS,P}/ · dft/         # per-stage working files (--tsopt / --thermo / --dft)
└─ _work/                                      # pipeline scratch (safe to rm -rf)
   ├─ models/ · scan/ · add_elem_info/ · fix_altloc/
   └─ path_opt/                                # raw MEP-engine output (path_search/ with --refine-path True)
```

TSOPT のみのモードでは MEP ステージがないため、`_work/path_search/` は存在せず、成果物は `segments/seg_01/` 下に置かれます。モードごとの完全な内訳は [all](all.md) を参照してください。

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

`summary.json` は、サブコマンドが `write_result_json` を呼ぶ限り（失敗パスを含む — 失敗エンベロープはスキーマバージョン + エラークラスチェーンを保持します）必ず存在することが保証されます。その横に書き出される `result.json` は内容が同一です。
