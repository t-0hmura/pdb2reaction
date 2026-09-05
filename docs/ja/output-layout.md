# 出力ディレクトリのレイアウト

各 `pdb2reaction` サブコマンドの出力ファイルと配置をまとめます。

## ファイル名の規約

| ファイル名 | 書き出し元 | 用途 |
|---|---|---|
| `summary.json` | 集約結果の書き込み処理まで到達した `all` / `path-search` | 集約ワークフローの正規 JSON エンベロープ（[JSON 出力リファレンス](json-output.md)）。早期の CLI 引数または入力の検証では作られない場合があります。 |
| `summary.json` | `--out-json`（デフォルト: `--no-out-json`）を指定し、個別計算・レポートの結果を書き出した場合。DFT の非収束時にも出力。早期失敗時はエラー情報だけの場合あり | `result.json` の互換用コピー。書き込み正常終了時は同一内容。`fix-altloc`、`add-elem-info`、`bond-summary` は対象外。 |
| `result.json` | 個別計算の `summary.json` と同じ条件 | 個別計算・レポートの正規 JSON。非収束時にも生成される場合があるが、早期の入力検証・インポートエラーでは存在しない場合あり。 |
| `run.log` | コマンドへ到達し、出力ディレクトリが作られた CLI / Colab 実行 | shell-safe な実行コマンドと、コマンド実行中の標準出力・標準エラー。早期の Click 検証、help、version、dry-run、出力先が単一ファイルのユーティリティでは生成しません。 |
| `summary.log` | `path-search`、`all` | 実行要約（セグメント／ステージごとに 1 行）。 |
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

通常は `--out-dir <path>`（または `-o`）で上書きします。`extract` だけは repeatable な `-o/--output <file>` を使います。

## スタンドアロン と `all` の比較

単独実行では `result_<subcmd>/` にファイルが並び、`segments/` や `_work/` はありません。`all` では、後処理の各段階が同じファイル構成で `segments/seg_NN/<subcmd>/` に配置されます。

- **`path-search` / `path-opt` はエンジンの例外です。** スタンドアロンで実行すると、それぞれの出力が成果物となります: `path-search` → `result_path_search/`（`summary.log`、`mep.pdb`、bridge入力時の`mep.cif`、`mep_trj.xyz`、`mep_plot.png`、`energy_diagram_MEP.png`）、`path-opt` → `result_path_opt/`（`final_geometries_trj.xyz`、`hei.xyz`）。`all` の内部では、その生のエンジン出力は `_work/path_opt/`（`--refine-path` 指定時は `_work/path_search/`）下のスクラッチとして扱われ、主要成果物（`mep.pdb`、bridge入力時の`mep.cif`、`mep_trj.xyz`、`--write-ref-merge` 指定時の確認用`mep_w_ref.pdb` / `.cif`、`energy_diagram_MEP.png`）のみがパイプラインのルートに配置されます。
したがって `all` のツリーには 3 つのゾーンがあります。

```text
result_all/
├─ summary.log · summary.json                 # ルートに書き出し
├─ mep.{pdb,cif} · mep_trj.xyz                       # MEP座標
├─ mep_w_ref.{pdb,cif}                               # 確認用座標composite（--write-ref-merge）
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
# 実行したコマンドに対応する正規のファイル名を選ぶ。
import json
from pathlib import Path

subcommand = "opt"  # 実行したコマンドに置き換える
out_dir = "result_opt"  # 実際の出力ディレクトリに置き換える
primary = "summary.json" if subcommand in {"all", "path-search"} else "result.json"
summary = json.loads((Path(out_dir) / primary).read_text())

if summary["status"] == "error":
    error_type = summary.get("error_type", "RuntimeError")
    raise RuntimeError(f"{error_type}: {summary['error']}")
```

正常終了した段階別コマンドが `summary.json` / `result.json` を書くのは、
`--out-json` 指定時だけです。捕捉した実行時エラーでは、フラグなしでも可能な範囲で
エラーエンベロープを書く場合がありますが、入力検証による終了や出力先の確定前に
失敗した場合は何も書かれないことがあります。書き出されたエンベロープは、
スキーマバージョンとステータス（エラー時はクラス階層を含む）を保持します。
