# `path-opt` サブコマンド

## 概要
`pdb2reaction path-opt` は、`--mep-mode` で選択されるpysisyphusのGrowing String法（GSM）またはDirect Max Flux（DMF）を使用して、2つのエンドポイント構造間の最小エネルギー経路（MEP）を探索します。UMAはすべてのイメージにエネルギー/勾配/ヘシアンを提供します。GSMがデフォルトの経路生成器であり、単一構造最適化は `light`（LBFGS）プリセットがデフォルトです。

## 使用法
```bash
pdb2reaction path-opt -i REACTANT.{pdb|xyz} PRODUCT.{pdb|xyz} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [-m MULT] \
                      [--workers N] [--workers-per-node N] \
                      [--mep-mode {gsm|dmf}] [--freeze-links {True\|False}] [--max-nodes N] [--max-cycles N] \
                      [--climb {True\|False}] [--dump {True\|False}] [--thresh PRESET] \
                      [--preopt {True\|False}] [--preopt-max-cycles N] [--opt-mode light|heavy] [--fix-ends {True\|False}] \
                      [--out-dir DIR] [--args-yaml FILE] \
                      [--convert-files {True\|False}] [--ref-pdb FILE]
```

## ワークフロー
1. **事前アライメント & 凍結解決**
   - 最初以降のすべてのエンドポイントは最初の構造にKabschアライメントされる
   - `--freeze-links=True`（デフォルト）のPDB入力では、リンク水素の親原子が検出され `freeze_atoms` にマージ

2. **ストリング成長とHEIエクスポート**
   - 経路が成長・精密化された後、最高エネルギー内部局所極大を探索
   - 最高エネルギーイメージ（HEI）は `.xyz` および `.pdb`（PDB参照が存在する場合）として書き込み

## CLIオプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH PATH` | 反応物と生成物構造 | 必須 |
| `-q, --charge INT` | 総電荷 | テンプレート/`0` |
| `--ligand-charge TEXT` | `-q` が省略された場合に使用される総電荷または残基名ごとのマッピング | _None_ |
| `-m, --multiplicity INT` | スピン多重度 | テンプレート/`1` |
| `--freeze-links {True\|False}` | PDBのみ: リンクH親を凍結 | `True` |
| `--max-nodes INT` | 内部ノード数（ストリングイメージ = `max_nodes + 2`） | `10` |
| `--mep-mode {gsm\|dmf}` | GSM（ストリングベース）またはDMF（ダイレクトフラックス）経路生成器を選択 | `gsm` |
| `--max-cycles INT` | オプティマイザーマクロイテレーション上限 | `300` |
| `--climb {True\|False}` | クライミングイメージ精密化を有効化 | `True` |
| `--dump {True\|False}` | MEP軌跡/リスタートをダンプ | `False` |
| `--opt-mode TEXT` | エンドポイント事前最適化用の単一構造オプティマイザー | `light` |
| `--convert-files {True\|False}` | PDB/Gaussian入力用のXYZ/TRJ → PDB/GJFコンパニオンをトグル | `True` |
| `--ref-pdb FILE` | XYZ/GJF入力用の参照PDBトポロジー | _None_ |
| `--out-dir TEXT` | 出力ディレクトリ | `./result_path_opt/` |
| `--thresh TEXT` | GSM/ストリングオプティマイザーの収束プリセットをオーバーライド | `gau` |
| `--args-yaml FILE` | YAMLオーバーライド（セクション `geom`、`calc`、`gs`、`opt`） | _None_ |
| `--preopt {True\|False}` | アライメント/MEP探索前に各エンドポイントを事前最適化 | `False` |
| `--fix-ends {True\|False}` | GSM成長/精密化中にエンドポイント構造を固定 | `False` |

## 出力
```
out_dir/
├─ final_geometries.trj        # XYZ経路
├─ final_geometries.pdb        # PDB参照が利用可能で変換が有効な場合
├─ hei.xyz                     # 最高エネルギーイメージ
├─ hei.pdb                     # PDB参照が利用可能な場合のHEI（変換有効時）
├─ hei.gjf                     # Gaussianテンプレートを使用して書き込まれたHEI（変換有効時）
└─ <オプティマイザーダンプ/リスタート>
```
