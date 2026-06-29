# `path-opt`

`pdb2reaction path-opt` は **ちょうど 2 つの構造**間の最小エネルギー経路（MEP）を、GSM（デフォルト）または DMF（`--mep-mode dmf`）で探索します。経路軌跡を書き出し、最高エネルギー画像（HEI）を TS 候補として出力します。HEI は *候補* に過ぎないため、[tsopt](tsopt.md)（虚振動数チェックを内蔵）→ [irc](irc.md) による接続性の確認が必須です。**2 つ以上の構造**を入力して反応領域だけを自動で精密化したい場合は、[path-search](path-search.md) を使用してください。

反応物と生成物の **2 端点**（R → P）が揃っており、再帰的な精密化なしで MEP の初期推定だけが必要な場面で使用します。ストリングベースの経路生成には GSM（デフォルト）を、Direct Max Flux 生成器には `--mep-mode dmf` で DMF を選択します。

MLIP バックエンド（デフォルト: UMA、`-b/--backend` で ORB ・ MACE ・ AIMNet2 も選択可能）で各イメージのエネルギー/勾配/ヘシアンを評価します。最適化の前に剛体アライメントを行い、ストリングの安定性を向上させます。`freeze_atoms` を指定した場合、RMSD フィットにはその原子群のみを使用しますが、変換自体は全原子に適用されます。

```{note}
**DMF モードでの凍結原子**は、GSM で使用される pysisyphus のハード座標凍結ではなく、`HarmonicFixAtoms`（k=300 eV/Å² の調和拘束）を使用します。そのため、DMF での凍結原子は参照位置からわずかに移動する可能性があり、GSM モードの剛体凍結とは挙動が異なります。
```

## 実行例

コマンド形式:

```bash
pdb2reaction path-opt -i REACTANT.{pdb|xyz} PRODUCT.{pdb|xyz} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m MULT] \
 [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] \
 [--workers N] [--workers-per-node N] \
 [--mep-mode {gsm|dmf}] [--freeze-links/--no-freeze-links] [--max-nodes N] [--max-cycles N] \
 [--climb/--no-climb] [--dump/--no-dump] [--thresh PRESET] [--thresh-stopt PRESET] \
 [--preopt/--no-preopt] [--preopt-max-cycles N] [--opt-mode grad|hess] [--fix-ends/--no-fix-ends] \
 [--show-config/--no-show-config] [--dry-run/--no-dry-run] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE]
```

2 端点間の MEP 探索:

```bash
pdb2reaction path-opt -i reactant.pdb product.pdb -q 0 -m 1 \
 --out-dir ./result_path_opt
```

MEP 探索前に端点を事前最適化する:

```bash
# MEP 探索前に端点を事前最適化する
pdb2reaction path-opt -i reactant.pdb product.pdb -q 0 -m 1 \
 --preopt --preopt-max-cycles 20000 --out-dir ./result_path_opt_preopt
```

GSM ではなく DMF モードで実行する:

```bash
# GSM ではなく DMF モードで実行する
pdb2reaction path-opt -i reactant.pdb product.pdb -q 0 -m 1 \
 --mep-mode dmf --max-nodes 12 --out-dir ./result_path_opt_dmf
```

```{note}
DMF モードは追加で `cyipopt` が必要です（`--mep-mode dmf` 実行前に conda-forge からインストールしてください）。`pydmf` は `pdb2reaction` の依存として同梱されています。デフォルトの `--dmf-backend gpu` は PyTorch/CUDA の `dmf.torch` バックエンドを使用します。GPU メモリ不足時は `--dmf-backend cpu`（`dmf`/NumPy）を指定してください。
```

キャップ親原子を凍結し、climb を切って短時間で確認するには `--freeze-links --no-climb` を追加します。

## 処理の流れ

1. **事前アライメント & 凍結解決**
 - 2 番目以降のエンドポイントは最初の構造に対して Kabsch アライメントされます。いずれかのエンドポイントで `freeze_atoms` が定義されている場合、RMSD フィットにはその原子のみを使用しますが、得られた変換は全原子に適用されます。
 - `--freeze-links` が有効な場合、キャップ水素の親原子は自動的に凍結されます（{ref}`キャップ水素と凍結原子 <ja-link-hydrogen-and-frozen-atoms>` を参照）。

2. **ストリング成長と HEI エクスポート**
 - 経路の成長・精密化後、内部ノード間の局所極大のうちエネルギーが最も高いものを優先的に選択します。内部の局所極大がない場合は内部ノードの最大値に、内部ノードもない場合は全体の最大値にフォールバックします。
 - 最高エネルギー画像（HEI）は `.xyz` として書き込まれます。PDB 参照がある場合は `.pdb`、Gaussian テンプレートがある場合は `.gjf` も出力します（いずれも `--convert-files` の設定に従います）。

## 出力

```text
out_dir/
├─ final_geometries_trj.xyz # XYZ経路（コメント行にエネルギーを保持）
├─ final_geometries.pdb # PDB 参照が利用可能で変換が有効な場合の全画像 PDB
├─ final_geometries.gjf # Gaussian テンプレート検出時の対応する Gaussian（変換有効時）
├─ hei.xyz # 最高エネルギー画像（コメント行にエネルギーを保持）
├─ hei.pdb # PDB 参照が利用可能な場合のHEI（変換有効時）
├─ hei.gjf # Gaussian テンプレートを使用して書き込まれたHEI（変換有効時）
├─ align_refine/ # 剛体アライメント/リファイン段階の中間ファイル（アライメント実行時）
└─ <オプティマイザーダンプ> # `--dump` 指定時の軌跡ダンプ（リスタート YAML は YAML の `dump_restart` 経由のみ）
```

主要な出力ファイル:

- `result_path_opt/final_geometries_trj.xyz`
- `result_path_opt/hei.xyz`
- `result_path_opt/hei.pdb`（PDB 変換が有効な場合）

コンソールには解決済み YAML ブロックが出力され、GSM/DMF の MEP 進行状況とタイミングが報告されます。

設定の優先順位は {ref}`CLI 規約: 設定の優先順位 <ja-configuration-precedence>` を参照してください。

## CLI オプション

完全なフラグ一覧は生成された [コマンドリファレンス](../reference/commands/index.md) を参照してください。以下の表は説明が必要なオプションのみを扱い、網羅的な一覧を手で複製しません。

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH PATH` | 反応物と生成物構造（`.pdb`/`.xyz`）。入力は 2 構造のみ。形式は `geom_loader` に準拠 | 必須 |
| `-q, --charge INT` | 総電荷（`calc.charge`）。`.gjf` 以外では `--ligand-charge` 導出が成功しない限り必須（PDB 入力または `--ref-pdb` 付き XYZ/GJF）。`.gjf` テンプレートがあればそれを使用し、電荷メタデータが無い `.gjf` 入力は `-q` が無いと中断。両方指定時は `-q` が優先。解決順序は {ref}`CLI 規約: 電荷の指定 <ja-charge-specification>` を参照 | テンプレート/導出がない限り必須 |
| `-l, --ligand-charge TEXT` | 総電荷または残基別マッピング（`-q` 省略時）。PDB 入力（または `--ref-pdb` 付き XYZ/GJF）で extract と同じ全系電荷導出を起動します | _None_ |
| `--workers`, `--workers-per-node` | MLIP 予測器の並列度（workers > 1 で解析ヘシアン無効; UMA バックエンドのみ; `workers_per_node` は並列予測器に渡されます）。診断上の注意は {ref}`ja-workers-fd-downgrade` を参照 | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度（`calc.spin`） | テンプレート/`1` |
| `--freeze-links/--no-freeze-links` | PDB 入力（または `--ref-pdb` 付き XYZ/GJF）: キャップ H 親を凍結（YAML とマージ）。詳細は [extract](extract.md) を参照 | `True` |
| `--freeze-atoms TEXT` | 凍結する原子の 1 始まりインデックスをカンマ区切りで明示的に指定（例: `'1,3,5'`）。`--freeze-links` と併用可、任意の入力形式に適用 | _None_ |
| `--max-nodes INT` | 内部ノード数。**GSM:** 総イメージ数 = `max_nodes + 2`（端点 2 つは固定）。**DMF:** チェーン上の*移動可能な*イメージ数（端点の暗黙的展開なし） | `20` |
| `--mep-mode {gsm\|dmf}` | GSM（ストリングベース）または DMF（Direct Max Flux）経路生成器を選択 | `gsm` |
| `--dmf-backend {cpu\|gpu}` | DMF 計算バックエンド（`--mep-mode dmf` 時のみ）: `gpu`（`dmf.torch`/CUDA）または `cpu`（`dmf`/NumPy）。GPU メモリ不足時は `cpu` で再実行 | `gpu` |
| `--max-cycles INT` | MEP 最適化サイクル上限（`stopt.max_cycles`、`stopt.stop_in_when_full`、`dmf.max_cycles` を同時設定） | `300` |
| `--climb/--no-climb` | クライミングイメージ精密化を有効化（Lanczos 接線も同時切替） | `True` |
| `--dump/--no-dump` | MEP 軌跡をダンプ（GSM/DMF）。リスタート YAML は YAML で有効化した場合のみ書き出されます | `False` |
| `--opt-mode TEXT` | エンドポイント事前最適化用の単一構造オプティマイザー（`grad` = L-BFGS、`hess` = RFO） | `grad` |
| `--convert-files/--no-convert-files` | PDB/Gaussian 入力用の XYZ/TRJ → 対応する PDB/GJF 出力の切り替え | `True` |
| `--ref-pdb FILE` | XYZ/GJF 入力用の参照 PDB トポロジー（XYZ 座標は保持し PDB 変換を有効化） | _None_ |
| `-o, --out-dir TEXT` | 出力ディレクトリ | `./result_path_opt/` |
| `--thresh TEXT` | エンドポイント事前最適化のみの収束プリセットを上書き（`opt.lbfgs/rfo.thresh`） | `gau` |
| `--thresh-stopt TEXT` | ストリングオプティマイザー（GSM 成長およびクライミング精密化）の収束プリセットを上書き（`stopt.thresh`; `gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） | `gau_loose` |
| `--config FILE` | 明示 CLI 指定より前に適用されるベース YAML | _None_ |
| `--show-config/--no-show-config` | 解決済み設定（YAML レイヤ情報を含む）を表示して実行継続 | `False` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP バックエンド | `uma` |
| `--solvent TEXT` | xTB 暗黙溶媒（例: `water`）。`none` で無効化 | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB 溶媒モデル | `alpb` |
| `--dry-run/--no-dry-run` | 実行せずに検証と実行計画表示のみを行う | `False` |
| `--preopt/--no-preopt` | アライメント/MEP 探索前に各エンドポイントを事前最適化（GSM/DMF）。 | `True` |
| `--preopt-max-cycles INT` | エンドポイント事前最適化サイクルの上限 | `10000` |
| `--fix-ends/--no-fix-ends` | GSM 成長/精密化中にエンドポイント構造を固定 | `True` |
| `--out-json/--no-out-json` | `out_dir` に機械可読な `result.json` を書き出す。スキーマは [JSON 出力スキーマ](json-output.md) を参照 | `False` |

## YAML 設定

### `path-opt` で使用される YAML セクション

完全なキー一覧は [YAML リファレンス](yaml-reference.md) を参照:

- [`geom`](yaml-reference.md#geom) — PDB 入力では `--freeze-links` が `freeze_atoms` にマージされます。
- [`calc`](yaml-reference.md#calc) — MLIP バックエンド設定。
- [`gs`](yaml-reference.md#gs) — Growing String 表現（GSM モード）。
- [`dmf`](yaml-reference.md#dmf) — Direct Max Flux + (C)FB-ENM 補間（DMF モード）。
- [`stopt`](yaml-reference.md#stopt) — StringOptimizer 設定。
- [`opt.lbfgs`](yaml-reference.md#lbfgs) / [`opt.rfo`](yaml-reference.md#rfo) — エンドポイント単一構造事前最適化。YAML が CLI の `--preopt-max-cycles` を上書きします。

### `path-opt` 固有のデフォルト

`path-opt` 経由で実行した場合、以下のキーが正規デフォルトと異なります:

```yaml
stopt:
 out_dir: ./result_path_opt/ # output directory (path-opt default)
opt:
 lbfgs:
   out_dir: ./result_path_opt/ # output directory (path-opt default)
 rfo:
   out_dir: ./result_path_opt/ # output directory (path-opt default)
```

## 終了コード

CLI 規約の {ref}`ja-exit-codes` を参照してください。

## 関連項目

- [path-search](path-search.md) — 自動精密化を伴う再帰的 MEP 探索（2+構造用）
- [tsopt](tsopt.md) — HEI を TS 候補として最適化（内部で虚振動数チェック済み）。続けて IRC で接続性を確認
- [extract](extract.md) — path-opt 入力用の活性部位モデル PDB を生成
- [all](all.md) — 一気通貫ワークフロー（デフォルトで単一パス path-opt を使用; `--refine-path True` で再帰的 path-search に切替。`--refine-path` フラグは `pdb2reaction all` にのみ属します — 定義は {ref}`ja-mep-search-options` を参照してください）
- [YAML リファレンス](yaml-reference.md) — `gs`、`dmf`、`stopt`、`opt` の完全な設定オプション
- [用語集](glossary.md) — MEP、GSM、DMF、HEI の定義
- [典型エラー別レシピ](recipes-common-errors.md) — 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) — 詳細な対処ガイド
