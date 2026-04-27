# `path-opt`

## 概要

> **要約:** **2 つの構造**（R → P）から、GSM（デフォルト）または DMF（`--mep-mode dmf`）で MEP を最適化します。経路軌跡を書き出し、最高エネルギー画像（HEI）を TS 候補として出力します。

### 要点
- **想定場面:** 反応物と生成物の **2 端点**が揃っていて、まず MEP の初期推定を得たい場合に使います。
- **手法:** デフォルトは GSM。`--mep-mode dmf` で DMF に切り替え可能。
- **主な出力:** `final_geometries_trj.xyz`（経路）と `hei.xyz`（HEI）。変換が有効なら `.pdb`/`.gjf` コンパニオンも生成。
- **デフォルト値:** `--opt-mode grad`（LBFGS）、`--climb`、`--max-nodes 20`、`--thresh gau`、`--thresh-stopt gau_loose`。
- **次にやること:** HEI は **TS 候補**です。`tsopt`（内部で虚振動数チェック済み、**1 つ** であること）→ `irc` で検証します。

`pdb2reaction path-opt` は 2 端点間の最小エネルギー経路（MEP）を探索し、最高エネルギー画像（HEI）を報告します。HEI は *候補* に過ぎないため、[tsopt](tsopt.md)（内部で虚振動数チェック済み）→ [irc](irc.md) による接続性の確認が必須です。**2 つ以上の構造**を入力して反応領域だけを自動で精密化したい場合は、[path-search](path-search.md) を使用してください。

> **`path-opt` と `path-search` の使い分け:** 端点が 2 構造だけで再帰精密化が不要な場合は `path-opt` を使います。2 構造以上を入力し、結合変化のある領域を自動で再帰精密化したい場合は `path-search` を使います。

MLIP バックエンド（デフォルト: UMA、`-b/--backend` で ORB・MACE・AIMNet2 も選択可能）で各イメージのエネルギー/勾配/ヘシアンを評価します。最適化の前に剛体アライメントを行い、ストリングの安定性を向上させます。`freeze_atoms` を指定した場合、RMSD フィットにはその原子群のみを使用しますが、変換自体は全原子に適用されます。

```{note}
**DMF モードでの凍結原子**は、GSM で使用される pysisyphus のハード座標凍結ではなく、`HarmonicFixAtoms`（k=300 eV/Å² の調和拘束）を使用します。そのため、DMF での凍結原子は参照位置からわずかに移動する可能性があり、GSM モードの剛体凍結とは挙動が異なります。
```


## 最小例

```bash
pdb2reaction path-opt -i reactant.pdb product.pdb -q 0 -m 1 \
 --out-dir ./result_path_opt
```

## 出力の見方

- `result_path_opt/final_geometries_trj.xyz`
- `result_path_opt/hei.xyz`
- `result_path_opt/hei.pdb`（PDB 変換が有効な場合）

## よくある例

1. MEP 探索前に端点を事前最適化する。

```bash
pdb2reaction path-opt -i reactant.pdb product.pdb -q 0 -m 1 \
 --preopt --preopt-max-cycles 20000 --out-dir ./result_path_opt_preopt
```

2. GSM ではなく DMF モードで実行する。

```bash
pdb2reaction path-opt -i reactant.pdb product.pdb -q 0 -m 1 \
 --mep-mode dmf --max-nodes 12 --out-dir ./result_path_opt_dmf
```

3. リンク親原子を凍結し、climb を切って短時間で確認する。

```bash
pdb2reaction path-opt -i reactant.pdb product.pdb -q 0 -m 1 \
 --freeze-links --no-climb --out-dir ./result_path_opt_fast
```

## 使用法
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

## ワークフロー
1. **事前アライメント & 凍結解決**
 - 2 番目以降のエンドポイントは最初の構造に対して Kabsch アライメントされます。いずれかのエンドポイントで `freeze_atoms` が定義されている場合、RMSD フィットにはその原子のみを使用しますが、得られた変換は全原子に適用されます。
 - `--freeze-links` が有効な場合、リンク水素の親原子は自動的に凍結されます（{ref}`リンク水素と凍結原子 <ja-link-hydrogen-and-frozen-atoms>` を参照）。

2. **ストリング成長とHEIエクスポート**
 - 経路の成長・精密化後、内部ノード間の局所極大のうちエネルギーが最も高いものを優先的に選択します。内部の局所極大がない場合は内部ノードの最大値に、内部ノードもない場合は全体の最大値にフォールバックします。
 - 最高エネルギー画像（HEI）は `.xyz` として書き込まれます。PDB 参照がある場合は `.pdb`、Gaussian テンプレートがある場合は `.gjf` も出力します（いずれも `--convert-files` の設定に従います）。

### 主要な挙動
- **エンドポイント**: 入力は2構造のみ。形式は `geom_loader` に準拠。PDB 入力（または `--ref-pdb` 付きXYZ/GJF）で軌跡/HEIのPDB 出力が有効。
- **電荷/スピン**: 電荷の解決順序の詳細は {ref}`CLI 規約: 電荷の指定 <ja-charge-specification>` を参照してください。
- **MEPセグメント**: `--max-nodes` は内部ノード数を制御します。GSM の場合、総画像数は `max_nodes + 2`（固定端点を含む）。DMF の場合、`max_nodes` はチェーン上の移動可能なイメージ数です。GSM成長およびクライミング精密化の収束プリセットは `--thresh-stopt` または `stopt.thresh`（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`）で指定します。
- **エンドポイント事前最適化**: `--thresh` は `--opt-mode` で選ばれた単一構造最適化（`opt.lbfgs.thresh` / `opt.rfo.thresh`）のみに適用されます。
- **クライミングイメージ**: `--climb` は標準のクライミング手順とLanczos接線リファインの両方を切り替え。
- **ダンプ**: `--dump` で StringOptimizer の `stopt.dump=True` に対応し、`out_dir` 内に軌跡ダンプを出力します。リスタート YAML は YAML で有効化した場合のみ書き出されます。
- **終了コード**: 終了コードは CLI 規約の {ref}`ja-exit-codes` を参照。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH PATH` | 反応物と生成物構造 | 必須 |
| `-q, --charge INT` | 総電荷（`calc.charge`）。`.gjf` 以外では `--ligand-charge` 導出が成功しない限り必須（PDB 入力または `--ref-pdb` 付きXYZ/GJF）。`.gjf` テンプレートがあればそれを使用し、電荷メタデータが無い `.gjf` 入力は `-q` が無いと中断。両方指定時は `-q` が優先 | テンプレート/導出がない限り必須 |
| `-l, --ligand-charge TEXT` | 残基別電荷マッピング（例: `GPP:-3,SAM:1`）。PDB の残基電荷から全系の電荷を自動導出します（手動計算不要）。`-q` 省略時に使用（PDB 入力、または `--ref-pdb` 付き XYZ/GJF） | _None_ |
| `--workers`, `--workers-per-node` | UMA予測器の並列度（workers > 1 で解析ヘシアン無効; `workers_per_node` は並列予測器に渡されます）。診断上の注意は {ref}`ja-workers-fd-downgrade` を参照。 | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度 | テンプレート/`1` |
| `--freeze-links/--no-freeze-links` | PDBのみ: リンクH親を凍結（YAMLとマージ）。詳細は [extract](extract.md) を参照 | `True` |
| `--freeze-atoms TEXT` | 凍結する原子の 1 始まりインデックスをカンマ区切りで明示的に指定（例: `'1,3,5'`）。`--freeze-links` と併用可、任意の入力形式に適用。 | _None_ |
| `--max-nodes INT` | 内部ノード数。**GSM:** 総イメージ数 = `max_nodes + 2`（端点 2 つは固定）。**DMF:** チェーン上の*移動可能な*イメージ数（端点の暗黙的展開なし）。 | `20` |
| `--mep-mode {gsm\|dmf}` | GSM（ストリングベース）またはDMF（ダイレクトフラックス）経路生成器を選択 | `gsm` |
| `--max-cycles INT` | オプティマイザーマクロイテレーション上限 | `300` |
| `--climb/--no-climb` | クライミングイメージ精密化を有効化 | `True` |
| `--dump/--no-dump` | MEP軌跡/リスタートをダンプ | `False` |
| `--opt-mode TEXT` | エンドポイント事前最適化用の単一構造オプティマイザー（`grad` = LBFGS、`hess` = RFO） | `grad` |
| `--convert-files/--no-convert-files` | PDB/Gaussian入力用のXYZ/TRJ → PDB/GJFコンパニオン出力の切り替え | `True` |
| `--ref-pdb FILE` | XYZ/GJF 入力用の参照 PDB トポロジー（XYZ 座標は保持し PDB 変換を有効化） | _None_ |
| `-o, --out-dir TEXT` | 出力ディレクトリ | `./result_path_opt/` |
| `--thresh TEXT` | エンドポイント事前最適化のみの収束プリセットを上書き（`opt.lbfgs/rfo.thresh`） | `gau` |
| `--thresh-stopt TEXT` | ストリングオプティマイザーの収束プリセットを上書き（`stopt.thresh`） | `gau_loose` |
| `--config FILE` | 明示 CLI 指定より前に適用されるベース YAML | _None_ |
| `--show-config/--no-show-config` | 解決済み設定（YAML レイヤ情報を含む）を表示して実行継続 | `False` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP バックエンド | `uma` |
| `--solvent TEXT` | xTB 暗黙溶媒（例: `water`）。`none` で無効化 | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB 溶媒モデル | `alpb` |
| `--dry-run/--no-dry-run` | 実行せずに検証と実行計画表示のみを行う。 | `False` |
| `--preopt/--no-preopt` | アライメント/MEP 探索前に各エンドポイントを事前最適化（GSM/DMF）。**スコープ依存デフォルト:** 単体の `path-opt` では `False`、**`pdb2reaction all` 経由では `True` に反転**されます（{ref}`mep-search-options` を参照）。 | `False` |
| `--preopt-max-cycles INT` | エンドポイント事前最適化サイクルの上限 | `10000` |
| `--fix-ends/--no-fix-ends` | GSM成長/精密化中にエンドポイント構造を固定 | `False` |
| `--out-json/--no-out-json` | `out_dir` に機械可読な `result.json` を書き出す。スキーマは [JSON 出力スキーマ](json-output.md) を参照。 | `False` |

## 出力
```
out_dir/
├─ final_geometries_trj.xyz # XYZ経路（コメント行にエネルギーを保持）
├─ final_geometries_trj.pdb # PDB 参照が利用可能で変換が有効な場合
├─ hei.xyz # 最高エネルギー画像
├─ hei.pdb # PDB 参照が利用可能な場合のHEI（変換有効時）
├─ hei.gjf # Gaussian テンプレートを使用して書き込まれたHEI（変換有効時）
├─ align_refine/ # 剛体アライメント/リファイン段階の中間ファイル（アライメント実行時）
└─ <オプティマイザーダンプ/リスタート>
```
コンソールには解決済みYAMLブロックが出力され、GSM/DMFのMEP進行状況とタイミングが報告されます。


設定の優先順位は {ref}`CLI 規約: 設定の優先順位 <ja-configuration-precedence>` を参照してください。

```{note}
正規定義は [YAML リファレンス](yaml-reference.md) および `pdb2reaction/defaults.py` を参照してください（齟齬がある場合はそちらが優先）。以下は `path-opt` 固有のデフォルト（例: `out_dir`）を含む抜粋です。
```

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

---

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド

- [path-search](path-search.md) — 自動精密化を伴う再帰的MEP 探索（2+構造用）
- [tsopt](tsopt.md) — HEI を TS 候補として最適化（内部で虚振動数チェック済み）。続けて IRC で接続性を確認
- [extract](extract.md) — path-opt入力用の活性部位モデルPDBを生成
- [all](all.md) — 一気通貫ワークフロー（デフォルトで再帰的 path-search を使用; `--refine-path False` で path-opt に切替。`--refine-path` フラグは `pdb2reaction all` にのみ属します — 定義は {ref}`mep-search-options` を参照してください）
- [YAML リファレンス](yaml-reference.md) — `gs`、`dmf`、`stopt`、`opt` の完全な設定オプション
- [用語集](glossary.md) — MEP、GSM、DMF、HEIの定義
