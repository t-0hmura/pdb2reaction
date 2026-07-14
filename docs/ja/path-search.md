# `path-search`

R → … → P の **2 構造以上**から、連続的な最小エネルギー経路（MEP）を構築します。自動精密化を含む連続 MEP を 1 本の軌跡にまとめたい場面で使用します。**2 端点だけ**で再帰精密化が不要な場合は、[path-opt](path-opt.md) の方がシンプルです。

経路は 2 つのエンジンのいずれかで生成します。GSM（デフォルト、`--mep-mode gsm`、string ベース）か、DMF（`--mep-mode dmf`、direct flux）です。共有結合変化が検出される領域のみを選択的に精密化します（`--refine-mode peak` は HEI±1 を最適化、`--refine-mode minima` は最寄り局所極小点へ外側探索、デフォルトは GSM で `peak`、DMF で `minima`）。解決済みのサブパスを連結して 1 本の軌跡にまとめ、各セグメントの最高エネルギー画像（HEI）を TS 候補として出力します（tsopt + IRC で検証）。

再帰的分解は結合変化の情報を使い、多段階反応の候補となるより狭い反応セグメントを提案します。これは幾何構造に基づくヒューリスティクであり、各セグメントが 1 つの素反応または 1 つの TS を含むことを証明しません。各 HEI を `tsopt`、虚振動 1 本の確認、IRC 接続性で必ず検証してください。複雑な機構では、入力中間体、スキャン仕様、収束閾値の手動調整が必要な場合があります。

## 実行例

コマンド形式:

```bash
pdb2reaction path-search -i R.pdb [-i I.pdb ...] -i P.pdb [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [--multiplicity 2S+1] \
 [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] \
 [--workers N] [--workers-per-node N] \
 [--mep-mode {gsm|dmf}] [--freeze-links/--no-freeze-links] [--thresh PRESET] [--thresh-stopt PRESET] \
 [--refine-mode {peak|minima}] \
 [--max-nodes N] [--max-cycles N] [--climb/--no-climb] \
 [--opt-mode grad|hess] [--dump/--no-dump] \
 [--out-dir DIR] [--preopt/--no-preopt] \
 [--align/--no-align] [--ref-full-pdb FILE...] [--ref-pdb FILE...] \
 [--convert-files/--no-convert-files] \
 [--show-config/--no-show-config] [--dry-run/--no-dry-run]
```

2 端点（反応物 → 生成物）:

```bash
pdb2reaction path-search -i reactant.pdb -i product.pdb -q 0 -m 1 \
 --out-dir ./result_path_search
```

中間体を明示して多段の経路を与える:

```bash
# 中間体を明示して多段の経路を与える
pdb2reaction path-search -i R.pdb -i IM1.pdb -i IM2.pdb -i P.pdb -q -1 -m 1 \
 --out-dir ./result_path_search_multi
```

テンプレート参照を使って全系マージ出力を有効化する:

```bash
# テンプレート参照を使って全系マージ出力を有効化する
pdb2reaction path-search -i R.pdb -i IM1.pdb -i P.pdb -q 0 -m 1 \
 --ref-full-pdb holo_template.pdb --out-dir ./result_path_search_merge
```

DMF + minima 精密化で探索する:

```bash
# DMF + minima 精密化で探索する
pdb2reaction path-search -i reactant.pdb -i product.pdb -q 0 -m 1 \
 --mep-mode dmf --refine-mode minima --out-dir ./result_path_search_dmf
```

## 処理の流れ

1. **ペアごとの初期セグメント（GSM/DMF）** – 各隣接入力（A→B）間で `GrowingString` または DMF を実行し、粗い MEP と最高エネルギー画像（HEI）を取得。
2. **HEI 周辺の局所緩和** – `refine-mode=peak` なら HEI±1、`refine-mode=minima` なら HEI 近傍の局所極小点を、選択した単一構造オプティマイザ（`opt-mode`）で精密化し `End1`/`End2` を得る。
   > **デフォルト:** `--refine-mode` 省略時は GSM では `peak`、DMF では `minima` が選択されます。
3. **ねじれ vs. 精密化の決定** – `End1` と `End2` 間に共有結合変化がなければ *ねじれ*（kink: 共有結合変化を伴わない構造変化区間。[用語集](glossary.md) 参照）とみなし、`search.kink_max_nodes` の線形ノードを挿入して個別最適化。結合変化がある場合は *反応セグメント*（端点間に共有結合変化が検出される区間。[用語集](glossary.md) 参照）として扱い、`End1` と `End2` 間に **精密化セグメント (GSM/DMF)** を起動して障壁を先鋭化。
4. **選択的再帰** – `(A→End1)` と `(End2→B)` の結合変化を `bond` しきい値で比較し、共有結合更新が残るサブ区間のみ再帰的に探索。再帰深度は `search.max_depth` で制限。
5. **スティッチング & ブリッジング** – 解決済みのサブパスを連結し、RMSD ≤ `search.stitch_rmsd_thresh` の重複エンドポイントを除去。RMSD ギャップが `search.bridge_rmsd_thresh` を超える場合は *ブリッジセグメント*（非隣接の中間体間を接続するセグメント。[用語集](glossary.md) 参照）を GSM/DMF で挿入。境界で結合変化が検出される場合はブリッジではなく新規の再帰セグメントで置換。
6. **アライメント & マージング（オプション）** – `--align`（デフォルト）で事前最適化構造を先頭入力へ剛体アライメントし、`freeze_atoms` を整合。`--ref-full-pdb` を指定すると活性部位モデル軌跡をフルサイズ PDB テンプレートへマージ（`--align` により先頭テンプレートの再利用が可能）。

結合変化の判定は `bond_changes.compare_structures` を用い、`bond` セクションのしきい値に従います。MLIP バックエンドは全構造で共有・再利用されます。

## 出力

```
out_dir/ (デフォルト:./result_path_search/)
├─ mep_trj.xyz # 主要 MEP 軌跡
├─ mep.pdb # PDB/mmCIF topology入力で変換有効時
├─ mep.cif # mmCIF/oversized-PDB入力。元IDを復元
├─ mep.gjf # Gaussian テンプレート検出時に対応する Gaussian
├─ mep_w_ref.pdb # マージされた全系MEP（参照topologyが必要）
├─ mep_w_ref.cif # bridge templateを元IDで復元
├─ mep_seg_XX_trj.xyz # セグメントごとの MEP 軌跡（XYZ）
├─ mep_seg_XX.pdb # セグメントごとのPDB（変換有効時）
├─ mep_seg_XX.cif # bridge入力のCIF
├─ mep_seg_XX.gjf # セグメントごとに対応する Gaussian（テンプレート検出時）
├─ mep_w_ref_seg_XX.pdb # 共有結合変化がある場合のマージ済みsegment
├─ mep_w_ref_seg_XX.cif # bridge templateのCIF
├─ hei_seg_XX.xyz # セグメントごとの最高エネルギー画像
├─ hei_seg_XX.pdb # HEI に対応する PDB（変換有効時）
├─ hei_seg_XX.cif # bridge入力のCIF
├─ hei_seg_XX.gjf # HEI に対応する Gaussian（テンプレート検出時）
├─ hei_mode_seg_XX.txt # HEI の energy-upwinding Cartesian 接線
├─ hei_w_ref_seg_XX.pdb # 全系コンテキストでマージされた HEI
├─ hei_w_ref_seg_XX.cif # bridge templateのCIF
├─ summary.json # すべての再帰セグメントの障壁と分類サマリー
├─ summary.log # 結果要約
├─ mep_plot.png # `trj2fig` で生成した ΔE プロファイル（kcal/mol、反応物基準）
├─ energy_diagram_MEP.png # MEP 状態エネルギーダイアグラムの静止画出力（反応物基準）
└─ seg_NNN_*/ # セグメントごとの GSM/DMF ダンプ、HEI スナップショット、kink/精密化の診断情報
```

結果は通常、次のファイルを開いて確認します。

- `result_path_search/mep_trj.xyz`
- `result_path_search/summary.json`
- `result_path_search/summary.log`
- `result_path_search/mep_plot.png`（プロット生成時）

- コンソールには確定済みの設定ブロック（`geom`, `calc`, `gs`, `stopt`, `opt.*`, `bond`, `search`）が出力されます。詳細は {ref}`ja-verbosity-levels` を参照してください。

## CLI オプション

完全なフラグ一覧は生成された [コマンドリファレンス](../reference/commands/index.md) を参照してください。以下の表は説明が必要なオプションのみを扱います。

表は目的ごとにグループ化しており、各グループ内では使用頻度の高いオプションを先に並べています。

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| **入力と電荷** | | |
| `-i, --input PATH` | 反応順序の 2 つ以上の構造（反応物 → 生成物）。各ファイルごとに `-i`/`--input` を繰り返すか、単一の `-i` の後ろに複数ファイルを並べる（例: `-i R.pdb -i IM1.pdb -i P.pdb` または `-i R.pdb IM1.pdb P.pdb`） | 必須 |
| `-q, --charge INT` | 総電荷。非 `.gjf` 入力では `--ligand-charge` の導出が成功しない限り必須。両方指定時は `-q` が優先 | テンプレート/導出が適用されない限り必須 |
| `-l, --ligand-charge TEXT` | 単一の整数（例: `-1`）または残基別マッピング（例: `GPP:-3,SAM:1`）から PDB/mmCIF トポロジーの全系電荷を導出。裸のXYZでは使用不可だが、`--ref-pdb`でトポロジーを付与すれば使用可。正しいGJFはヘッダーの電荷・多重度を自動継承する | _None_ |
| `-m, --multiplicity INT` | スピン多重度（2S+1） | `.gjf` テンプレート値または `1` |
| **バックエンドと計算** | | |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP バックエンド | `uma` |
| `--workers`, `--workers-per-node` | UMA 予測器の並列度（`workers_per_node` は並列予測器へ転送）。`workers > 1` と明示的な解析 Hessian は併用不可。{ref}`ja-workers-analytical-error` を参照 | `1`, `1` |
| `--solvent TEXT` | xTB 暗黙溶媒（例: `water`）。`none` で無効化 | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB 溶媒モデル | `alpb` |
| **活性領域の凍結** | | |
| `--freeze-links/--no-freeze-links` | PDB 活性部位モデル読み込み時、キャップ水素の親原子を凍結。詳細は [extract](extract.md) を参照 | `True` |
| `--freeze-atoms TEXT` | 凍結する原子の 1 始まりインデックスをカンマ区切りで明示的に指定（例: `'1,3,5'`）。`--freeze-links` と併用可、任意の入力形式に適用 | _None_ |
| **MEP 探索** | | |
| `--mep-mode {gsm\|dmf}` | セグメント生成器: GSM（string）または DMF（direct flux） | `gsm` |
| `--dmf-backend {cpu\|gpu}` | DMF 計算バックエンド（`--mep-mode dmf` 時のみ）: `gpu`（`dmf.torch`/CUDA）または `cpu`（`dmf`/NumPy）。GPU メモリ不足時は `cpu` で再実行 | `gpu` |
| `--preopt/--no-preopt` | 選択された単一構造オプティマイザ（L-BFGS/RFO）で MEP 探索前に各エンドポイントを事前最適化。 | `True` |
| `--max-nodes INT` | MEP セグメントごとの内部ノード（GSM string image または DMF image） | `20` |
| `--max-cycles INT` | 最大 MEP 最適化サイクル（GSM/DMF） | `300` |
| `--climb/--no-climb` | GSM セグメントのクライミングイメージを有効化（ブリッジは無効） | `True` |
| **精密化** | | |
| `--refine-mode {peak\|minima}` | 精密化シード: `peak` は HEI±1、`minima` は HEI から最寄り局所極小点へ外側探索。未指定時は GSM で `peak`、DMF で `minima` | _Auto_ |
| `--opt-mode TEXT` | HEI±1/ねじれノード用の単一構造オプティマイザ（`grad`=L-BFGS、`hess`=RFO）。同じトークンが `tsopt` では Dimer / RS-P-RFO へ対応する点については {ref}`ja-opt-mode-semantics` を参照してください | `grad` |
| **収束閾値** | | |
| `--thresh TEXT` | 単一構造最適化のみの収束プリセットを上書き（`opt.lbfgs/rfo.thresh`） | `gau` |
| `--thresh-stopt TEXT` | ストリングオプティマイザの収束プリセットを上書き（`stopt.thresh`） | `gau_loose` |
| **マージとアライメント** | | |
| `--align/--no-align` | 探索前にすべての入力を最初の構造にアライメント | `True` |
| `--ref-full-pdb PATH...` | フルサイズテンプレート PDB（入力と同数。`--align` があれば先頭のみ再利用可） | _None_ |
| `--ref-pdb PATH...` | XYZ/GJF入力の全系マージに使うactive-site PDB/mmCIF参照（入力と同数・同順） | _None_ |
| **出力と設定** | | |
| `-o, --out-dir TEXT` | 出力ディレクトリ | `./result_path_search/` |
| `--dump/--no-dump` | MEP（GSM/DMF）と単一構造軌跡をダンプ。リスタート YAML は YAML で有効化した場合のみ書き出されます | `False` |
| `--convert-files/--no-convert-files` | XYZ/TRJ → PDB/CIF/GJFを切り替え。bridge入力は元IDのCIFを追加し、XYZ/GJFは参照topologyなしではPDBを生成しません。 | `True` |
| `--config FILE` | 明示 CLI 指定より前に適用されるベース YAML | _None_ |
| `--show-config/--no-show-config` | 解決済み設定（YAML レイヤ情報を含む）を表示して実行継続 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに検証と実行計画表示のみを行う | `False` |

設定の優先順位は {ref}`CLI 規約: 設定の優先順位 <ja-configuration-precedence>` を参照してください。

## YAML 設定

YAML ルートはマッピングでなければなりません。共通セクションは [YAML リファレンス](yaml-reference.md) を再利用します: `geom`/`calc` は単一構造設定を反映し（PDB/mmCIF トポロジーでは `--freeze-links` が `geom.freeze_atoms` を補強します。詳細は {ref}`キャップ水素と凍結原子 <ja-link-hydrogen-and-frozen-atoms>` を参照）、`stopt` は `path-opt`（[path-opt.md](path-opt.md)）に記載の StringOptimizer 設定を継承します。

`bond` と `search` は `path-search` の再帰ロジックの中核であり、ここで詳述します。`gs`、`dmf`、`stopt`、`opt.lbfgs`、`opt.rfo` は `path-search` 固有の `out_dir` 上書きのみ再掲します。

`bond` は MLIP ベースの結合変化検出パラメータで、[scan](scan.md) の bond セクションと共通の `device`, `bond_factor`, `margin_fraction`, `delta_fraction` を持ちます。

### `path-search` 固有の上書き

```yaml
stopt:
 out_dir: ./result_path_search/ # path-search 上書き（正準デフォルト: ./result_path_opt/）
opt:
 lbfgs:
   out_dir: ./result_path_search/ # path-search 上書き（正準デフォルト: ./result_opt/）
 rfo:
   out_dir: ./result_path_search/ # path-search 上書き（正準デフォルト: ./result_opt/）
```

`bond` / `search`:

```yaml
bond:
 device: auto # MLIP device for bond analysis
 bond_factor: 1.2 # covalent-radius scaling
 margin_fraction: 0.05 # tolerance margin for comparisons
 delta_fraction: 0.05 # minimum relative change to flag bonds
search:
 max_depth: 10 # recursion depth limit
 stitch_rmsd_thresh: 0.0001 # RMSD threshold for stitching segments
 bridge_rmsd_thresh: 0.0001 # RMSD threshold for bridging nodes
 max_nodes_segment: 20 # max nodes per segment
 max_nodes_bridge: 5 # max nodes per bridge
 kink_max_nodes: 3 # max nodes for kink optimizations
 max_seq_kink: 2 # max sequential kinks
 refine_mode: null # optional refinement strategy (auto-chooses when null)
```

## 注記

- 入力は 2 つ以上が必須。満たさない場合、`-i/--input` の "invalid value" エラーで終了します。

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) — 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) — 詳細な対処ガイド
- [path-opt](path-opt.md) — 単一パス MEP 最適化（再帰的精密化なし）
- [tsopt](tsopt.md) — HEI を遷移状態として最適化
- [extract](extract.md) — path-search 入力用の活性部位モデル PDB を生成
- [all](all.md) — 一気通貫ワークフロー（デフォルトで単一パス path-opt を使用; `--refine-path` で再帰的 path-search に切替）
- [YAML リファレンス](yaml-reference.md) — `gs`、`dmf`、`bond`、`search` の完全な設定オプション
- [用語集](glossary.md) — MEP、GSM、DMF、HEI の定義
