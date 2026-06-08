# `path-search`

**2 構造以上**から、GSM（デフォルト）または DMF（`--mep-mode dmf`）で連続的な最小エネルギー経路（MEP）を構築します。共有結合変化が検出される領域のみを選択的に精密化し、解決済みのサブパスを連結して 1 本の軌跡にまとめ、各セグメントの最高エネルギー画像（HEI）を TS 候補として出力します（tsopt + IRC で検証）。再帰的分解により多段階反応を自動検出し、各素反応ステップの詳細な MEP を構築しますが、複雑な多段階反応の機構を満足な経路として得るには、入力中間体やスキャン仕様、収束閾値の調整など手動での試行錯誤が必要になることがあります。

## 使いどころ

- R → … → P の **2 構造以上**を入力とし、自動精密化を含む連続 MEP を 1 本の軌跡にまとめたい場面。
- セグメント生成器として `--mep-mode gsm`（デフォルト、string ベース）または `--mep-mode dmf`（direct flux）を選びます。
- 精密化シードとして `--refine-mode peak`（HEI±1 を最適化）または `--refine-mode minima`（最寄り局所極小点へ外側探索）を選びます（未指定時は GSM で `peak`、DMF で `minima`）。
- **2 端点だけ**で再帰精密化が不要な場合は、[path-opt](path-opt.md) の方がシンプルです。

## 実行例

```bash
pdb2reaction path-search -i reactant.pdb product.pdb -q 0 -m 1 \
 --out-dir ./result_path_search
```

```bash
# 中間体を明示して多段の経路を与える
pdb2reaction path-search -i R.pdb IM1.pdb IM2.pdb P.pdb -q -1 -m 1 \
 --out-dir ./result_path_search_multi
```

```bash
# テンプレート参照を使って全系マージ出力を有効化する
pdb2reaction path-search -i R.pdb IM1.pdb P.pdb -q 0 -m 1 \
 --ref-full-pdb holo_template.pdb --out-dir ./result_path_search_merge
```

```bash
# DMF + minima リファインで探索する
pdb2reaction path-search -i reactant.pdb product.pdb -q 0 -m 1 \
 --mep-mode dmf --refine-mode minima --out-dir ./result_path_search_dmf
```

## 入力

コマンド形式:

```bash
pdb2reaction path-search -i R.pdb [I.pdb ...] P.pdb [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [--multiplicity 2S+1] \
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

| 入力 | 必須 | 備考 |
| --- | --- | --- |
| `-i, --input` | はい | 反応順序の 2 つ以上の構造（反応物 → 生成物）。すべてのファイルを単一の `-i`/`--input` の後ろに並べて指定。 |
| `-q, --charge` | 条件付き | 総電荷。非 `.gjf` 入力では `--ligand-charge/-l` の導出が成功しない限り必須（PDB 入力）。両方指定時は `-q` が `--ligand-charge/-l` より優先。 |
| `-l, --ligand-charge` | 任意 | スカラー整数のリガンド総電荷、または残基別マッピング（例: `GPP:-3,SAM:1`）。`-q` 省略時に使用（PDB 入力のみ。XYZ/GJF は `-q` 必須）。 |
| `-m, --multiplicity` | いいえ | スピン多重度（2S+1）。既定は `.gjf` テンプレート値または `1`。 |
| `--ref-pdb` | XYZ/GJF マージ時 | 入力が XYZ/GJF の場合に最終的な全系マージで用いる活性部位モデル参照 PDB（入力と同数・同順）。 |

`--convert-files` が有効（デフォルト）な場合、参照 PDB があれば軌跡の `.pdb` コンパニオンを、Gaussian テンプレートがあれば HEI スナップショットの `.gjf` コンパニオンを生成します。XYZ/GJF 入力では `--ref-pdb` が最終的な全系マージで用いるポケット参照 PDB（入力と同数・同順）を提供し、`--ref-full-pdb` によりフルテンプレートへのマージが可能です（XYZ/GJF 入力では主軌跡の PDB コンパニオンは生成されません）。

## 処理の流れ

1. **ペアごとの初期セグメント（GSM/DMF）** – 各隣接入力（A→B）間で `GrowingString` または DMF を実行し、粗い MEP と最高エネルギー画像（HEI）を取得。
2. **HEI 周辺の局所緩和** – `refine-mode=peak` なら HEI±1、`refine-mode=minima` なら HEI 近傍の局所極小点を、選択した単一構造オプティマイザー（`opt-mode`）で精密化し `End1`/`End2` を得る。
   > **デフォルト:** `--refine-mode` 省略時は GSM では `peak`、DMF では `minima` が選択されます。
3. **ねじれ vs. 精密化の決定** – `End1` と `End2` 間に共有結合変化がなければ *ねじれ*（kink: 共有結合変化を伴わない構造変化区間。[用語集](glossary.md) 参照）とみなし、`search.kink_max_nodes` の線形ノードを挿入して個別最適化。結合変化がある場合は *反応セグメント*（端点間に共有結合変化が検出される区間。[用語集](glossary.md) 参照）として扱い、`End1` と `End2` 間に **精密化セグメント (GSM/DMF)** を起動して障壁を鋭利化。
4. **選択的再帰** – `(A→End1)` と `(End2→B)` の結合変化を `bond` しきい値で比較し、共有結合更新が残るサブ区間のみ再帰的に探索。再帰深度は `search.max_depth` で制限。
5. **スティッチング & ブリッジング** – 解決済みのサブパスを連結し、RMSD ≤ `search.stitch_rmsd_thresh` の重複エンドポイントを除去。RMSD ギャップが `search.bridge_rmsd_thresh` を超える場合は *ブリッジセグメント*（非隣接の中間体間を接続するセグメント。[用語集](glossary.md) 参照）を GSM/DMF で挿入。境界で結合変化が検出される場合はブリッジではなく新規の再帰セグメントで置換。
6. **アライメント & マージング（オプション）** – `--align`（デフォルト）で事前最適化構造を先頭入力へ剛体アライメントし、`freeze_atoms` を整合。`--ref-full-pdb` を指定すると活性部位モデル軌跡をフルサイズ PDB テンプレートへマージ（`--align` により先頭テンプレートの再利用が可能）。

結合変化の判定は `bond_changes.compare_structures` を用い、`bond` セクションのしきい値に従います。MLIP バックエンドは全構造で共有・再利用されます。

## 出力

```
out_dir/ (デフォルト:./result_path_search/)
├─ mep_trj.xyz # 主要 MEP 軌跡
├─ mep.pdb # 入力がPDB テンプレートで変換が有効な場合のPDB コンパニオン
├─ mep.gjf # Gaussian テンプレート検出時の Gaussian コンパニオン
├─ mep_w_ref.pdb # マージされた全系MEP（参照 PDB/テンプレートが必要）
├─ mep_seg_XX_trj.xyz # セグメントごとの MEP 軌跡（XYZ）
├─ mep_seg_XX.pdb # セグメントごとの PDB コンパニオン（変換有効時）
├─ mep_seg_XX.gjf # セグメントごとの Gaussian コンパニオン（テンプレート検出時）
├─ mep_w_ref_seg_XX.pdb # 共有結合変化がある場合のマージされたセグメントごとのパス
├─ hei_seg_XX.xyz # セグメントごとの最高エネルギー画像
├─ hei_seg_XX.pdb # HEI PDB コンパニオン（変換有効時）
├─ hei_seg_XX.gjf # HEI Gaussian コンパニオン（テンプレート検出時）
├─ hei_w_ref_seg_XX.pdb # 全系コンテキストでマージされた HEI（参照 PDB が必要）
├─ summary.json # すべての再帰セグメントの障壁と分類サマリー
├─ summary.log # 結果要約
├─ mep_plot.png # `trj2fig` で生成した ΔE プロファイル（kcal/mol、反応物基準）
├─ energy_diagram_MEP.png # MEP 状態エネルギーダイアグラムの静的エクスポート（反応物基準）
└─ seg_NNN_*/ # セグメントごとの GSM/DMF ダンプ、HEI スナップショット、kink/精密化の診断情報
```

結果は通常、次のファイルを開いて確認します。

- `result_path_search/mep_trj.xyz`
- `result_path_search/summary.json`
- `result_path_search/summary.log`
- `result_path_search/mep_plot.png`（プロット生成時）

- コンソールには確定済みの設定ブロック（`geom`, `calc`, `gs`, `stopt`, `opt.*`, `bond`, `search`）が出力されます。詳細は {ref}`ja-verbosity-levels` を参照してください。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH...` | 反応順序の 2 つ以上の構造（反応物 → 生成物）。すべてのファイルを単一の `-i`/`--input` の後ろに並べて指定 | 必須 |
| `-q, --charge INT` | 総電荷。非 `.gjf` 入力では `--ligand-charge` の導出が成功しない限り必須。両方指定時は `-q` が優先 | テンプレート/導出が適用されない限り必須 |
| `-l, --ligand-charge TEXT` | スカラー整数（例: `-1`）でリガンド総電荷を指定するか、残基別マッピング（例: `GPP:-3,SAM:1`）で PDB 残基電荷から全系の電荷を導出。`-q` 省略時に使用（PDB 入力のみ。XYZ/GJF は `-q` 必須） | _None_ |
| `--workers`, `--workers-per-node` | MLIP 予測器の並列度（workers > 1 で解析ヘシアン無効; UMA バックエンドのみ; `workers_per_node` は並列予測器に渡されます）。診断上の注意は {ref}`ja-workers-fd-downgrade` を参照 | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度（2S+1） | `.gjf` テンプレート値または `1` |
| `--freeze-links/--no-freeze-links` | PDB 活性部位モデル読み込み時、リンク水素の親原子を凍結。詳細は [extract](extract.md) を参照 | `True` |
| `--freeze-atoms TEXT` | 凍結する原子の 1 始まりインデックスをカンマ区切りで明示的に指定（例: `'1,3,5'`）。`--freeze-links` と併用可、任意の入力形式に適用 | _None_ |
| `--max-nodes INT` | MEP セグメントごとの内部ノード（GSM string image または DMF image） | `20` |
| `--max-cycles INT` | 最大 MEP 最適化サイクル（GSM/DMF） | `300` |
| `--climb/--no-climb` | GSM セグメントのクライミングイメージを有効化（ブリッジは無効） | `True` |
| `--opt-mode TEXT` | HEI±1/ねじれノード用の単一構造オプティマイザー（`grad`=L-BFGS、`hess`=RFO）。同じトークンが `tsopt` では Dimer / RS-I-RFO へ対応する点については {ref}`ja-opt-mode-semantics` を参照してください | `grad` |
| `--mep-mode {gsm\|dmf}` | セグメント生成器: GSM（string）または DMF（direct flux） | `gsm` |
| `--refine-mode {peak\|minima}` | 精密化シード: `peak` は HEI±1、`minima` は HEI から最寄り局所極小点へ外側探索。未指定時は GSM で `peak`、DMF で `minima` | _Auto_ |
| `--dump/--no-dump` | MEP（GSM/DMF）と単一構造軌跡をダンプ。リスタート YAML は YAML で有効化した場合のみ書き出されます | `False` |
| `--convert-files/--no-convert-files` | PDB/Gaussian 入力の XYZ/TRJ → PDB/GJF コンパニオンを切り替え | `True` |
| `-o, --out-dir TEXT` | 出力ディレクトリ | `./result_path_search/` |
| `--thresh TEXT` | 単一構造最適化のみの収束プリセットを上書き（`opt.lbfgs/rfo.thresh`） | `gau` |
| `--thresh-stopt TEXT` | ストリングオプティマイザーの収束プリセットを上書き（`stopt.thresh`） | `gau_loose` |
| `--config FILE` | 明示 CLI 指定より前に適用されるベース YAML | _None_ |
| `--show-config/--no-show-config` | 解決済み設定（YAML レイヤ情報を含む）を表示して実行継続 | `False` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP バックエンド | `uma` |
| `--solvent TEXT` | xTB 暗黙溶媒（例: `water`）。`none` で無効化 | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB 溶媒モデル | `alpb` |
| `--dry-run/--no-dry-run` | 実行せずに検証と実行計画表示のみを行う | `False` |
| `--preopt/--no-preopt` | 選択された単一構造オプティマイザー（L-BFGS/RFO）で MEP 探索前に各エンドポイントを事前最適化。 | `True` |
| `--align/--no-align` | 探索前にすべての入力を最初の構造にアライメント | `True` |
| `--ref-full-pdb PATH...` | フルサイズテンプレート PDB（入力と同数。`--align` があれば先頭のみ再利用可） | _None_ |
| `--ref-pdb PATH...` | 入力が XYZ/GJF の場合に最終的な全系マージで用いるポケット参照 PDB（入力と同数・同順） | _None_ |

設定の優先順位は {ref}`CLI 規約: 設定の優先順位 <ja-configuration-precedence>` を参照してください。

## YAML 設定

YAML ルートはマッピングでなければなりません。共通セクションは [YAML リファレンス](yaml-reference.md) を再利用します: `geom`/`calc` は単一構造設定を反映し（PDB 入力では `--freeze-links` が `geom.freeze_atoms` を補強します。詳細は {ref}`リンク水素と凍結原子 <ja-link-hydrogen-and-frozen-atoms>` を参照）、`stopt` は `path-opt`（[path-opt.md](path-opt.md)）に記載の StringOptimizer 設定を継承します。

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

## 注意事項

- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。
- 入力は 2 つ以上が必須。満たさない場合、`-i/--input` の "invalid value" エラーで終了します。
- 複数テンプレートを渡す場合は `--ref-full-pdb` をファイルごとに繰り返して指定します。`--align` が有効な場合、マージでは先頭テンプレートのみが再利用されます。
- MLIP バックエンドは全構造で共有・再利用されます。
- `--dump` が有効な場合、MEP（GSM/DMF）と単一構造最適化の軌跡が出力されます。リスタート YAML は YAML で `dump_restart` を有効にした場合のみ書き出されます。

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
- [path-opt](path-opt.md) — 単一パス MEP 最適化（再帰的精密化なし）
- [tsopt](tsopt.md) — HEI を遷移状態として最適化
- [extract](extract.md) — path-search 入力用の活性部位モデル PDB を生成
- [all](all.md) — 一気通貫ワークフロー（デフォルトで再帰的 path-search を使用; `--refine-path False` で path-opt に切替）
- [YAML リファレンス](yaml-reference.md) — `gs`、`dmf`、`bond`、`search` の完全な設定オプション
- [用語集](glossary.md) — MEP、GSM、DMF、HEI の定義
