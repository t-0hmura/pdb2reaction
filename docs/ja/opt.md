# `opt`

このコマンドは pysisyphus の L-BFGS（`lbfgs`）または RFOptimizer（`rfo`）を用い、MLIP（デフォルト: UMA、`-b/--backend` で ORB・MACE・AIMNet2 も選択可能）のエネルギー・勾配・Hessian で単一構造を局所極小点へ最適化します。距離拘束や虚振動数モードのフラット化も任意で併用できます。入力は PDB/mmCIF、XYZ、GJF の単一 geometry とし、trajectory は使用する frame を `.xyz` へ抽出してから渡します。L-BFGS 最小化には `--opt-mode grad`（alias `lbfgs`、デフォルト）、RFOptimizer には `--opt-mode hess`（alias `rfo`）を選択します。

## 実行例

コマンド形式:

```bash
pdb2reaction opt -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m MULT] \
 [-b/--backend uma|orb|mace|aimnet2] \
 [--opt-mode grad|hess|lbfgs|rfo] [--flatten/--no-flatten] [--freeze-links/--no-freeze-links] \
 [--dist-freeze '[(i,j,target_Å),...]'] [--one-based|--zero-based] \
 [--bias-k K_eV_per_Å²] [--dump/--no-dump] [-o/--out-dir DIR] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE]
```

基本的な最小化:

```bash
pdb2reaction opt -i input.pdb -q 0 -m 1 --out-dir ./result_opt
```

収束を厳しくして軌跡ダンプを保存する:

```bash
pdb2reaction opt -i input.pdb -q 0 -m 1 --thresh gau_tight --dump \
 --out-dir ./result_opt_tight
```

調和距離拘束を追加する。例では `--bias-k 20.0`（目標距離付近でゆるく誘導する弱い拘束）を使っていますが、`bias.k` のデフォルトは 300 eV·Å⁻² で、最適化中に拘束を支配的にしたい場合はデフォルト値の方が適しています:

```bash
pdb2reaction opt -i input.pdb -q 0 -m 1 \
 --dist-freeze '[(1,5,2.0)]' --bias-k 20.0 --out-dir ./result_opt_rest
# 2-tuple 形式: 原子 1-5 間の距離を現在値に固定 --dist-freeze '[(1,5)]'
```

RFO モードを明示して実行する:

```bash
pdb2reaction opt -i input.pdb -q 0 -m 1 --opt-mode hess \
 --out-dir ./result_opt_hess
```

## 処理の流れ

- **オプティマイザ**: `--opt-mode grad`（alias: `lbfgs`、デフォルト）→ L-BFGS、`--opt-mode hess`（alias: `rfo`）→ RFOptimizer。サブコマンド別のトークン→アルゴリズム対応は {ref}`ja-opt-mode-semantics` を参照。
  > **命名規則の注意:** CLI は `grad|lbfgs` および `hess|rfo` を受け付けます。YAML では `lbfgs` または `rfo` を直接指定してください。
- **Flatten loop**: `--flatten` を有効にすると、最適化後に虚振動数モードのフラット化ループを実行します。`opt` では各反復で検出された虚振動数モードをすべてフラット化し、虚振動数が残らなくなるか内部ループ上限に達するまで繰り返します。PHVAの固有値解析は、凍結anchorを尊重する constrained 処理に固定されています。
- **拘束**: `--dist-freeze` は Python リテラルタプル `(i, j, target_Å)` を解釈します（`target_Å` は目標距離、単位は Å）。3 番目の要素を省略すると開始距離を拘束します。`--bias-k` はグローバル調和強度（eV·Å⁻²）を設定します。インデックスはデフォルトで 1 始まりですが、`--zero-based` で 0 始まりに切り替えられます。
- **電荷/スピン解決**: 電荷の解決順序の詳細は {ref}`CLI 規約: 電荷の指定 <ja-charge-specification>` を参照してください。
- **凍結原子**: `--freeze-links` が有効な場合、キャップ水素の親原子は自動的に凍結されます（{ref}`キャップ水素と凍結原子 <ja-link-hydrogen-and-frozen-atoms>` を参照）。
- **ダンプ & 変換**: `--dump` は `opt.dump=True` を反映し `optimization_trj.xyz` を出力します。変換が有効な場合、PDB入力はPDB companion、mmCIF／oversized-PDB入力はPDBと元IDを復元したCIFを出力します。`opt.dump_restart` を有効にするとリスタート YAML が出力されます。
- **終了コード**: 終了コードは CLI 規約の {ref}`ja-exit-codes` を参照。

## 出力

```
out_dir/
├─ final_geometry.xyz # 常に出力
├─ final_geometry.pdb # PDB/mmCIF topology入力、変換有効時
├─ final_geometry.cif # mmCIF/oversized-PDB入力、変換有効時
├─ final_geometry.gjf # Gaussian テンプレートが検出され変換が有効な場合
├─ optimization_trj.xyz # ダンプが有効な場合のみ
├─ optimization.pdb # 軌跡のPDB変換（topology入力、変換有効時）
├─ optimization.cif # bridge入力。元IDを復元した軌跡
└─ restart*.yml # opt.dump_restart 設定時のリスタートファイル（任意）
```
コンソールには解決済みの `geom`/`calc`/`opt`/`lbfgs`/`rfo` ブロックとサイクル進行、総実行時間が出力されます。

最適化後は、主な成果物として次を確認します。

- `result_opt/final_geometry.xyz`
- `result_opt/final_geometry.pdb`（PDB 入力かつ変換有効時）
- `result_opt/optimization_trj.xyz`（`--dump` 有効時）

設定の優先順位は {ref}`CLI 規約: 設定の優先順位 <ja-configuration-precedence>` を参照してください。

## CLI オプション

完全なフラグ一覧は生成された [コマンドリファレンス](../reference/commands/index.md) にあります。以下の表は説明が必要なオプションのみを扱い、`--backend-model`・`--precision` などを含む網羅的な一覧はここでは重複して記載していません。

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | PDB/mmCIF、XYZ、GJF の単一 geometry。trajectory は使用する frame を `.xyz` へ抽出してから指定 | 必須 |
| `-q, --charge INT` | 総電荷。`.gjf` テンプレートまたは `--ligand-charge`（PDB/mmCIF 入力または `--ref-pdb` 付き XYZ/GJF）が提供しない限り必須。両方指定時は `-q` が優先 | テンプレート/導出が適用されない限り必須 |
| `-l, --ligand-charge TEXT` | 単一の整数（例: `-1`）でリガンド総電荷を指定するか、残基別マッピング（例: `GPP:-3,SAM:1`）で PDB/mmCIF 残基電荷から全系の電荷を導出。`-q` 省略時に使用（PDB/mmCIF 入力、または `--ref-pdb` 付き XYZ/GJF） | _None_ |
| `--workers INT` | UMA 予測器の並列度。`workers > 1` と明示的な解析 Hessian は併用できないため、`workers = 1` または有限差分を使用。{ref}`ja-workers-analytical-error` を参照 | `1` |
| `--workers-per-node INT` | ノードあたりのワーカー数。並列予測器に渡されます | `1` |
| `-m, --multiplicity INT` | スピン多重度（2S+1）。`.gjf` テンプレートまたは `1` にフォールバック | テンプレート/`1` |
| `--dist-freeze TEXT` | 調和拘束用の `(i,j,target_Å)` タプルを記述する Python リテラル文字列（繰り返し指定可） | _None_ |
| `--one-based/--zero-based` | `--dist-freeze` インデックスを 1 始まり（デフォルト）または 0 始まりとして解釈 | `True` |
| `--bias-k FLOAT` | すべての `--dist-freeze` タプルに適用される調和バイアス強度（eV·Å⁻²） | `300` |
| `--freeze-links/--no-freeze-links` | キャップ水素の親原子の凍結を切り替え（PDB/mmCIF 入力、または `--ref-pdb` 付き XYZ/GJF） | `True` |
| `--freeze-atoms TEXT` | 凍結する原子の 1 始まりインデックスをカンマ区切りで明示的に指定（例: `'1,3,5'`）。`--freeze-links` と併用可、任意の入力形式に適用 | _None_ |
| `--max-cycles INT` | 最適化反復の上限 | `10000` |
| `--opt-mode TEXT` | 最適化モード: `grad`（`lbfgs`）または `hess`（`rfo`）。`lbfgs`/`rfo` も指定可。サブコマンド別の対応表（`opt` は L-BFGS/RFO、`tsopt` は Dimer/RS-P-RFO）は {ref}`ja-opt-mode-semantics` を参照 | `grad` |
| `--flatten/--no-flatten` | 最適化後の虚振動数モードフラット化ループを有効/無効化 | `False` |
| `--reject-uphill/--no-reject-uphill` | `hess` モードで RFO の上り坂試行ステップ拒否を明示的に有効化（許容値 `1e-4` Hartree、低エネルギー形状へロールバックして trust radius を縮小）。`grad`/`lbfgs` モードでは無効。emergency trust floor 到達時は、非収束停止を報告する前に保持構造を通常の収束条件で最終確認 | `False` |
| `--dump/--no-dump` | 軌跡ダンプ（`optimization_trj.xyz`）を出力 | `False` |
| `--convert-files/--no-convert-files` | PDB/mmCIF topology入力用の XYZ/TRJ → PDB/CIF、および Gaussian template用の XYZ → GJF を切り替え | `True` |
| `--ref-pdb FILE` | XYZ/GJF入力に使用する参照PDBまたはmmCIF topology | _None_ |
| `-o, --out-dir TEXT` | すべてのファイルの出力ディレクトリ | `./result_opt/` |
| `--thresh TEXT` | 収束プリセットの上書き（`gau_loose`、`gau`、`gau_tight`、`gau_vtight`、`baker`、`never`） | `gau` |
| `--config FILE` | ベース YAML 設定ファイル | _None_ |
| `--show-config/--no-show-config` | 実行前に解決済み YAML レイヤ情報を表示 | `False` |
| `--out-json/--no-out-json` | `out_dir` に機械可読な `result.json` を書き出す。スキーマは [JSON 出力スキーマ](json-output.md) を参照 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに設定検証と実行計画表示のみ行う | `False` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP バックエンド | `uma` |
## YAML 設定

共有セクションは [YAML リファレンス](yaml-reference.md) を再利用し、変更が必要な値だけを調整します。`geom`、`calc`、`opt`、オプティマイザ固有の `lbfgs`/`rfo` ブロックは正規のキーとデフォルトを使用します。代表的な最小構成:

```yaml
geom:
  coord_type: cart        # or `dlc` for delocalized internal coordinates
  freeze_atoms: []        # 1-based frozen indices; merged with CLI link detection
calc:
  charge: 0               # mirrors the CLI option; defaults from `.gjf` when present
  spin: 1
opt:
  thresh: gau
  max_cycles: 10000
  out_dir: ./result_opt/  # opt-specific default
```

### `geom`
- `coord_type`（`"cart"`）: デカルト座標 vs `"dlc"` 非局在化内部座標
- `freeze_atoms`（`[]`）: 1 始まりの凍結原子インデックス。CLI のキャップ検出結果と自動的にマージされます

### `calc`
- MLIP バックエンド設定（`model`、`task_name`、デバイス選択、近傍半径、Hessian 形式など）
- `charge`/`spin` は CLI オプションに対応（`.gjf` がある場合はテンプレート値がデフォルト）

### `opt`
L-BFGS と RFO の両方で使用される共有オプティマイザ制御:
- `thresh` プリセット（Gaussian 系または Baker 系）。`baker` は4列
  （max/rms force、max/rms step）に加えて `|delta E| < 1e-6` を要求し、
  文献の Baker 基準より厳しい設定です。プリセット名は
  `pdb2reaction/core/defaults.py`（`THRESH_CHOICES`）に定義されます。
- `max_cycles`、`print_every`（`100`）、`min_step_norm`（`1e-8`）、`assert_min_step`、収束切り替え（`rms_force` など）、RMSD ベースの `converge_to_geom_rms_thresh`、`overachieve_factor`、`check_eigval_structure`、`line_search`。
- 平坦なエネルギー地形による停止（`energy_plateau`、`energy_plateau_thresh`、`energy_plateau_window`）— デフォルトは無効で、`--stop-plateau` で opt-in します。有効時は直近ステップのエネルギーレンジが平坦化した場合、収束扱いにせず `stalled` として停止します（MLIP の力のノイズで力ベース収束に到達できない場合に有効。下の注記を参照）。
- ダンプ/管理項目（`dump`、`dump_restart`、`prefix`、`out_dir`）。

### `lbfgs`
`opt` を L-BFGS 固有の設定で拡張: `keep_last`、`beta`、`gamma_mult`、`max_step`、`control_step`、`double_damp`、およびオプションの正則化パラメータ `mu_reg`/`max_mu_reg_adaptions`

### `rfo`
`opt` を RFOptimizer 固有の設定で拡張: 信頼領域サイジング（`trust_radius`、`trust_min`、`trust_max`、`trust_update`）、`max_energy_incr`、Hessian 管理（`hessian_update`、`hessian_init`、`hessian_recalc`、`hessian_recalc_adapt`、`small_eigval_thresh`）、マイクロイテレーション制御（`alpha0`、`max_micro_cycles`、`rfo_overlaps`）、DIIS ヘルパー（`gdiis`、`gediis`、閾値、`gdiis_test_direction`）、および `adapt_step_func`

### opt 固有のデフォルト

`opt` サブコマンドは `opt.out_dir`（および `lbfgs.out_dir` / `rfo.out_dir`）を `./result_opt/` に設定します。

`geom`、`calc`、`opt`、`lbfgs`、`rfo` の完全な YAML スキーマは [YAML リファレンス](yaml-reference.md) を参照してください。

## 注記

```{note}
**平坦なエネルギー地形による停止（opt-in、デフォルト無効）。** `--stop-plateau`
（YAML `energy_plateau: true`）を指定すると、直近 `--stop-plateau-window`
ステップのエネルギーレンジ（max − min）が `--stop-plateau-thresh`
（デフォルト `1×10⁻⁴ au ≈ 0.06 kcal/mol`、50 ステップ）を
下回った時点で、optimizerは収束扱いにせず `stalled` として停止します。これにより、backend/model/system依存の
force noise/flatnessが選択したforce閾値への到達を妨げる場合でも、無駄なcycleを
消費せずに停止できます。エネルギーの平坦化は停留点の証拠ではないため収束扱いにはならず、
`--max-cycles` が常に実質的な上限です。なお chain-of-states optimizer
（イメージごとのエネルギー配列を保持するもの）ではこの停止はスキップされます。
```

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
- [tsopt](tsopt.md) — 極小ではなく遷移状態（鞍点）を最適化
- [freq](freq.md) — 最適化が極小に達したことを確認する振動解析
- [extract](extract.md) — 最適化前に活性部位モデル（バインディングポケット） PDB を生成
- [all](all.md) — 端点を事前最適化する一気通貫ワークフロー
- [YAML リファレンス](yaml-reference.md) — `opt`、`lbfgs`、`rfo` の完全な設定オプション
- [用語集](glossary.md) — L-BFGS、RFO の定義
