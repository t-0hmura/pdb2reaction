# `scan3d`

調和拘束と機械学習原子間ポテンシャル（MLIP）緩和により、3 距離 `(d₁, d₂, d₃)` のグリッドスキャンを行い、その距離空間上の 3次元ポテンシャルエネルギー分布（マップ）を生成します。この分布を得たいとき、または既存の `surface.csv` を再プロットしたいときに使用します。

コマンドの呼び出し方は 2 通りあります。新規スキャンを実行するには、ターゲットを `-s/--scan-lists` で指定します（YAML/JSON ファイルパスが推奨、またはインライン Python リテラル）。エネルギーを再評価せずに既存の `surface.csv` を再プロットするだけなら、`--csv` で渡します。スキャン中、`scan3d` は d₁ → d₂ → d₃ の順にループをネストし、対応する調和拘束をかけて各格子点を緩和します。

デフォルトのオプティマイザは L-BFGS（`--opt-mode grad`）です。RFOptimizer が必要な場合は `--opt-mode hess` を指定してください。

XYZ/GJF 入力では、`--ref-pdb` で参照 PDB/mmCIF トポロジーを指定すると、XYZ 座標を保持したまま PDB/CIF/GJF companion を生成できます。

## 実行例

コマンド形式:

```bash
pdb2reaction scan3d [-i INPUT.{pdb|xyz|trj|...}] [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m MULT] \
 [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] \
 [-s/--scan-lists scan3d.yaml | '[(i,j,lowÅ,highÅ), (i,j,lowÅ,highÅ), (i,j,lowÅ,highÅ)]'] [options] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE] [--csv PATH]
```

推奨: YAML/JSON spec ファイル。

```bash
# 推奨: YAML/JSON spec
cat > scan3d.yaml << 'YAML'
one_based: true
pairs:
 - ["TYR,285,CA", "SAM,309,C10", 1.30, 3.10]
 - ["TYR,285,CB", "SAM,309,C11", 1.20, 3.20]
 - ["TYR,285,CG", "SAM,309,C12", 1.10, 3.00]
YAML
pdb2reaction scan3d -i input.pdb -q 0 -s scan3d.yaml
```

代替: インライン Python リテラル。

```bash
# 代替: Python リテラル
pdb2reaction scan3d -i input.pdb -q 0 \
 -s '[("TYR,285,CA","SAM,309,C10",1.30,3.10),("TYR,285,CB","SAM,309,C11",1.20,3.20),("TYR,285,CG","SAM,309,C12",1.10,3.00)]'
```

L-BFGS 緩和、内側軌跡ダンプ、HTML 等値面プロット。

```bash
# L-BFGS 緩和、内側軌跡ダンプ、HTML 等値面プロット
pdb2reaction scan3d -i input.pdb -q 0 \
 -s '[("TYR,285,CA","SAM,309,C10",1.30,3.10),("TYR,285,CB","SAM,309,C11",1.20,3.20),("TYR,285,CG","SAM,309,C12",1.10,3.00)]' \
 --max-step-size 0.20 --dump -o ./result_scan3d/ --opt-mode grad \
 --preopt --baseline min
```

既存 `surface.csv` からのプロットのみ（新規エネルギー評価をスキップ）。

```bash
# 既存 surface.csv からのプロットのみ（スキャンしない）
pdb2reaction scan3d --csv ./result_scan3d/surface.csv --zmin -10 --zmax 40 -o ./result_scan3d/
```

## 処理の流れ

1. `geom_loader` で構造を読み込み、CLI または Gaussian テンプレートから電荷とスピンを解決します。`--preopt` の場合は無バイアスの事前最適化を実行します。`-q` が省略され `--ligand-charge` が与えられている場合、構造は酵素--基質複合体として扱われ、PDB/mmCIF 入力（または `--ref-pdb` 付き XYZ/GJF）で `extract.py` の電荷サマリーから総電荷を導出します。
2. `-s/--scan-lists`を3つの4要素tupleへ解析します。3-field selectorは順不同で、重複する残基名／番号には位置固定`CHAIN:RESNAME:RESSEQ[ICODE]:ATOM`を使います。`h = --max-step-size`で各距離の線形gridを生成し、参照構造の距離に近い値から走査します。このため CSV の index は距離の昇順ではなく走査順です。
3. 外側ループで `d1[i]` を走査し、**d₁ 拘束のみ**を適用して緩和します。近い d₁ 値の既存構造から開始します。
4. 中間ループで `d2[j]` を走査し、**d₁ + d₂ 拘束**を適用して緩和します。近い (d₁, d₂) の構造から開始します。
5. 内側ループで `d3[k]` を走査し、**3 つの拘束すべて**を適用して緩和します。バイアスを除去したエネルギーを測定し、構造と収束フラグを書き出します。
6. 完了後に `surface.csv`（カラム: `i,j,k,d1_A,d2_A,d3_A,energy_hartree,bias_converged,is_preopt,energy_kcal,d1_label,d2_label,d3_label`）を組み立て、`--baseline {min|first}` で kcal/mol の基準を設定し、`--zmin/--zmax` に従った 3D RBF 補間等値面図 `scan3d_density.html` を生成します。丸め後の距離タグが別の点と重なる場合、後の構造ファイル名には 0 始まりの格子 index `_grid_III_JJJ_KKK` が付きます。`--csv` が指定された場合、この可視化ステップのみを実行します。プロットには `d1_A`、`d2_A`、`d3_A` と Hartree または kcal/mol の energy 列が必要で、事前最適化行・明示的な非収束行・非有限行を除外します。4 点以上の非共面 usable point が 3 軸すべてを張る必要があります。

## 出力

主要な成果物は `surface.csv`、`grid/` 配下の各点の構造、そして等値面図 `scan3d_density.html` です。

```
out_dir/ (デフォルト:./result_scan3d/)
├─ surface.csv # グリッドメタデータ（i=j=k=-1 の参照行を含む場合あり）
├─ scan3d_density.html # 3D エネルギー等値面の可視化
├─ grid/point_i###_j###_k###.xyz # 各グリッド点の緩和構造（Å×100 タグ）
├─ grid/point_i###_j###_k###.pdb # 変換有効時に対応する PDB
├─ grid/point_i###_j###_k###.cif # bridge入力。元IDを復元
├─ grid/point_i###_j###_k###.gjf # テンプレートがある場合に対応する Gaussian
├─ grid/preopt_i###_j###_k###.xyz # スキャン開始前の構造（--preopt の場合は最適化済み）
└─ grid/inner_path_d1_###_d2_###_trj.xyz # --dump の場合のみ（参照トポロジーと変換があれば PDB/CIF companion も生成）
```

グリッド点の構造は `Å×100` タグを使用するため、`point_i130_j310_k200.xyz` は d₁=1.30, d₂=3.10, d₃=2.00 Å に対応します。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| **入力と電荷** | | |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル（PDB / XYZ / TRJ / GJF） | `--csv` 未指定時は必須 |
| `-q, --charge INT` | 新規 scan の総電荷（CLI > template/`--ligand-charge`）。両方指定時は `-q` が優先。plot-only `--csv` mode では不要 | 新規 scan で template/導出がない場合は必須 |
| `-l, --ligand-charge TEXT` | 新規 scan で、単一整数または残基別 mapping から PDB/mmCIF の総電荷を導出。`-q` 省略時に使用。plot-only `--csv` mode では使用しない | _None_ |
| `-m, --multiplicity INT` | スピン多重度 2S+1。`.gjf` テンプレートがあれば継承し、未指定時は `1` | `.gjf` テンプレート値または `1` |
| **バックエンドと計算** | | |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP バックエンド | `uma` |
| `--workers`, `--workers-per-node` | UMA 予測器の並列度（`workers_per_node` は並列予測器へ転送）。`workers > 1` と明示的な解析 Hessian は併用不可。{ref}`ja-workers-analytical-error` を参照 | `1`, `1` |
| `--solvent TEXT` | xTB 暗黙溶媒（例: `water`）。`none` で無効化 | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB 溶媒モデル | `alpb` |
| **活性領域の凍結** | | |
| `--freeze-links/--no-freeze-links` | PDB/mmCIF トポロジー入力時にキャップ水素の親原子を凍結 | `True` |
| `--freeze-atoms TEXT` | 凍結する原子の 1 始まりインデックスをカンマ区切りで明示的に指定（例: `'1,3,5'`）。`--freeze-links` と併用可、任意の入力形式に適用 | _None_ |
| **スキャンターゲット** | | |
| `-s, --scan-lists TEXT` | スキャンターゲット: YAML/JSON スペックファイルパス（推奨）または **単一**のインライン Python リテラルで 3 つの4 要素タプル `(i,j,lowÅ,highÅ)` を指定。`i`/`j` は整数インデックスまたは PDB セレクタ | `--csv` 未指定時に必須 |
| `--one-based/--zero-based` | `(i, j)` のインデックスを 1 始まり/0 始まりとして解釈 | `True` |
| `--print-parsed/--no-print-parsed` | `-s/--scan-lists` 解釈後のペア情報を表示 | `False` |
| `--max-step-size FLOAT` | 各距離の 1 増分あたりの最大変化量（Å）。グリッド密度を決定 | `0.20` |
| **緩和** | | |
| `--bias-k FLOAT` | 調和バイアス強度 `k`（eV·Å⁻²） | `300` |
| `--opt-mode TEXT` | `grad` → L-BFGS、`hess` → RFOptimizer | `grad` |
| `--relax-max-cycles INT` | 各バイアス緩和の最大最適化サイクル数。明示値は YAML `opt.max_cycles` を上書き | `10000` |
| `--thresh TEXT` | 収束プリセットの上書き（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） | `baker` |
| `--preopt/--no-preopt` | スキャン前に無バイアス最適化を実行 | `False` |
| **マージとアラインメント** | | |
| `--ref-pdb FILE` | XYZ/GJF 入力時の参照 PDB/mmCIF トポロジー（XYZ 座標を保持） | _None_ |
| `--convert-files/--no-convert-files` | 入力トポロジー／テンプレートに応じた XYZ/TRJ → PDB/CIF/GJF companion 生成を切り替え | `True` |
| **出力と設定** | | |
| `-o, --out-dir TEXT` | グリッドとプロットの出力ディレクトリ | `./result_scan3d/` |
| `--csv PATH` | 既存の `surface.csv` を読み込みプロットのみ実行（新規スキャンなし）。指定時は `-i/--input` と `-s/--scan-lists` が任意になります | _None_ |
| `--dump/--no-dump` | 各 (d₁, d₂) ペアの `inner_path_d1_###_d2_###_trj.xyz` を保存 | `False` |
| `--baseline {min,first}` | Hartree 列がある energy の基準を global minimum または最初の eligible grid point に設定。kcal-only の `--csv` input は与えられた zero を保持 | `min` |
| `--zmin FLOAT`, `--zmax FLOAT` | 等値面の色範囲（kcal/mol） | 自動 |
| `--out-json/--no-out-json` | `out_dir` に機械可読な `result.json` を書き出す。スキーマは [JSON 出力スキーマ](json-output.md) を参照 | `False` |
| `--config FILE` | ベース YAML 設定ファイル（最初に適用） | _None_ |

## YAML 設定

### 共有 YAML セクション
- `geom`, `calc`, `opt`, `lbfgs`, `rfo`: [YAML リファレンス](yaml-reference.md) と同じキーを使用しますが、run-scoped の `opt.dump` と optimizer `out_dir` は無視されます。軌跡出力は `--dump`、output directory は command-owned の `-o/--out-dir` で制御します。

```yaml
geom:
 coord_type: cart # coordinate type: cartesian vs dlc internals
 freeze_atoms: [] # 1-based frozen atoms merged with CLI/link detection
calc:
 charge: 0 # total charge (CLI/template override)
 spin: 1 # spin multiplicity 2S+1
 model: uma-s-1p2 # uma-s-1p2 | uma-m-1p1
 device: auto # MLIP device selection
opt:
 thresh: baker # convergence preset (default: baker)
 max_cycles: 10000 # optimizer cycle cap
 dump: false # optimizer dumps (scan trajectories are controlled by --dump)
lbfgs:
 max_step: 0.3 # maximum step length
rfo:
 trust_radius: 0.10 # trust-region radius
bias:
 k: 300.0 # harmonic bias strength (eV·Å⁻²)
```

明示した `--relax-max-cycles` は YAML `opt.max_cycles` を上書きします。省略時は YAML が優先され、いずれもなければデフォルト `10000` です。

### セクション `bias`
- `k`（`300`）: 調和バイアス強度（eV·Å⁻²）。

## 注記

- `scan3d` はちょうど **3 つ**の4 要素タプル `(i, j, low_Å, high_Å)` を受け付けます（YAML/JSON では `pairs` キー、インラインでは単一リテラル）。`scan` と異なり、リテラルは **1 つだけ**を受け付けます（複数ステージは非対応）。YAML/JSON ファイル書式、インライン Python リテラル構文、原子セレクタ、クォート規則については {ref}`CLI 規約: スキャンリスト仕様 <ja-scan-list-spec>` を参照してください。
- 3D グリッドは点数が急激に増加するため、まず `--max-step-size` を大きくするか範囲を狭めることを検討してください。
- 計算エンジンは MLIP バックエンド（デフォルト: UMA）で、1D/2D スキャンと同じ `HarmonicBiasCalculator` を再利用します。
- Å 単位の制限値は内部で Bohr に変換され、L-BFGS ステップや RFO 信頼半径の制御に使われます。最適化の一時ファイルはテンポラリディレクトリに配置されます。
- `--baseline` はデフォルトでグローバル最小値を基準としてゼロにします。`--baseline first` は `(i,j,k)=(0,0,0)` が eligible ならその点を使い、対象外なら eligible minimum へ fallback します。`energy_hartree` がなく `energy_kcal` だけの plot-only CSV では、与えられた zero を保持します。
- 3D 可視化は 50×50×50 グリッドでの RBF 補間と、半透明の段階的等値面を使用します（断面表示はありません）。
- `--freeze-links` はユーザー指定の `freeze_atoms` にキャップ水素親原子をマージし、抽出された活性部位モデルの境界を固定します。
- `-s/--scan-lists` の解釈結果を確認したい場合は `--print-parsed` を追加してください。
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

## 関連項目
- [scan](scan.md) — 1D 結合距離スキャン
- [scan2d](scan2d.md) — 2D 距離グリッドスキャン
- [opt](opt.md) — スキャン前後の単一構造最適化
- [all](all.md) — 一気通貫ワークフロー
- [典型エラー別レシピ](recipes-common-errors.md) — 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) — 詳細な対処ガイド
