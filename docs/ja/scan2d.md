# `scan2d`

調和拘束と MLIP 緩和により、2 距離 `(d₁, d₂)` のグリッドスキャンを行い、`(d₁, d₂)` 上の 2D ポテンシャル面を得ます。TS 領域の特定や、MEP 精密化前の反応ランドスケープの可視化に使います。`scan2d` は `--max-step-size` に基づいて両軸の線形グリッドを作成し、各格子点を対応する拘束付きで緩和して、可視化用にバイアスなしの MLIP エネルギーを記録します。入力は 1 つの構造 + `-s/--scan-lists scan2d.yaml`（推奨）、またはちょうど 2 つの 4 要素タプルを含む `-s/--scan-lists` の **単一** インラインリテラルです。デフォルトのバックエンドは UMA で、`-b/--backend` で他のバックエンドも選択できます。L-BFGS の代わりに RFOptimizer を使う場合は `--opt-mode hess` を指定してください。

XYZ/GJF 入力では、`--ref-pdb` で参照 PDB/mmCIF トポロジーを指定すると、XYZ 座標を保持したまま PDB/CIF/GJF companion を生成できます。

## 実行例

YAML スペックファイルを使った最小実行:

```bash
pdb2reaction scan2d -i input.pdb -q 0 -s scan2d.yaml -o ./result_scan2d/
```

コマンド形式:

```bash
pdb2reaction scan2d -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m MULT] \
 [-b/--backend uma|orb|mace|aimnet2] \
 [-s/--scan-lists scan2d.yaml | '[(i,j,lowÅ,highÅ), (i,j,lowÅ,highÅ)]'] [options] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE]
```

推奨: YAML/JSON スペックファイル:

```bash
cat > scan2d.yaml << 'YAML'
one_based: true
pairs:
 - ["TYR,285,CA", "SAM,309,C10", 1.30, 3.10]
 - ["TYR,285,CB", "SAM,309,C11", 1.20, 3.20]
YAML
pdb2reaction scan2d -i input.pdb -q 0 -s scan2d.yaml
```

代替: インライン Python リテラル:

```bash
pdb2reaction scan2d -i input.pdb -q 0 \
 -s '[("TYR,285,CA","SAM,309,C10",1.30,3.10),("TYR,285,CB","SAM,309,C11",1.20,3.20)]'
```

L-BFGS、内側軌跡ダンプ、Plotly 出力:

```bash
pdb2reaction scan2d -i input.pdb -q 0 \
 -s '[("TYR,285,CA","SAM,309,C10",1.30,3.10),("TYR,285,CB","SAM,309,C11",1.20,3.20)]' \
 --max-step-size 0.20 --dump -o ./result_scan2d/ --opt-mode grad \
 --preopt --baseline min
```

### スキャンリスト仕様

`scan2d` はちょうど **2 つ**の4 要素タプル `(i, j, low_Å, high_Å)` を受け付けます（YAML/JSON では `pairs` キー、インラインでは単一リテラル）。`scan` と異なり、リテラルは **1 つだけ**を受け付けます（複数ステージは非対応）。

YAML/JSON ファイル書式、インライン Python リテラル構文、原子セレクタ、クォート規則については
{ref}`CLI 規約: スキャンリスト仕様 <ja-scan-list-spec>` を参照してください。

## 処理の流れ

1. `geom_loader` で入力構造をロードし、電荷とスピンを解決します。`--preopt` の場合は無バイアスの事前最適化を実行します。`-q` が省略され `--ligand-charge` がある場合、構造は酵素--基質複合体として扱われ、PDB/mmCIF 入力（または `--ref-pdb` 付き XYZ/GJF）では `extract.py` の電荷サマリーから総電荷を導出します。事前最適化構造は `grid/preopt_iDDD_jDDD.*`（`DDD = round(d × 100)` Å）に保存され、`surface.csv` には `i = j = -1` のエントリとしてバイアスなしエネルギーが記録されます。
2. `-s/--scan-lists`を2つの4要素tupleへ解析します。3-field selectorは順不同で、重複する残基名／番号には位置固定`CHAIN:RESNAME:RESSEQ[ICODE]:ATOM`を使います。線形gridは`ceil(|high − low| / h) + 1`点（両端を含む）です。
3. 外側ループで `d1[i]`（近い順）を走査します。各値で **d₁ 拘束のみ**を適用して緩和し、その構造をスナップショットとして保存します。次に内側ループで `d2[j]` を走査し、**d₁ と d₂ の両拘束**を適用して、最も近い既収束構造から緩和を開始します。
4. 各 `(i, j)` について、`<out-dir>/grid/point_iDDD_jDDD.xyz`（`DDD = round(d × 100)` Å。例 `d1=1.30 Å, d2=3.10 Å` → `point_i130_j310.xyz`）に構造を保存し、バイアス収束の可否を記録し、バイアスを除去した MLIP エネルギーを評価します。丸め後のタグが別の点と重なる場合、後のファイル名には 0 始まりの格子 index `_grid_III_JJJ` が付きます。`--dump` の場合、外側ループごとの内側軌跡が `inner_path_d1_###_trj.xyz`（`###` は外側ステップ index）として保存されます。
5. すべての点を走査したら、`i,j,d1_A,d2_A,energy_hartree,bias_converged,is_preopt,energy_kcal,d1_label,d2_label` の列を持つ `<out-dir>/surface.csv` を作成します。任意の参照行は `i = j = -1`、`is_preopt = true` です。`--baseline {min|first}` で kcal/mol の基準をシフトします。`--baseline first` は再並べ替え後の最初の eligible 格子点（`i = j = 0`）を使い、その点が対象外なら eligible point の最小値へ fallback します。eligible point が少なくとも 3 つあり、重複せず非共線の場合だけ、`scan2d_map.png`（2D contour）と `scan2d_landscape.html`（3D surface）を生成します。条件を満たさない場合も `surface.csv` は残り、plot は生成せず正常に扱います。`--zmin/--zmax` でカラースケールを固定できます。

## 出力

実行後は `surface.csv`、`grid/` 配下の各点の構造、`scan2d_map.png` / `scan2d_landscape.html` のプロットを確認してください。

- `result_scan2d/surface.csv`
- `result_scan2d/grid/point_iDDD_jDDD.xyz` (`DDD = round(d × 100)` Å。例: `d1=1.30 Å, d2=3.10 Å` → `point_i130_j310.xyz`)
- `result_scan2d/scan2d_map.png` と `result_scan2d/scan2d_landscape.html`

```text
out_dir/ (デフォルト:./result_scan2d/)
├─ surface.csv # 構造化グリッド表
├─ scan2d_map.png # 2D コンター（Kaleido 必須; PNG 出力に失敗すると実行が停止）
├─ scan2d_landscape.html # 3D サーフェス可視化
├─ grid/point_iDDD_jDDD.xyz # DDD = round(d × 100) Å (例: d1=1.30 Å, d2=3.10 Å → point_i130_j310.xyz)
├─ grid/point_iDDD_jDDD.pdb # 変換有効時に対応する PDB
├─ grid/point_iDDD_jDDD.cif # bridge入力。元IDを復元
├─ grid/point_iDDD_jDDD.gjf # テンプレートがある場合に対応する Gaussian
├─ grid/preopt_iDDD_jDDD.xyz # 事前最適化構造（--preopt が True の場合）、DDD = round(d × 100)
├─ grid/preopt_iDDD_jDDD.pdb # 変換有効時に対応する PDB
├─ grid/preopt_iDDD_jDDD.cif # bridge入力
├─ grid/preopt_iDDD_jDDD.gjf # テンプレートがある場合に対応する Gaussian
└─ grid/inner_path_d1_###_trj.xyz # --dump の場合のみ（参照トポロジーと変換があれば PDB/CIF companion も生成）
```

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル | 必須 |
| `-q, --charge INT` | 総電荷（CLI > テンプレート/`--ligand-charge`）。両方指定時は `-q` が優先 | テンプレート/導出がない場合は必須 |
| `-l, --ligand-charge TEXT` | 単一の整数（例: `-1`）でリガンド総電荷を指定するか、残基別マッピング（例: `GPP:-3,SAM:1`）で PDB/mmCIF 残基電荷から全系の電荷を導出。`-q` 省略時に使用（PDB/mmCIF 入力、または `--ref-pdb` 付き XYZ/GJF） | _None_ |
| `--workers`, `--workers-per-node` | UMA 予測器の並列度（`workers_per_node` は並列予測器へ転送）。`workers > 1` と明示的な解析 Hessian は併用不可。{ref}`ja-workers-analytical-error` を参照 | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度 2S+1。`.gjf` テンプレートがあれば継承し、未指定時は `1` | `.gjf` テンプレート値または `1` |
| `-s, --scan-lists TEXT` | スキャンターゲット: YAML/JSON スペックファイルパス（推奨）または **単一**のインライン Python リテラルで 2 つの4 要素タプル `(i,j,lowÅ,highÅ)` を指定。`i`/`j` は整数インデックスまたは PDB セレクタ（`'TYR,285,CA'`） | 必須 |
| `--one-based/--zero-based` | `(i, j)` のインデックスを 1 始まり/0 始まりとして解釈 | `True` |
| `--print-parsed/--no-print-parsed` | `-s/--scan-lists` 解釈後のペア情報を表示 | `False` |
| `--max-step-size FLOAT` | 各距離の 1 増分あたりの最大変化量（Å）。グリッド密度を決定 | `0.20` |
| `--bias-k FLOAT` | 調和バイアス強度 `k`（eV·Å⁻²） | `300` |
| `--relax-max-cycles INT` | 各バイアス緩和の最大最適化サイクル数。明示値は YAML `opt.max_cycles` を上書き | `10000` |
| `--opt-mode TEXT` | `grad` → L-BFGS、`hess` → RFOptimizer | `grad` |
| `--freeze-links/--no-freeze-links` | PDB/mmCIF トポロジー入力時にキャップ水素の親原子を凍結 | `True` |
| `--freeze-atoms TEXT` | 凍結する原子の 1 始まりインデックスをカンマ区切りで明示的に指定（例: `'1,3,5'`）。`--freeze-links` と併用可、任意の入力形式に適用 | _None_ |
| `--dump/--no-dump` | 外側ループごとの `inner_path_d1_###_trj.xyz` を保存 | `False` |
| `--convert-files/--no-convert-files` | 入力トポロジー／テンプレートに応じた XYZ/TRJ → PDB/CIF/GJF companion 生成を切り替え | `True` |
| `--ref-pdb FILE` | XYZ/GJF 入力時の参照 PDB/mmCIF トポロジー（XYZ 座標を保持） | _None_ |
| `-o, --out-dir TEXT` | 出力ディレクトリ | `./result_scan2d/` |
| `--thresh TEXT` | 収束プリセットの上書き（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） | `baker` |
| `--config FILE` | ベース YAML 設定ファイル（最初に適用） | _None_ |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP バックエンド | `uma` |
| `--preopt/--no-preopt` | スキャン前に無バイアス最適化を実行 | `False` |
| `--baseline {min,first}` | kcal/mol の基準をグローバル最小値または最初の格子点に設定 | `min` |
| `--zmin FLOAT`, `--zmax FLOAT` | カラースケールの下限/上限（kcal/mol） | 自動 |
| `--out-json/--no-out-json` | `out_dir` に機械可読な `result.json` を書き出す。スキーマは [JSON 出力スキーマ](json-output.md) を参照 | `False` |

## YAML 設定

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

### 共有 YAML セクション
- `geom`, `calc`, `opt`, `lbfgs`, `rfo`: [YAML リファレンス](yaml-reference.md) と同じキーを使用しますが、run-scoped の `opt.dump` と optimizer `out_dir` は無視されます。スキャン軌跡は `--dump`、出力先は command-owned の `-o/--out-dir` で制御します。

`opt` の詳細は [YAML リファレンス](yaml-reference.md) を参照してください。

## 注意事項
- 計算エンジンは MLIP バックエンド（デフォルト: UMA、`-b/--backend` で切替可能）で、1D スキャンと同じ `HarmonicBiasCalculator` を再利用します。
- Å 単位の制限値は内部で Bohr に変換され、L-BFGS ステップや RFO 信頼半径の制御に使われます。最適化の一時ファイルはテンポラリディレクトリに配置されます。
- バイアスはエネルギー記録前に除去されるため、`surface.csv` を下流のフィッティングや可視化スクリプトにそのまま利用できます。
- `--freeze-links` はユーザー指定の `freeze_atoms` にキャップ水素親原子をマージし、抽出された活性部位モデルの境界を固定します。
- 明示した `--relax-max-cycles` は YAML `opt.max_cycles` を上書きします。省略時は YAML が優先され、いずれもなければデフォルト `10000` です。

## 関連項目
- [scan](scan.md) -- 1D 結合距離スキャン
- [scan3d](scan3d.md) -- 3D 距離グリッドスキャン
- [opt](opt.md) -- スキャン前後の単一構造最適化
- [all](all.md) -- 一気通貫ワークフロー
- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
