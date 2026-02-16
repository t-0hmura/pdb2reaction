# `scan3d`

## 概要

> **要約:** 調和拘束と UMA 緩和により、3 距離（d₁, d₂, d₃）のグリッドスキャンを行います。`--scan-lists` に 3 つの四つ組 `(i, j, lowÅ, highÅ)` を含む **1 つのリテラル**を与えるか、`--csv` で既存の `surface.csv` を可視化のみ実行できます。

### 要点
- **入力:** 1 つの構造 + `--scan-lists` の **単一**リテラル（四つ組は 3 つ）。`--csv` 指定時はプロットのみで実行可能。
- **訪問順:** 事前最適化構造に近い値が先に走査されるよう、各軸が並べ替えられます。
- **エネルギー:** 記録されるエネルギーは **バイアスを除去して**評価されるため、格子点間で直接比較できます。
- **主な出力:** `surface.csv`、`grid/` 配下の各点の構造、HTML の等値面図（`scan3d_density.html`）。
- **注意:** 3D グリッドは点数が急増します。まず `--max-step-size` を大きくするか範囲を狭めることを検討してください。

`scan3d` は d₁ → d₂ → d₃ の順にループをネストし、対応する拘束をかけて各格子点を緩和します。デフォルトは LBFGS（`--opt-mode light`）で、RFOptimizer が必要な場合は `--opt-mode heavy` を指定してください。

XYZ/GJF 入力では、`--ref-pdb` で参照 PDB トポロジーを指定すると、XYZ 座標を保持したまま PDB/GJF へのフォーマット対応変換が可能になります。

## 使用法
```bash
pdb2reaction scan3d [-i INPUT.{pdb|xyz|trj|...}] [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [-m MULT] \
                    [--scan-lists '[(i,j,lowÅ,highÅ), (i,j,lowÅ,highÅ), (i,j,lowÅ,highÅ)]'] [options] \
                    [--convert-files {True\|False}] [--ref-pdb FILE] [--csv PATH]
```
注: `-i/--input` と `--scan-lists` は `--csv` が指定されていない限り必須です。

### 例
```bash
# 3距離の最小スキャン
pdb2reaction scan3d -i input.pdb -q 0 \
    --scan-lists '[("TYR,285,CA","MMT,309,C10",1.30,3.10),("TYR,285,CB","MMT,309,C11",1.20,3.20),("TYR,285,CG","MMT,309,C12",1.10,3.00)]'

# LBFGS、内側軌跡ダンプ、HTML等値面
pdb2reaction scan3d -i input.pdb -q 0 \
    --scan-lists '[("TYR,285,CA","MMT,309,C10",1.30,3.10),("TYR,285,CB","MMT,309,C11",1.20,3.20),("TYR,285,CG","MMT,309,C12",1.10,3.00)]' \
    --max-step-size 0.20 --dump True --out-dir ./result_scan3d/ --opt-mode light \
    --preopt True --baseline min

# 既存surface.csvからのプロットのみ（スキャンしない）
pdb2reaction scan3d --csv ./result_scan3d/surface.csv --zmin -10 --zmax 40 --out-dir ./result_scan3d/
```

## `--scan-lists` の書式

`--scan-lists` は **1 つの Python リテラル文字列**として CLI に渡します。シェルのクォート処理に注意が必要です。

### 基本構造

リテラルは、ちょうど **3 つ**の四つ組 `(原子1, 原子2, 下限Å, 上限Å)` を含む Python リストです。

```
--scan-lists '[(原子1, 原子2, 下限Å, 上限Å), (原子3, 原子4, 下限Å, 上限Å), (原子5, 原子6, 下限Å, 上限Å)]'
```

- リスト全体を **シングルクォート** `'...'` で囲みます（シェルがカッコや空白を解釈しないようにするため）。
- 各四つ組は 1 つのスキャン軸を定義し、`原子1`–`原子2` 間の距離を `下限Å` から `上限Å` まで走査します。
- `scan` と異なり、リテラルは **1 つだけ**を受け付けます（複数ステージは非対応）。

### 原子の指定方法

原子は**整数インデックス**または **PDB セレクタ文字列**で指定します。

| 方法 | 例 | 備考 |
| --- | --- | --- |
| 整数インデックス | `(1, 5, 1.30, 3.10)` | デフォルトは 1 始まり（`--one-based True`） |
| PDB セレクタ | `("TYR,285,CA", "MMT,309,C10", 1.30, 3.10)` | 残基名、残基番号、原子名の三つ組 |

PDB セレクタのトークンは、カンマ `,`、スペース、スラッシュ `/`、バッククォート `` ` ``、バックスラッシュ `\` のいずれでも区切れます。トークンの順序も任意です。

```bash
# 以下はすべて同じ原子を指定:
"TYR,285,CA"
"TYR 285 CA"
"TYR/285/CA"
"285,TYR,CA"   # 順序は自由
```

### クォートの規則

```bash
# 正しい: シングルクォートでリストを囲み、セレクタはダブルクォート
--scan-lists '[("TYR,285,CA","MMT,309,C10",1.30,3.10),("TYR,285,CB","MMT,309,C11",1.20,3.20),("TYR,285,CG","MMT,309,C12",1.10,3.00)]'

# 正しい: 整数インデックスなら内側のダブルクォートは不要
--scan-lists '[(1, 5, 1.30, 3.10), (2, 8, 1.20, 3.20), (3, 12, 1.10, 3.00)]'

# 非推奨: ダブルクォートで外側を囲むとエスケープが必要
--scan-lists "[(\"TYR,285,CA\",\"MMT,309,C10\",1.30,3.10), ...]"
```

## ワークフロー
1. `geom_loader` で構造を読み込み、CLI または Gaussian テンプレートから電荷とスピンを解決します。`--preopt True` の場合は無バイアスの事前最適化を実行します。`-q` が省略され `--ligand-charge` が与えられている場合、構造は酵素--基質複合体として扱われ、PDB 入力（または `--ref-pdb` 付き XYZ/GJF）で `extract.py` の電荷サマリーから総電荷を導出します。
2. 単一の `--scan-lists` リテラル（デフォルト 1 始まり、`--one-based False` で 0 始まり）を 3 つの四つ組に解析します。PDB 入力では、各原子指定に整数インデックスまたは `'TYR,285,CA'` のようなセレクタ文字列を使用できます。区切りは空白・カンマ・スラッシュ・バッククォート・バックスラッシュのいずれも可で、トークン順序は任意です。`h = --max-step-size` で各距離の線形グリッドを生成し、開始距離に近い値が先に走査されるよう並べ替えます。
3. 外側ループで `d1[i]` を走査し、**d₁ 拘束のみ**を適用して緩和します。近い d₁ 値の既存構造から開始します。
4. 中間ループで `d2[j]` を走査し、**d₁ + d₂ 拘束**を適用して緩和します。近い (d₁, d₂) の構造から開始します。
5. 内側ループで `d3[k]` を走査し、**3 つの拘束すべて**を適用して緩和します。バイアスを除去したエネルギーを測定し、構造と収束フラグを書き出します。
6. 完了後に `surface.csv` を組み立て、`--baseline {min|first}` で kcal/mol の基準を設定し、`--zmin/--zmax` に従った 3D RBF 補間等値面図 `scan3d_density.html` を生成します。`--csv` が指定された場合、この可視化ステップのみを実行します。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル | `--csv` 未指定時は必須 |
| `-q, --charge INT` | 総電荷（CLI > テンプレート/`--ligand-charge`）。両方指定時は `-q` が優先 | テンプレート/導出がない場合は必須 |
| `--ligand-charge TEXT` | `-q` 省略時に使う総電荷または残基名ごとのマッピング。PDB 入力（または `--ref-pdb` 付き XYZ/GJF）で電荷導出を有効化 | _None_ |
| `--workers`, `--workers-per-node` | UMA 予測器の並列度（workers > 1 で解析ヘシアン無効; `workers_per_node` は並列予測器へ転送） | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度 2S+1。`.gjf` テンプレートがあれば継承し、未指定時は `1` | `.gjf` テンプレート値または `1` |
| `--scan-lists, --scan-list TEXT` | **単一**の Python リテラルで 3 つの四つ組 `(i,j,lowÅ,highÅ)` を指定。`i`/`j` は整数インデックスまたは PDB セレクタ | `--csv` 未指定時は必須 |
| `--one-based {True\|False}` | `(i, j)` のインデックスを 1 始まり/0 始まりとして解釈 | `True` |
| `--max-step-size FLOAT` | 各距離の 1 増分あたりの最大変化量（Å）。グリッド密度を決定 | `0.20` |
| `--bias-k FLOAT` | 調和バイアス強度 `k`（eV·Å⁻²） | `300` |
| `--relax-max-cycles INT` | 各バイアス緩和の最大最適化サイクル数。YAML で `opt.max_cycles` が指定されていない場合に使用 | `10000` |
| `--opt-mode TEXT` | `light` → LBFGS、`heavy` → RFOptimizer | `light` |
| `--freeze-links {True\|False}` | PDB 入力時にリンク水素の親原子を凍結 | `True` |
| `--dump {True\|False}` | 各 (d₁, d₂) ペアの `inner_path_d1_###_d2_###.trj` を保存 | `False` |
| `--convert-files {True\|False}` | PDB/Gaussian 入力で XYZ/TRJ → PDB/GJF 変換を切り替え | `True` |
| `--ref-pdb FILE` | XYZ/GJF 入力時の参照 PDB トポロジー（XYZ 座標を保持） | _None_ |
| `--out-dir TEXT` | 出力ディレクトリ | `./result_scan3d/` |
| `--csv PATH` | 既存の `surface.csv` を読み込みプロットのみ実行（新規スキャンなし） | _None_ |
| `--thresh TEXT` | 収束プリセットの上書き（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） | `baker` |
| `--args-yaml FILE` | YAML による上書き（`geom`, `calc`, `opt`, `lbfgs`, `rfo`, `bias`） | _None_ |
| `--preopt {True\|False}` | スキャン前に無バイアス最適化を実行 | `True` |
| `--baseline {min,first}` | kcal/mol の基準をグローバル最小値または `(i,j,k)=(0,0,0)` に設定 | `min` |
| `--zmin FLOAT`, `--zmax FLOAT` | 等値面の色範囲（kcal/mol） | 自動 |

### 共有 YAML セクション
- `geom`, `calc`, `opt`, `lbfgs`, `rfo`: [YAML リファレンス](yaml-reference.md) と同じキーを使用します。`opt.dump` は YAML で設定可能ですが、軌跡出力は `--dump` で制御します。

`opt` の詳細は [docs/opt.md](opt.md) を参照してください。

## YAML 設定（`--args-yaml`）
最小例（詳細は {ref}`opt <ja-yaml-configuration-args-yaml>` を参照）:

```yaml
geom:
  coord_type: cart           # coordinate type: cartesian vs dlc internals
  freeze_atoms: []           # 0-based frozen atoms merged with CLI/link detection
calc:
  charge: 0                  # total charge (CLI/template override)
  spin: 1                    # spin multiplicity 2S+1
  model: uma-s-1p1           # UMA model tag
  device: auto               # UMA device selection
opt:
  thresh: baker              # convergence preset (default: baker)
  max_cycles: 10000          # optimizer cycle cap
  dump: false                # optimizer dumps (scan trajectories are controlled by --dump)
  out_dir: ./result_scan3d/  # output directory
lbfgs:
  max_step: 0.3              # maximum step length
  out_dir: ./result_scan3d/  # LBFGS-specific output directory
rfo:
  trust_radius: 0.1          # trust-region radius
  out_dir: ./result_scan3d/  # RFO-specific output directory
bias:
  k: 300.0                  # harmonic bias strength (eV·Å⁻²)
```

`--relax-max-cycles` は**明示的に指定され**、かつ YAML で `opt.max_cycles` が設定されていない場合にのみ適用されます（デフォルト `10000`）。

### セクション `bias`
- `k`（`300`）: 調和バイアス強度（eV·Å⁻²）。

## 出力
```
out_dir/ (デフォルト: ./result_scan3d/)
├─ surface.csv                     # グリッドメタデータ（i=j=k=-1 の参照行を含む場合あり）
├─ scan3d_density.html             # 3D エネルギー等値面の可視化
├─ grid/point_i###_j###_k###.xyz   # 各グリッド点の緩和構造（Å×100 タグ）
├─ grid/point_i###_j###_k###.pdb   # 変換有効時の PDB コンパニオン
├─ grid/point_i###_j###_k###.gjf   # テンプレートがある場合の Gaussian コンパニオン
├─ grid/preopt_i###_j###_k###.xyz  # スキャン開始前の構造（--preopt True の場合は最適化済み）
└─ grid/inner_path_d1_###_d2_###.trj # --dump True の場合のみ（変換有効時は .pdb/.gjf も生成）
```

## 注意事項
- `uma_pysis` 経由の UMA が唯一の計算バックエンドであり、1D/2D スキャンと同じ `HarmonicBiasCalculator` を再利用します。
- Å 単位の制限値は内部で Bohr に変換され、LBFGS ステップや RFO 信頼半径の制御に使われます。最適化の一時ファイルはテンポラリディレクトリに配置されます。
- `--baseline` はデフォルトでグローバル最小値を基準としてゼロにします。`--baseline first` は `(i,j,k)=(0,0,0)` の格子点を基準にします。
- 3D 可視化は 50x50x50 グリッドでの RBF 補間と、半透明の段階的等値面を使用します（断面表示はありません）。
- `--freeze-links` はユーザー指定の `freeze_atoms` にリンク水素親原子をマージし、抽出ポケットを固定します。
- 電荷は Gaussian テンプレートがあれば継承されます。`.gjf` 以外の入力では `-q/--charge` が必須ですが、`--ligand-charge` がある場合は例外です（PDB 入力、または `--ref-pdb` 付き XYZ/GJF）。明示的な `-q` は常に最優先されます。**多重度は `.gjf` テンプレートがあれば継承され、未指定時は `1` になります。**
