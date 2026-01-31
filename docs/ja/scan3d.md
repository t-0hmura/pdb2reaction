# `scan3d`

## 概要
`scan3d` は、UMA 計算機と調和拘束を使った**3距離のグリッドスキャン**を実行します。`--scan-lists` には 3つの四つ組 `(i, j, lowÅ, highÅ)` を含む**1つのリテラル**を与えます。`--max-step-size` で各距離の線形グリッドを作り、（事前最適化された）開始構造に近い値が先に訪れられるよう並べ替えます。ループは外側 d₁ → 中間 d₂ → 内側 d₃ の順です。各グリッド点は対応する拘束をかけて緩和され、**バイアスを除去した**エネルギーが記録されます。`surface.csv` を事前に用意して、再計算せずに可視化だけを行うことも可能です。デフォルトの `--opt-mode` は **light**（LBFGS）です。RFOptimizerを使用する場合は `--opt-mode heavy` を指定してください。
XYZ/GJF入力では、`--ref-pdb` が参照 PDB トポロジーを提供し、XYZ座標は保持されるため、PDB/GJF変換が可能になります。

> 3D図の視認性を調整したい場合は、スキャン完了後にCSVを読み込むモードにして `--zmin` / `--zmax` を変更する運用を推奨します。

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

## ワークフロー
1. `geom_loader` で構造を読み込み、CLIまたはGaussianテンプレートから電荷/スピンを解決し、`--preopt True` の場合は無バイアスの事前最適化を実行します。`-q` が省略され `--ligand-charge` が与えられている場合、構造は酵素–基質複合体として扱われ、PDB 入力（または `--ref-pdb` 付きXYZ/GJF）で `extract.py` の電荷サマリーから総電荷を導出します。
2. 単一の `--scan-lists` リテラル（デフォルト1始まり、`--one-based False` で0始まり）を3つの四つ組に解析します。PDB 入力では、各原子指定は整数インデックスまたは `'TYR,285,CA'` のようなセレクタ文字列を使用できます。区切りは空白/カンマ/スラッシュ/バッククォート/バックスラッシュで、トークン順は任意です。`h = --max-step-size` で各距離の線形グリッドを生成し、開始距離に近い値が先に訪れられるよう並べ替えます。
3. 外側ループで `d1[i]` を回し、**d₁拘束のみ**を適用して緩和します。近い d₁ 値の既存構造から開始します。
4. 中間ループで `d2[j]` を回し、**d₁ + d₂拘束**を適用して緩和します。近い (d₁, d₂) 構造から開始します。
5. 内側ループで `d3[k]` を回し、**3拘束**を適用して緩和します。バイアスを外したエネルギーを測定し、構造と収束フラグを書き出します。
6. 完了後に `surface.csv` を組み立て、`--baseline {min|first}` で kcal/mol の基準を設定し、`--zmin/--zmax` を尊重した 3D RBF 等値面図 `scan3d_density.html` を生成します。`--csv` が指定された場合、この可視化のみを実行します。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル | `--csv` 未指定時は必須 |
| `-q, --charge INT` | 総電荷（CLI > テンプレート/`--ligand-charge`）。両方指定時は `-q` が優先 | テンプレート/導出がない場合は必須 |
| `--ligand-charge TEXT` | `-q` が省略された場合に使う総電荷または残基名ごとのマッピング。PDB 入力（または `--ref-pdb` 付きXYZ/GJF）で電荷導出を有効化 | _None_ |
| `--workers`, `--workers-per-node` | UMA予測器の並列度（workers > 1 で解析ヘシアン無効; `workers_per_node` は並列予測器へ転送） | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度 2S+1 | `1` |
| `--scan-lists TEXT` | **単一**のPythonリテラルで3つの四つ組 `(i,j,lowÅ,highÅ)` を指定。`i`/`j` は整数インデックスまたは PDB セレクタ | `--csv` 未指定時は必須 |
| `--one-based {True\|False}` | `(i, j)` のインデックス解釈 | `True` |
| `--max-step-size FLOAT` | 1距離あたりの最大増分（Å）。グリッド密度を制御 | `0.20` |
| `--bias-k FLOAT` | 調和バイアス強度 `k`（eV·Å⁻²）。`bias.k` を上書き | `100` |
| `--relax-max-cycles INT` | 各バイアス緩和の最大最適化サイクル。`opt.max_cycles` を上書き | `10000` |
| `--opt-mode TEXT` | `light` → LBFGS、`heavy` → RFOptimizer | `light` |
| `--freeze-links {True\|False}` | PDB 入力でリンクHの親を凍結 | `True` |
| `--dump {True\|False}` | `inner_path_d1_###_d2_###.trj` を保存 | `False` |
| `--convert-files {True\|False}` | PDB/Gaussian入力の XYZ/TRJ → PDB/GJF 変換を切り替え | `True` |
| `--ref-pdb FILE` | XYZ/GJF入力時の参照 PDB トポロジー（XYZ座標を保持） | _None_ |
| `--out-dir TEXT` | 出力ディレクトリ | `./result_scan3d/` |
| `--csv PATH` | 既存 `surface.csv` を読み込みプロットのみ実行（新規スキャンなし） | _None_ |
| `--thresh TEXT` | 収束プリセット上書き（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） | `baker` |
| `--args-yaml FILE` | YAML 上書き（`geom`, `calc`, `opt`, `lbfgs`, `rfo`, `bias`） | _None_ |
| `--preopt {True\|False}` | スキャン前に無バイアス最適化を実行 | `True` |
| `--baseline {min,first}` | kcal/mol の基準をグローバル最小 or `(i,j,k)=(0,0,0)` に設定 | `min` |
| `--zmin FLOAT`, `--zmax FLOAT` | 等値面の色範囲（kcal/mol） | 自動 |

### 共有YAMLセクション
- `geom`, `calc`, `opt`, `lbfgs`, `rfo`: [YAML リファレンス](yaml-reference.md) と同じキー。`opt.dump` は `False` に強制され、軌跡制御はCLIで行います。

`opt` の詳細は [docs/opt.md](opt.md) を参照してください。

## YAML 設定（`--args-yaml`）
最小例（詳細は {ref}`opt <yaml-configuration-args-yaml>` を参照）:

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
  dump: false                # trajectory dumping disabled (CLI controls dumping)
  out_dir: ./result_scan3d/  # output directory
lbfgs:
  max_step: 0.3              # maximum step length
  out_dir: ./result_scan3d/  # LBFGS-specific output directory
rfo:
  trust_radius: 0.1          # trust-region radius
  out_dir: ./result_scan3d/  # RFO-specific output directory
bias:
  k: 100.0                  # harmonic bias strength (eV·Å⁻²)
```

`--relax-max-cycles` は明示的に指定された場合のみ `opt.max_cycles` を上書き。未指定の場合はYAMLの `opt.max_cycles` が使われます（デフォルト `10000`）。

### セクション `bias`
- `k`（`100`）: 調和バイアス強度（eV·Å⁻²）。`--bias-k` で上書き。

## 出力
```
out_dir/ (デフォルト: ./result_scan3d/)
├─ surface.csv                     # グリッドメタデータ（i=j=k=-1 の参照行を含む場合あり）
├─ scan3d_density.html             # 3Dエネルギー等値面の可視化
├─ grid/point_i###_j###_k###.xyz   # 各グリッド点の緩和構造（Å×100 タグ）
├─ grid/point_i###_j###_k###.pdb   # 変換有効時のPDBコンパニオン
├─ grid/point_i###_j###_k###.gjf   # テンプレートがある場合のGaussianコンパニオン
├─ grid/preopt_i###_j###_k###.xyz  # スキャン開始前の構造（--preopt True の場合は最適化済み）
└─ grid/inner_path_d1_###_d2_###.trj # --dump True の場合のみ（変換有効時は .pdb/.gjf も生成）
```

## 注意事項
- UMA via `uma_pysis` が唯一の計算機で、1D/2Dスキャンと同じ `HarmonicBiasCalculator` を再利用します。
- Å単位の制限は内部でBohrに変換され、LBFGSステップやRFO信頼半径を制御します。最適化の一時ファイルはテンポラリディレクトリに置かれます。
- `--baseline` はデフォルトでグローバル最小に対してゼロ化し、`--baseline first` は `(i,j,k)=(0,0,0)` の点を基準にします。
- 3D可視化は 50×50×50 グリッドでのRBF補間と、半透明の段階的等値面を使用します（断面表示はしません）。
- `--freeze-links` は `freeze_atoms` とリンクH親原子をマージし、抽出ポケットを固定します。
- 電荷はテンプレートがあればそれを優先。`.gjf` 以外の入力では `-q/--charge` が必須ですが、`--ligand-charge` がある場合は例外（PDB 入力、または `--ref-pdb` 付きXYZ/GJF）。明示的な `-q` が常に優先され、多重度は省略時に `1` です。