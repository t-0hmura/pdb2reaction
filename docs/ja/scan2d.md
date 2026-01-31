# `scan2d`

## 概要
`scan2d` は、2つの距離 (d₁, d₂) のグリッドスキャンを調和拘束とUMA最適化で実行します。1つの `--scan-lists` リテラルに2つの四つ組 `(i, j, lowÅ, highÅ)` を与えると、`--max-step-size` に基づいて両軸の線形グリッドが作成されます。その後、**（事前最適化された）構造に最も近い点を先に訪れる**ように各軸が並べ替えられます。各グリッド点は緩和され、プロット可能なCSV/図出力が生成されます。`surface.csv` のエネルギーは常に**バイアスなし**で評価されるため、グリッド点を直接比較できます。最適化は `--opt-mode light`（デフォルト）の LBFGS、または `--opt-mode heavy` の RFOptimizer を使用します。
XYZ/GJF入力では、`--ref-pdb` が参照 PDB トポロジーを提供し、XYZ座標は保持されるため、PDB/GJF変換が可能になります。

## 使用法
```bash
pdb2reaction scan2d -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [-m MULT] \
                    --scan-lists '[(i,j,lowÅ,highÅ), (i,j,lowÅ,highÅ)]' [options]
                    [--convert-files {True\|False}] [--ref-pdb FILE]
```

### 例
```bash
# 2距離の最小スキャン
pdb2reaction scan2d -i input.pdb -q 0 \
    --scan-lists '[("TYR,285,CA","MMT,309,C10",1.30,3.10),("TYR,285,CB","MMT,309,C11",1.20,3.20)]'

# LBFGS + 内側軌跡ダンプ + Plotly出力
pdb2reaction scan2d -i input.pdb -q 0 \
    --scan-lists '[("TYR,285,CA","MMT,309,C10",1.30,3.10),("TYR,285,CB","MMT,309,C11",1.20,3.20)]' \
    --max-step-size 0.20 --dump True --out-dir ./result_scan2d/ --opt-mode light \
    --preopt True --baseline min
```

## ワークフロー
1. `geom_loader` で入力構造をロードし、電荷/スピンを解決し、`--preopt True` の場合は無バイアスの事前最適化を実行します。`-q` が省略され `--ligand-charge` がある場合、構造は酵素–基質複合体として扱われ、PDB 入力（または `--ref-pdb` 付きXYZ/GJF）では `extract.py` の電荷サマリーから総電荷を導出します。事前最適化構造は `grid/preopt_i###_j###.*` に保存され、`surface.csv` には `i = j = -1` で無バイアスエネルギーが記録されます。
2. 単一の `--scan-lists` リテラルを2つの四つ組に解析し、インデックスを正規化（デフォルトは1始まり）。PDB 入力では、各エントリは整数インデックスまたは `'TYR,285,CA'` のようなセレクタ文字列を指定できます。区切りは空白/カンマ/スラッシュ/バッククォート/バックスラッシュで、トークン順序は任意（フォールバックは resname, resseq, atom）。線形グリッドは `N = ceil(|high − low| / h)`（`h = --max-step-size`）で構成し、長さ0の範囲は1点に縮退します。その後、各軸は事前最適化構造に最も近い距離が `i = 0` / `j = 0` になるよう並べ替えます。
3. `d1[i]`（近い順）ごとに外側ループを回します。各値で **d₁拘束のみ**を適用して緩和し、その構造をスナップショットします。
4. 内側ループで `d2[j]` を回し、**d₁ と d₂ の両拘束**を適用して、最も近い既収束構造から開始して緩和します。
5. 各 `(i, j)` で、`<out-dir>/grid/point_i###_j###.xyz` に構造を保存し、バイアス収束の可否を記録し、バイアスを外したUMAエネルギーを評価します。`--dump True` の場合、外側ループごとの内側軌跡が `inner_path_d1_###.trj` に保存されます。
6. すべての点を訪れたら、`<out-dir>/surface.csv` を作成し、`--baseline {min|first}` で kcal/mol の基準を設定します。`--baseline first` では、基準は `(low₁, low₂)` ではなく再並べ替え後の最初の格子点（`i = j = 0`）になります。`scan2d_map.png`（2Dコンター）と `scan2d_landscape.html`（3Dサーフェス）を `<out-dir>/` に生成します。`--zmin/--zmax` でカラースケールを固定できます。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル | 必須 |
| `-q, --charge INT` | 総電荷（CLI > テンプレート/`--ligand-charge`）。両方指定時は `-q` が優先 | テンプレート/導出がない場合は必須 |
| `--ligand-charge TEXT` | `-q` が省略された場合に使う総電荷または残基名ごとのマッピング。PDB 入力（または `--ref-pdb` 付きXYZ/GJF）で電荷導出を有効化 | _None_ |
| `--workers`, `--workers-per-node` | UMA予測器の並列度（workers > 1 で解析ヘシアン無効; `workers_per_node` は並列予測器へ転送） | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度 2S+1（CLI > テンプレート > 1） | `.gjf` テンプレート値または `1` |
| `--scan-lists TEXT` | **単一**のPythonリテラルで2つの四つ組 `(i,j,lowÅ,highÅ)` を指定。`i`/`j` は整数インデックスまたは PDB セレクタ（`'TYR,285,CA'`） | 必須 |
| `--one-based {True\|False}` | `(i, j)` のインデックス解釈 | `True` |
| `--max-step-size FLOAT` | 1距離あたりの最大増分（Å）。グリッド密度を制御 | `0.20` |
| `--bias-k FLOAT` | 調和バイアス強度 `k`（eV·Å⁻²）。`bias.k` を上書き | `100` |
| `--relax-max-cycles INT` | 各バイアス緩和の最大最適化サイクル。`opt.max_cycles` を上書き | `10000` |
| `--opt-mode TEXT` | `light` → LBFGS、`heavy` → RFOptimizer | `light` |
| `--freeze-links {True\|False}` | PDB 入力でリンクHの親を凍結 | `True` |
| `--dump {True\|False}` | 外側ループごとの `inner_path_d1_###.trj` を保存 | `False` |
| `--convert-files {True\|False}` | PDB/Gaussian入力の XYZ/TRJ → PDB/GJF 変換を切り替え | `True` |
| `--ref-pdb FILE` | XYZ/GJF入力時の参照 PDB トポロジー（XYZ座標を保持） | _None_ |
| `--out-dir TEXT` | 出力ディレクトリ | `./result_scan2d/` |
| `--thresh TEXT` | 収束プリセット上書き（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） | `baker` |
| `--args-yaml FILE` | YAML 上書き（`geom`, `calc`, `opt`, `lbfgs`, `rfo`, `bias`） | _None_ |
| `--preopt {True\|False}` | スキャン前に無バイアス最適化を実行 | `True` |
| `--baseline {min,first}` | kcal/mol の基準をグローバル最小 or 最初の格子点に設定 | `min` |
| `--zmin FLOAT`, `--zmax FLOAT` | カラースケール下限/上限（kcal/mol） | 自動 |

### 共有YAMLセクション
- `geom`, `calc`, `opt`, `lbfgs`, `rfo`: [YAML リファレンス](yaml-reference.md) と同じキー。`opt.dump` は `False` に強制され、軌跡制御はCLIで行います。

### セクション `bias`
- `k`（`100`）: 調和バイアス強度（eV·Å⁻²）。`--bias-k` で上書き。

## 出力
```
out_dir/ (デフォルト: ./result_scan2d/)
├─ surface.csv                # 構造化グリッド表
├─ scan2d_map.png             # 2Dコンター（Kaleido必須; PNG出力に失敗すると停止）
├─ scan2d_landscape.html      # 3Dサーフェス可視化
├─ grid/point_i###_j###.xyz   # 各(i, j)の緩和構造
├─ grid/point_i###_j###.pdb   # 変換有効時のPDBコンパニオン
├─ grid/point_i###_j###.gjf   # テンプレートがある場合のGaussianコンパニオン
└─ grid/inner_path_d1_###.trj # --dump True の場合のみ（PDB入力時は .pdb にも変換）
```

## 注意事項
- UMA via `uma_pysis` が唯一の計算機で、1Dスキャンと同じ `HarmonicBiasCalculator` を再利用します。
- Å単位の制限は内部でBohrに変換され、LBFGSステップやRFO信頼半径を制御します。最適化の一時ファイルはテンポラリディレクトリに置かれます。
- バイアスは最終エネルギー評価前に必ず除去されるため、`surface.csv` を下流のフィッティングや可視化に再利用できます。
- `--freeze-links` は `freeze_atoms` とリンクH親原子をマージし、抽出ポケットを固定します。
- 電荷はテンプレートがあればそれを優先。`.gjf` 以外の入力では `-q/--charge` が必須ですが、`--ligand-charge` がある場合は例外（PDB 入力、または `--ref-pdb` 付きXYZ/GJF）。明示的な `-q` が常に優先され、多重度は省略時に `1` です。

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
  out_dir: ./result_scan2d/  # output directory
lbfgs:
  max_step: 0.3              # maximum step length
  out_dir: ./result_scan2d/  # LBFGS-specific output directory
rfo:
  trust_radius: 0.1          # trust-region radius
  out_dir: ./result_scan2d/  # RFO-specific output directory
bias:
  k: 100.0                  # harmonic bias strength (eV·Å⁻²)
```

`opt` の詳細は [docs/opt.md](opt.md) を参照してください。
`--relax-max-cycles` は明示的に指定された場合のみ `opt.max_cycles` を上書き。未指定の場合はYAMLの `opt.max_cycles` が使われます（デフォルト `10000`）。