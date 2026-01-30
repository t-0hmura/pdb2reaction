# `opt` サブコマンド

## 概要
`pdb2reaction opt` は、UMAがエネルギー、勾配、ヘシアンを提供しながら、pysisyphus LBFGS（"light"）またはRFOptimizer（"heavy"）エンジンで単一構造の構造最適化を実行します。入力構造は `.pdb`、`.xyz`、`.trj`、または `geom_loader` でサポートされる任意の形式が可能です。設定は**組み込みデフォルト → CLIオーバーライド → `--args-yaml` オーバーライド**の順序で適用され（YAMLが最も優先）、軽量なデフォルトを維持しながら選択的にオプションをオーバーライドできます。オプティマイザープリセットは現在LBFGSベースの**`light`**モードがデフォルトです。

開始構造がPDBまたはGaussianテンプレートの場合、フォーマット対応変換は最適化された構造を `.pdb`（PDB入力）および `.gjf`（Gaussianテンプレート）コンパニオンにミラーリングします（`--convert-files {True\|False}` で制御、デフォルトで有効）。

## 使用法
```bash
pdb2reaction opt -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [-m MULT] \
                 [--opt-mode light|heavy] [--freeze-links {True\|False}] \
                 [--dist-freeze '[(i,j,target_A), ...]'] [--one-based {True\|False}] \
                 [--bias-k K_eV_per_A2] [--dump {True\|False}] [--out-dir DIR] \
                 [--max-cycles N] [--thresh PRESET] [--args-yaml FILE] \
                 [--convert-files {True\|False}] [--ref-pdb FILE]
```

## ワークフロー
- **オプティマイザー**: `--opt-mode light`（デフォルト）→ L-BFGS; `--opt-mode heavy` → 信頼領域制御付きRational Function Optimizer
- **拘束**: `--dist-freeze` はPythonリテラルタプル `(i, j, target_A)` を消費; 3番目の要素を省略すると開始距離を拘束。`--bias-k` はグローバル調和強度（eV·Å⁻²）を設定
- **電荷/スピン解決**: CLI `-q/-m` は `.gjf` テンプレートメタデータをオーバーライドし、それは `calc` デフォルトをオーバーライド

## CLIオプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる入力構造 | 必須 |
| `-q, --charge INT` | 総電荷。`.gjf` テンプレートまたは `--ligand-charge` が提供しない限り必須 | テンプレート/導出が適用されない限り必須 |
| `--ligand-charge TEXT` | `-q` が省略された場合に使用される総電荷または残基名ごとのマッピング | _None_ |
| `-m, --multiplicity INT` | スピン多重度（2S+1）。`.gjf` テンプレートまたは `1` にフォールバック | テンプレート/`1` |
| `--dist-freeze TEXT` | 調和拘束用の `(i,j,target_A)` タプルを記述するPythonリテラルとして解析される文字列 | _None_ |
| `--one-based {True\|False}` | `--dist-freeze` インデックスを1始まり（デフォルト）または0始まりとして解釈 | `True` |
| `--bias-k FLOAT` | すべての `--dist-freeze` タプルに適用される調和バイアス強度（eV·Å⁻²） | `10.0` |
| `--freeze-links {True\|False}` | リンク水素親凍結をトグル（PDB入力のみ） | `True` |
| `--max-cycles INT` | 最適化反復のハードリミット | `10000` |
| `--opt-mode TEXT` | オプティマイザー選択: `light`（LBFGS）または `heavy`（RFO） | `light` |
| `--dump {True\|False}` | 軌跡ダンプ（`optimization.trj`）を出力 | `False` |
| `--convert-files {True\|False}` | PDB入力用のXYZ/TRJ → PDBコンパニオンおよびGaussianテンプレート用のXYZ → GJFコンパニオンを有効/無効化 | `True` |
| `--ref-pdb FILE` | 入力がXYZ/GJFの場合に使用する参照PDBトポロジー | _None_ |
| `--out-dir TEXT` | すべてのファイルの出力ディレクトリ | `./result_opt/` |
| `--thresh TEXT` | 収束プリセットのオーバーライド（`gau_loose`、`gau`、`gau_tight`、`gau_vtight`、`baker`、`never`） | `gau` |
| `--args-yaml FILE` | YAMLオーバーライドを提供（セクション `geom`、`calc`、`opt`、`lbfgs`、`rfo`） | _None_ |

## 出力
```
out_dir/
├─ final_geometry.xyz          # 常に書き込み
├─ final_geometry.pdb          # 入力がPDBで変換が有効な場合のみ
├─ final_geometry.gjf          # Gaussianテンプレートが検出され変換が有効な場合
├─ optimization.trj            # ダンプが有効な場合のみ
├─ optimization.pdb            # 軌跡のPDB変換（PDB入力、変換有効時）
└─ restart*.yml                # opt.dump_restartが設定されている場合のオプションのリスタート
```

## YAML設定（`--args-yaml`）
YAML値はCLIをオーバーライドし、CLIはデフォルトをオーバーライドします。

### `geom`
- `coord_type`（`"cart"`）: デカルト座標 vs `"dlc"` 非局在化内部座標
- `freeze_atoms`（`[]`）: 0始まりの凍結インデックス; CLIリンク検出と自動マージ

### `calc`
- UMA設定（`model`、`task_name`、デバイス選択、近傍半径、ヘシアン形式など）
- `charge`/`spin` はCLIオプションをミラー

### `opt`
LBFGSとRFOの両方で使用される共有オプティマイザー制御:
- `thresh` プリセット、`max_cycles`、`print_every`、`min_step_norm`、収束トグルなど

### `lbfgs`
L-BFGS固有で `opt` を拡張: `keep_last`、`beta`、`gamma_mult`、`max_step`、`control_step`、`double_damp`

### `rfo`
RFOptimizerフィールドで `opt` を拡張: 信頼領域サイジング、ヘシアン管理、マイクロイテレーション制御、DIISヘルパー
