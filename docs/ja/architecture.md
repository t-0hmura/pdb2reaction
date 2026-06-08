# アーキテクチャ: pdb2reaction

---

## 1. 概要

`pdb2reaction` は、活性部位クラスターモデルに対して **純 MLIP による酵素反応経路解析** を実行する Python 製 CLI です。PDB と基質名を起点に、活性部位クラスターを切り出し、切断された結合をリンク水素でキャップし、MLIP ポテンシャル上で RS-I-RFO TS 最適化による Hessian ベースの TS 探索を行い、反応経路を生成します（extract → MEP → tsopt → IRC → freq → dft）。


同梱された 2 つの fork（`pysisyphus/`、`thermoanalysis/`）は、リポジトリ最上位に repo 内部モジュールとして配置されています。これらは意図的に上流の PyPI 配布版では **ありません**。本パッケージと並べて PyPI から再インストールすると、ローカルの拡張が静かに壊れます。§6 を参照してください。

---

## 2. 階層構造（6 つの物理ディレクトリ）

### 2.1 階層テーブル

| layer | dir | responsibility | may depend on |
|---|---|---|---|
| **L1 Interface** | `pdb2reaction/cli/` | Click root group、共有 option-decorator ファクトリ（`common_options.py`）、`--help-advanced`、bool flag 正規化、サブコマンドリゾルバ | `workflows/`, `core/` |
| **L2 Application** | `pdb2reaction/workflows/` | サブコマンドごとのオーケストレーション; ステージランナー 1 ファイルにつき 1 つ（`all.py`, `path_search.py`, `tsopt.py`, `extract.py`, `irc.py`, `freq.py`, `dft.py`, …） | `domain/`, `backends/`, `io/`, `core/` |
| **L3 Domain** | `pdb2reaction/domain/` | 化学的に意味を持つヘルパーロジック（結合変化検出、結合サマリ、元素情報伝播） | `core/` |
| **L4a Infra (MLIP)** | `pdb2reaction/backends/` | MLIP バックエンドディスパッチャ + バックエンドごとのアダプタ（UMA / Orb / MACE / AIMNet2）+ xTB ALPB デルタ補正 | `core/` |
| **L4b Infra (I/O)** | `pdb2reaction/io/` | 出力レイアウト、サマリ、軌跡、PDB 修正、エネルギーダイアグラム、Hessian キャッシュ | `core/` |
| **L5 Foundation** | `pdb2reaction/core/` | defaults（単一真実源）、utils（PDB / XYZ / plot ヘルパー）、将来の `errors.py` / `types.py` / `logging.py` | (none) |
| (bundle, not a layer) | `<repo>/pysisyphus/`, `<repo>/thermoanalysis/` | repo 内部 fork（optimiser / thermochemistry） | (sibling, layer-external) |

**依存方向（一方向）**: `L1 → L2 → {L3, L4} → L5`。この方向性ルールは CI のマーカーカバレッジ（`.github/scripts/check_engineering_markers.py`）で強制されます。同梱 fork は階層グラフの外側に位置し、絶対パッケージパス（`from pysisyphus.X import Y`）を介してどの階層からでも import できます。

### 2.2 パッケージツリーの ASCII マップ

```
pdb2reaction/ [GH: t-0hmura/pdb2reaction]
├── pyproject.toml packages.find = ["pdb2reaction*", ...] (glob, frozen)
├── README.md / CONTRIBUTING.md / CHANGELOG.md
├── docs/
│ ├── architecture.md ← this file
│ └──... (Sphinx site, unchanged)
├── pdb2reaction/ ← package body, 6-layer physical dir
│ ├── __init__.py PEP 562 lazy: _LAZY_SYMBOLS / _LAZY_MODULES + __getattr__
│ ├── __main__.py `from pdb2reaction.cli.app import cli`
│ ├── _version.py / py.typed
│ │
│ ├── cli/ # === L1 Interface ===
│ │ ├── app.py Click group + _LAZY_SUBCOMMANDS registry (absolute paths)
│ │ ├── common_options.py @add_print_every_option / @add_irc_pos_def_option / @add_precision_option / @add_coord_type_option / @add_ml_charge_spin_options
│ │ ├── decorators.py resolve_yaml_sources / load_merged_yaml_cfg / _write_error_json
│ │ ├── help_pages.py --help-advanced pager
│ │ ├── bool_compat.py --flag / --no-flag normalisation
│ │ └── default_group.py subcommand resolver, lazy module import
│ │
│ ├── workflows/ # === L2 Application ===
│ │ ├── all.py full pipeline orchestrator (extract → … → DFT)
│ │ ├── path_search.py / path_opt.py MEP search / COS wrapper
│ │ ├── tsopt.py / freq.py / irc.py / dft.py per-stage runners
│ │ ├── opt.py / sp.py / scan.py / scan2d.py /
│ │ │ scan3d.py / scan_common.py geometry opt / single point / scans
│ │ ├── extract.py active-site extraction CLI
│ │ ├── restraints.py restraint helpers
│ │ └── align_freeze.py Kabsch + frozen-subset rmsd
│ │
│ ├── domain/ # === L3 Domain ===
│ │ ├── bond_changes.py R↔P bond detection
│ │ ├── bond_summary.py post-IRC diagnostic
│ │ └── add_elem_info.py PDB element column normaliser
│ │
│ ├── backends/ # === L4a Infra (MLIP) ===
│ │ ├── __init__.py backend dispatch + registry
│ │ ├── base.py MLIPCalculator protocol
│ │ ├── uma.py / orb.py / mace.py / aimnet2.py per-backend adapters
│ │ ├── solvent.py xTB ALPB implicit-solvent helper
│ │ └── xtb_alpb_correction.py xTB ALPB delta correction
│ │
│ ├── io/ # === L4b Infra (I/O) ===
│ │ ├── summary.py summary.json / summary.md writer
│ │ ├── energy_diagram.py Plotly diagram
│ │ ├── trj2fig.py trajectory → PNG / SVG / PDF / HTML / CSV
│ │ ├── pdb_fix.py altloc resolution
│ │ └── hessian_cache.py in-memory Hessian cache
│ │
│ └── core/ # === L5 Foundation ===
│   ├── defaults.py C1 single source of truth for every default
│   └── utils.py PDB / XYZ / plot helpers
│
├── tests/ smoke / unit
├── .github/ workflows/ + scripts/ (docs-quality lint helpers; CI-only)
└── (repo-top sibling, layer-external bundled forks)
 pysisyphus/ ~90 file, repo-internal fork (slimmed; CLI driver + QM backends + wavefunction + dead optimisers / IRC / NEB variants removed)
 thermoanalysis/ 5 file, repo-internal fork
```

### 2.3 階層ごとの責務詳細

**L1 `cli/`**（約 6 ファイル）。Click コマンドを構築し argv をパースするのはこの階層だけです。`app.py` は root `Click.Group` と `_LAZY_SUBCOMMANDS` レジストリを保持します。各エントリは **絶対モジュールパス**（`pdb2reaction.workflows.all`, `pdb2reaction.io.trj2fig`, …）を使うため、リゾルバは `default_group.py` 自体の所在に依存しません。`common_options.py` はサブコマンド間で共有される option-decorator ファクトリ（`@add_print_every_option`, `@add_irc_pos_def_option`, `@add_precision_option`, `@add_coord_type_option`, `@add_ml_charge_spin_options`）を集約します。サブコマンド本体はこれらのデコレータを `@click.pass_context` の上に積み重ね、`--help` テキストを揃えて保ちます。

**L2 `workflows/`**（17 ファイル）。サブコマンド 1 つにつき 1 ファイル。各ファイルは `cli` という名前の単一の `@click.command()` とそのプライベートヘルパーを所有します。大きなステージランナー（`all.py` = 4,979 LOC, `path_search.py` = 2,790 LOC, `tsopt.py` = 2,063 LOC, `extract.py` = 2,159 LOC）は、現在のレイアウトでは単一ファイルのまま残されています。

**L3 `domain/`**。`torch` / `numpy` / `pysisyphus.constants`（数値バックエンド）を import してよい化学的に意味を持つヘルパーロジックですが、MLIP ランタイム（`fairchem`, `orb_models`, `mace`, `aimnet`）を import しては **いけません**。`# DOMAIN_PURE` モジュール docstring マーカーと `.github/scripts/check_engineering_markers.py` がこの deny list を強制します。domain ヘルパーはどの L2 ステージランナーからでも再利用可能です。

**L4a `backends/`**（約 8 ファイル）。MLIP バックエンドディスパッチャ（`__init__.py` + `base.py`）と、サポートする各 MLIP につき 1 つのアダプタ（`uma.py`, `orb.py`, `mace.py`, `aimnet2.py`）。`solvent.py` と `xtb_alpb_correction.py` は xTB ALPB 暗黙溶媒デルタ補正（オプトインの MLIP ラッパー）を担います。`pdb2reaction` は純 MLIP のクラスターモデルパッケージであり、こちら側には ML/MM ONIOM 計算機コアも MM 力場カップリングも存在しません。

**L4b `io/`**。出力側の I/O に関する事項: ステージごとのサマリライタ、エネルギーダイアグラム、軌跡レンダリング、PDB altloc 修正、インメモリ Hessian キャッシュ。`io/` は `workflows/` に依存しません。出力フォーマットはここが所有し、ステージランナーが消費します。

**L5 `core/`**。最下層。`defaults.py` はすべての CLI デフォルトの **単一真実源** です。どこか他の場所に数値を追加する前に、まずここを grep してください。`utils.py` は PDB / XYZ / プロットヘルパーの約 2,560-LOC の寄せ集めです。

### 2.4 遅延 import の仕組み（概念図）

```text
External consumer                        Package root             Layer dir
---------------------------------------  ----------------------   ---------

from pdb2reaction.core.utils import x ──► (direct dotted import) ──► pdb2reaction/core/utils.py
import pdb2reaction.io.trj2fig        ──► (direct dotted import) ──► pdb2reaction/io/trj2fig.py

from pdb2reaction import <Symbol>     ──► pdb2reaction/__init__.py
                                          __getattr__("<Symbol>")
                                          └─► _LAZY_SYMBOLS["<Symbol>"]
                                              = "pdb2reaction.<layer>.<module>"
                                              └─► importlib.import_module(...)

from pdb2reaction import <module>     ──► pdb2reaction/__init__.py
        (= module attr)                   __getattr__("<module>")
                                          └─► _LAZY_MODULES["<module>"]
                                              = "pdb2reaction.<layer>.<module>"
                                              └─► importlib.import_module(...) returns module

pdb2reaction myaction                 ──► pdb2reaction/cli/app.py
                                          _LAZY_SUBCOMMANDS["myaction"]
                                          = ("pdb2reaction.workflows.myaction", "cli", "...")
                                          └─► importlib.import_module(absolute path)
                                              └─► getattr(module, "cli") → Click command
```

遅延 import 互換性の 2 階層と CLI ディスパッチ:

1. **Root シンボル属性**（`from pdb2reaction import <Symbol>`）— `pdb2reaction/__init__.py:_LAZY_SYMBOLS` + PEP 562 `__getattr__` が処理します。シンボルは初回アクセス時に layer-dir パスから読み込まれ、`pdb2reaction` import 時の import コストはゼロのまま保たれます。
2. **Root モジュール属性**（`from pdb2reaction import <module>`）— `_LAZY_MODULES` が処理します。`__getattr__` は `importlib.import_module` を介してモジュールオブジェクト自体を返します。`pdb2reaction` は現在、消費されるモジュール属性パスを 0 件しか持ちません（レジストリは空です。root 属性アクセスは将来の拡張のために予約されています）。

CLI サブコマンドリゾルバ（`cli/app.py:_LAZY_SUBCOMMANDS`）は **絶対** モジュールパス（例: `"pdb2reaction.workflows.all"`）を使うため、`default_group.py` を `cli/` に移動してもサブコマンド探索が静かに壊れることはありません（レジストリはもはや `__package__` に依存しません）。

---

## 3. Fresh-eyes 5 ステップナビゲーション（合計 ≈ 40 分）

リポジトリを初めて開くコントリビュータは、この道筋を上から下へ辿ってください。各ステップで 1 つの関心事を片付けます。

| step | minutes | open | what you learn |
|------|---------|------|-----------------|
| 1 | 3 | [`README.md`](https://github.com/t-0hmura/pdb2reaction/blob/main/README.md) | 1 段落のエレベーターピッチ + 単一コマンド使用法 |
| 2 | 5 | このファイル（`docs/architecture.md`）§2 + §4 | 6 階層のディレクトリツリー、依存方向、各関心事の所在 |
| 3 | 5 | [`pdb2reaction/cli/app.py`](../../pdb2reaction/cli/app.py) | Click root group、`_LAZY_SUBCOMMANDS` レジストリ（≈ 18 エントリ）、絶対パス解決 |
| 4 | 20 | [`pdb2reaction/workflows/all.py`](../../pdb2reaction/workflows/all.py)（4,147 LOC, 流し読み） | 1 つの完全なサブコマンドを上から下まで; `extract → MEP → tsopt → IRC → freq → dft` を追う |
| 5 | 7 | [`CONTRIBUTING.md`](https://github.com/t-0hmura/pdb2reaction/blob/main/CONTRIBUTING.md) §3 + §4 | 5 つの add-a-X レシピ + 「触るな」という隠れた制約 |

ステップ 5 の後は、§4 のファイル索引を辿ることで他のどのファイルでも読めます。本パッケージは意図的に **各階層内でフラット** です。`pdb2reaction/<layer>/` の下にネストしたパッケージは存在しないため、2 ディレクトリより深く辿る必要は決してありません。

---

## 4. ファイル索引 — 「この関心事はどこにあるか?」

### 4.1 CLI / エントリ（L1 `cli/`）

| concern | file |
|---|---|
| Click root group + サブコマンドディスパッチ | `pdb2reaction/cli/app.py` |
| サブコマンドリゾルバ（遅延 import） | `pdb2reaction/cli/default_group.py` |
| `python -m pdb2reaction` エントリ | `pdb2reaction/__main__.py` |
| 共有 option decorator ファクトリ | `pdb2reaction/cli/decorators.py` |
| `--help-advanced` ページャ | `pdb2reaction/cli/help_pages.py` |
| Bool flag 互換（`--flag` / `--no-flag` + value style） | `pdb2reaction/cli/bool_compat.py` |
| 共有 option-decorator ファクトリ（`--print-every`, `--irc-pos-def`, `--precision`, `--coord-type`, `--charge / --ligand-charge / --multiplicity`） | `pdb2reaction/cli/common_options.py` |

### 4.2 ワークフローステージランナー（L2 `workflows/`）

| concern | file |
|---|---|
| 全パイプラインオーケストレータ | `pdb2reaction/workflows/all.py` |
| 構造最適化（LBFGS / RFO） | `pdb2reaction/workflows/opt.py` |
| 1D / 2D / 3D スキャン + 共有 | `pdb2reaction/workflows/scan{,2d,3d,_common}.py` |
| MEP 探索（GSM） | `pdb2reaction/workflows/path_search.py` |
| MEP optimiser コア（pysisyphus COS） | `pdb2reaction/workflows/path_opt.py` |
| TS 最適化（RSIRFO + Bofill + macro/micro） | `pdb2reaction/workflows/tsopt.py` |
| 振動解析（PHVA + UMA active block） | `pdb2reaction/workflows/freq.py` |
| IRC 積分（macro / micro） | `pdb2reaction/workflows/irc.py` |
| 一点 DFT（gpu4pyscf サブプロセス） | `pdb2reaction/workflows/dft.py` |
| 活性部位抽出（クラスターキャップ） | `pdb2reaction/workflows/extract.py` |
| 拘束ヘルパー | `pdb2reaction/workflows/restraints.py` |
| Kabsch / frozen-subset アラインメント | `pdb2reaction/workflows/align_freeze.py` |

### 4.3 化学ヘルパー（L3 `domain/`）

| concern | file |
|---|---|
| R↔P 結合変化検出 | `pdb2reaction/domain/bond_changes.py` |
| IRC 後の結合サマリ | `pdb2reaction/domain/bond_summary.py` |
| PDB 元素列正規化 | `pdb2reaction/domain/add_elem_info.py` |

### 4.4 MLIP バックエンド（L4a `backends/`）

| concern | file |
|---|---|
| バックエンドディスパッチ + レジストリ | `pdb2reaction/backends/__init__.py` |
| `MLIPCalculator` プロトコル + base | `pdb2reaction/backends/base.py` |
| バックエンドごとのアダプタ | `pdb2reaction/backends/{uma, orb, mace, aimnet2}.py` |
| xTB ALPB 暗黙溶媒ヘルパー | `pdb2reaction/backends/solvent.py` |
| xTB ALPB デルタ補正 | `pdb2reaction/backends/xtb_alpb_correction.py` |

add-a-backend レシピは [Backends](backends.md) を参照してください。

### 4.5 I/O（L4b `io/`）

| concern | file |
|---|---|
| `summary.json` / `summary.md` ライタ | `pdb2reaction/io/summary.py` |
| Plotly エネルギーダイアグラム | `pdb2reaction/io/energy_diagram.py` |
| 軌跡 → PNG / SVG / PDF / HTML / CSV | `pdb2reaction/io/trj2fig.py` |
| PDB altloc 解決 | `pdb2reaction/io/pdb_fix.py` |
| インメモリ Hessian キャッシュ（実行ごとの TTL） | `pdb2reaction/io/hessian_cache.py` |
| 調和拘束のセットアップ | `pdb2reaction/workflows/restraints.py`（L2 ステージヘルパー） |

### 4.6 Foundation（L5 `core/`）

| concern | file |
|---|---|
| **すべての CLI デフォルト（単一真実源）** | `pdb2reaction/core/defaults.py` |
| PDB / XYZ / plot ヘルパー | `pdb2reaction/core/utils.py` |
| `-v` / `-vv` logging 配線 | `pdb2reaction/core/logging.py` |

### 4.7 repo 内部の同梱 fork

| dir | role | divergent files (do NOT replace with upstream) |
|---|---|---|
| `pysisyphus/` | optimiser / TS / IRC エンジン | `irc/IRC.py`（オプトインの `require_pos_def_hessian` PSD 収束ガード）、`optimizers/hessian_updates.py`（advanced index 上での Bofill scatter、GPU OOM 回避のための CPU 専用 `bofill_update` パス）、`tsoptimizers/{RSIRFOptimizer,RSPRFOptimizer,TRIM,TSHessianOptimizer}.py`、`calculators/{Calculator,Dimer}.py`、`_array.py`（torch/numpy バックエンド shim） |
| `thermoanalysis/` | thermochemistry（ΔG, ZPE, 分配関数） | `QCData.py`（上流との branding 差分） |

touch 制限の境界については各ディレクトリの `README.md` を参照してください。

---

## 5. 隠れた制約（パッチ前に必読）

### 5.1 化学ルール（grep レシピ）

正確性に直結する 3 つのルールが `backends/`、`workflows/`、`core/defaults.py` に散在しています。これらは smoke テストでは検出 **されません**。ここでの静かな drift は反応経路の精度を壊します。インラインの `# CHEMISTRY-RULE:N` マーカーと `# DOMAIN_PURE` モジュール docstring マーカーがルールを識別し、`.github/scripts/check_engineering_markers.py` が CI でマーカーの完全性を強制します。

編集前にすべての化学ルールを見つけるには:

```bash
# List all rule sites in the repo (host file + line)
grep -rnE '# CHEMISTRY-RULE:[0-9]+' pdb2reaction/

# List every # DOMAIN_PURE marker (= chemistry-rule host modules)
grep -rn '# DOMAIN_PURE' pdb2reaction/
```

3 つのルール（マーカー ID は連続していません）は次の通りです:

| marker | rule | host file |
|---|---|---|
| 4 | gpu4pyscf `rks_lowmem` triple-guard | `pdb2reaction/workflows/dft.py` |
| 5 | def2 family auto-ECP injection | `pdb2reaction/workflows/dft.py` |
| 7 | `bofill_update` advanced-indexing scatter | `pdb2reaction/workflows/tsopt.py` |

これらのいずれかを編集するには、`[CHEMISTRY-RULE:N]` コミットプレフィックスと HEAVY 層の numerical-golden ゲート通過が必要です（`CONTRIBUTING.md` §1.1 を参照）。先に DFT ペア（#4, #5）を学び、その後で TS scatter ルール（#7）を学んでください。

### 5.2 VRAM 管理の不変条件（`del` チェーンをリファクタしない）

IRC / TSopt / Freq ステージは、CUDA メモリを解放するためにステージ間で GPU 常駐オブジェクト（`calc`, `geom`, `hess`）を明示的に `del` します。`all` ワークフローはさらにステージ境界で `gc.collect()` を実行します。**これらの `del` / `gc.collect()` 文をリファクタで取り除かないでください**。大きな活性部位モデルを伴う長時間の `all` ジョブは、これらがないと OOM します。

### 5.3 同梱 fork: 上流を並べてインストールしない

同梱された `pysisyphus/` と `thermoanalysis/` パッケージは **fork** です。`pip install pysisyphus` または `pip install thermoanalysis` を本パッケージの隣に再インストールすると、次が静かに壊れます:

- `pysisyphus/irc/IRC.py` — 初期変位のメモリ衛生 + オプトインの `require_pos_def_hessian` kwarg
- `pysisyphus/optimizers/hessian_updates.py` — advanced index 上での Bofill scatter、GPU OOM 回避のための CPU 専用 `bofill_update` パス
- `pysisyphus/tsoptimizers/TSHessianOptimizer.py` — RSIRFO kwargs（fork 間でホストパッケージの import パスが分岐）
- `pysisyphus/calculators/{Calculator,Dimer}.py` — GPU 対応バックエンドフック（30 以上の QM バックエンドは削除済み。抽象 base + Dimer TS calculator のみ残存）
- `pysisyphus/_array.py` — `optimizers/hessian_updates.py`（および徐々に他のホットパスファイル）で使われる `get_xp` / `_outer` / `_dot` / `_eigh` shim
- `thermoanalysis/QCData.py` — 上流との branding / I/O 差分

### 5.4 `pyproject.toml` の配列は 0-diff

`[tool.setuptools.packages.find].include` と `dependencies` は、本リリース中は **0-diff 配列** として扱われます。`include` glob（`pdb2reaction*`）は新しい層サブパッケージをすでに自動探索します。`vendor/` や `internal/` のコンテナディレクトリを追加したり、新しいランタイム依存をピン留めしたりするとインストール契約が壊れ、リリーススコープで禁止されています。Reflow / コメント編集は問題ありません。**配列の内容** は凍結されています。

### 5.5 `_LAZY_SUBCOMMANDS` レジストリは絶対パスを使う必要がある

`pdb2reaction/cli/app.py:_LAZY_SUBCOMMANDS` はすべてのサブコマンドを **絶対** モジュールパスで解決します。いずれかのエントリを相対 dotted import（`".all"` など）に戻すと、`default_group.py` が移動した際にサブコマンド探索が静かに壊れます。リゾルバの `__package__` がパッケージルートから drift するためです。内部設計ノートを参照してください。

---

## 6. 同梱 fork（repo 内部）

`pdb2reaction` はリポジトリ最上位に **2 つ** の repo 内部モジュールを同梱しています:

| dir | upstream PyPI? | purpose | scope of edits allowed |
|---|---|---|---|
| `pysisyphus/` | NO — fork、`pip install pysisyphus` を並べないこと | optimiser, TS, IRC, COS, calculators | 本リリースラインでは annotation のみ（docstring + 型ヒント）; ロジック編集は禁止 |
| `thermoanalysis/` | NO — fork（branding 差分） | ΔG, ZPE, 分配関数, `QCData` | `pysisyphus/` と同じ |

各ディレクトリは分岐ファイルと touch 制限の境界を列挙した独自の `README.md` を持ちます。階層モデルから見ると、これらの fork は L1..L5 グラフの **外側** に位置します。どの階層も絶対パッケージパス（`from pysisyphus.X import Y`）でこれらを import でき、`L1 → L2 → {L3, L4} → L5` の方向を壊しません。

---

## 7. 推奨される深掘りの読み順

Fresh-eyes ツアー（§3）の後は、この深さ優先の読み順に従ってください:

1. `pdb2reaction/core/defaults.py` — デフォルト値テーブルを頭に入れる; 下流のすべてはここから読みます。
2. `pdb2reaction/cli/app.py` — Click root + `_LAZY_SUBCOMMANDS` レジストリ。
3. `pdb2reaction/workflows/all.py` — 1 つの完全なパイプラインを上から下まで。
4. `pdb2reaction/workflows/extract.py` — 活性部位クラスターキャップ。
5. `pdb2reaction/backends/__init__.py` + `base.py` — MLIP ディスパッチャとバックエンドごとのアダプタ契約。
6. `pdb2reaction/workflows/tsopt.py` — RS-I-RFO + Bofill scatter（CHEMISTRY-RULE:7）。
7. `pdb2reaction/workflows/freq.py` — クラスターモデル上での振動解析。
8. `pdb2reaction/workflows/irc.py` — VRAM 衛生 + IRC 積分。
9. `pdb2reaction/workflows/dft.py` — gpu4pyscf による一点 DFT（CHEMISTRY-RULE:4 + :5）。
10. `pdb2reaction/core/utils.py` — 共有 PDB / XYZ / plot ヘルパー。
