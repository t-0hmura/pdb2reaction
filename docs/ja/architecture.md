# アーキテクチャ: pdb2reaction

## 1. 概要

pdb2reaction のコードを直す人向けに、パッケージの層、ファイルの置き場、直す前に守る制約をまとめたページです。直したあとは、CONTRIBUTING の [gate cycle](https://github.com/t-0hmura/pdb2reaction/blob/main/CONTRIBUTING.md#11-gate-cycle) の検査を走らせてください。計算を実行するだけなら、[はじめに](getting-started.md)から読んでください。

`pdb2reaction` は、活性部位のクラスターモデルで酵素反応の経路を解析する Python の CLI です。構造と経路の段には組み込みの MLIP か自作の ASE calculator を使い、PySCF/GPU4PySCF による DFT の一点計算も追加できます。

---

## 2. 階層構造（6 つの計算層）

任意の MCP サーバー `pdb2reaction/mcp/` は 6 層の外にあります。中身は `server.py`・`_tools.py`・`_runner.py` です。[構成と使い方](mcp_server.md)を参照してください。

### 2.1 階層テーブル

6 つの層は L1〜L5 で、L4 は L4a と L4b に分かれます。

| 層 | ディレクトリ | 責務 | 依存してよい先 |
|---|---|---|---|
| **L1 Interface** | `pdb2reaction/cli/` | Click root group、共有 option-decorator ファクトリ（`common_options.py`）、`--help-advanced`、bool flag 正規化、サブコマンドリゾルバ | `workflows/`, `core/` |
| **L2 Application** | `pdb2reaction/workflows/` | サブコマンドごとのオーケストレーションと、段で共有する helper（`all.py`, `path_search.py`, `tsopt.py`, `extract.py`, `irc.py`, `freq.py`, `dft.py`, …） | `domain/`, `backends/`, `io/`, `core/` |
| **L3 Domain** | `pdb2reaction/domain/` | 化学的に意味を持つヘルパーロジック（結合変化検出、結合サマリ、元素情報伝播） | `core/` |
| **L4a Infra (MLIP)** | `pdb2reaction/backends/` | MLIP バックエンドディスパッチャ + バックエンドごとのアダプタ（UMA / Orb / MACE / AIMNet2） | `core/` |
| **L4b Infra (I/O)** | `pdb2reaction/io/` | 出力レイアウト、サマリ、軌跡、PDB 修正、エネルギーダイアグラム、Hessian キャッシュ | `core/` |
| **L5 Foundation** | `pdb2reaction/core/` | defaults（共有の既定値の主な出典）、utils（structure / coordinate / plot helper）、logging、output、result の公開 | (none、設計意図) |
| (bundle, not a layer) | `<repo>/pysisyphus/`, `<repo>/thermoanalysis/` | repo 内部 fork（optimizer / thermochemistry） | (sibling, layer-external) |

**依存の向き（設計目標）**: `L1 → L2 → {L3, L4} → L5`。表の最後の列がこの決まりです。今はこれを破る import があります: `workflows/* → cli`、`cli/common_options.py → backends`、`core/utils.py → domain`・`io`・`backends`・`cli` と `core/defaults.py → backends`、`domain/add_elem_info.py`・`domain/bond_summary.py` といくつかの `io/` のモジュール → `cli`、`io/charge.py`・`io/structure_formats.py → domain`、`io/trj2fig.py → backends`。同梱の fork は層の外にあり、どの層からも `from pysisyphus.X import Y` の形で import できます。

CI が検査するのはこの向きの一部だけです。どちらのスクリプトもリポジトリの root で `python` で実行し、終了コードが 0 なら合格です。

- `.github/scripts/check_import_graph.py` は、`pdb2reaction/*` のモジュール間の import の循環、`core/` と `domain/` から `workflows/*` への import、`pysisyphus/` から `pdb2reaction` への import を禁止します。
- `.github/scripts/check_engineering_markers.py` は、`# CHEMISTRY-RULE` と `# DOMAIN_PURE` のマーカー（§5.1）と、MLIP の SDK の `fairchem`、`orb_models`、`mace`、`aimnet` を `backends/` の下でだけ import していることを検査します。

### 2.2 パッケージツリーの ASCII マップ

```
pdb2reaction/ [GH: t-0hmura/pdb2reaction]
├── pyproject.toml packages.find = ["pdb2reaction*", ...] (glob, frozen)
├── README.md / CONTRIBUTING.md / CHANGELOG.md
├── docs/
│ ├── ja/architecture.md ← このファイル
│ └──... (Sphinx ドキュメントサイト)
├── pdb2reaction/ ← package body, 6-layer physical dir
│ ├── __init__.py PEP 562 lazy: _LAZY_SYMBOLS / _LAZY_MODULES + __getattr__
│ ├── __main__.py `from .cli import cli`
│ ├── _version.py / py.typed
│ │
│ ├── cli/ # === L1 Interface ===
│ │ ├── app.py Click group + _LAZY_SUBCOMMANDS registry (absolute paths)
│ │ ├── common_options.py @add_print_every_option / @add_irc_pos_def_option / @add_precision_option / @add_coord_type_option / @add_ml_charge_spin_options
│ │ ├── decorators.py resolve_yaml_sources / load_merged_yaml_cfg / _write_error_json
│ │ ├── help_pages.py --help-advanced pager
│ │ ├── bool_compat.py --flag / --no-flag normalization
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
│ │ ├── add_elem_info.py PDB element column normalizer
│ │ └── residue_data.py 残基・イオン電荷table
│ │
│ ├── backends/ # === L4a Infra (MLIP) ===
│ │ ├── __init__.py backend dispatch + registry
│ │ ├── base.py MLIPCalculator protocol
│ │ ├── custom.py custom ASE-calculator adapter
│ │ ├── _determinism.py deterministic reduction shim
│ │ ├── pyscf_dft.py 任意の PySCF/GPU4PySCF adapter
│ │ └── uma.py / orb.py / mace.py / aimnet2.py MLIP adapters
│ │
│ ├── io/ # === L4b Infra (I/O) ===
│ │ ├── summary.py summary.json / summary.log writer
│ │ ├── energy_diagram.py Plotly diagram
│ │ ├── trj2fig.py trajectory → PNG / JPEG / SVG / PDF / HTML / CSV
│ │ ├── pdb_fix.py altloc resolution
│ │ ├── altloc.py altloc identity／選択helper
│ │ ├── charge.py 残基を考慮した電荷解決
│ │ ├── structure_formats.py PDB/mmCIF 変換 + identifier 復元
│ │ └── hessian_cache.py in-memory Hessian cache
│ │
│ └── core/ # === L5 Foundation ===
│   ├── defaults.py 共有の既定値の主な出典
│   ├── dft_settings.py DFT の設定の解決（CHEMISTRY-RULE:4）
│   ├── logging.py -v / --verbose の配線
│   ├── utils.py PDB / XYZ / plot helpers
│   ├── output.py / result_commit.py output／result ownership
│   └── pes_composition.py energy成分合成
│
├── tests/ smoke / unit
├── .github/ workflows/ + scripts/ (CI、release、engineering、documentation checks)
└── (repo-top sibling, layer-external bundled forks)
 pysisyphus/ 98 file、repo 内部 fork（ここで使う数値計算の機能だけ）
 thermoanalysis/ 7 file、repo 内部 fork
```

### 2.3 階層ごとの責務詳細

**L1 `cli/`** は、root group、遅延 import の registry、argv の正規化、段階的な help、共有の option-decorator を受け持ちます。実際の `@click.command()` とコマンド固有のオプションは workflow や utility のモジュールの側にあり、root はそれらを見つけて呼び出します。`_LAZY_SUBCOMMANDS` の各 entry は **絶対モジュールパス**を使います。

**L2 `workflows/`** には、計算コマンドのモジュールと、`scan_common.py`、`restraints.py`、`align_freeze.py` などの段で共有する helper があります。コマンドのモジュールはふつう `cli` という名前の `@click.command()` を 1 つ公開しますが、utility のコマンドは `io/` と `domain/` にもあるので、「サブコマンドごとに必ず 1 ファイル」「workflow のファイルはすべてコマンド」という決まりはありません。

**L3 `domain/`**。`torch` / `numpy` / `pysisyphus.constants`（数値バックエンド）を import してよい化学的に意味を持つヘルパーロジックですが、MLIP の SDK を import しては **いけません**。domain ヘルパーはどの L2 ステージランナーからでも再利用可能です。

**L4a `backends/`**。MLIP バックエンドディスパッチャと、サポートする各 MLIP につき 1 つのアダプタがあります。ディスパッチャは自作の ASE calculator の受け口も持ちます。

**L4b `io/`**。出力側の I/O を受け持ちます。段ごとの summary の書き出し、エネルギーダイアグラム、軌跡の描画、PDB の altloc の修正、PDB/mmCIF の変換とテンプレートの identifier の復元、メモリ上の Hessian キャッシュです。出力の形式はここが持ち、段のランナーがそれを使います。foundation 以外への import は §2.1 に挙げています。

**L5 `core/`** は最下層です。`defaults.py` は、共有される数値と CLI の既定値の **主な出典** です。まずここを grep し、そのあと path engine の選択のような、理由のあるコマンド固有の既定値も確かめます。`utils.py` には、設定・構造・座標・plot の共有の helper があります。

### 2.4 遅延 import の仕組み（概念図）

```text
External consumer                        Package root             Layer dir
---------------------------------------  ----------------------   ---------

from pdb2reaction.core.utils import x ──► (direct dotted import) ──► pdb2reaction/core/utils.py
import pdb2reaction.io.trj2fig        ──► (direct dotted import) ──► pdb2reaction/io/trj2fig.py

pdb2reaction myaction                 ──► pdb2reaction/cli/app.py
                                          _LAZY_SUBCOMMANDS["myaction"]
                                          = ("pdb2reaction.workflows.myaction", "cli", "...")
                                          └─► importlib.import_module(absolute path)
                                              └─► getattr(module, "cli") → Click command
```

`pdb2reaction/__init__.py` には、root から再 export するための PEP 562 の registry `_LAZY_SYMBOLS` と `_LAZY_MODULES` もあります。どちらも空なので、シンボルは上の 2 行のように層のモジュールから import します。

---

## 3. 初めて読む人のための 5 ステップナビゲーション（合計 ≈ 40 分）

リポジトリを初めて開くコントリビュータは、この道筋を上から下へ辿ってください。

| ステップ | 分 | 開くもの | 分かること |
|------|---------|------|-----------------|
| 1 | 3 | [`README.md`](https://github.com/t-0hmura/pdb2reaction/blob/main/README.md) | 1 段落のエレベーターピッチ + 単一コマンド使用法 |
| 2 | 5 | このファイル（`docs/ja/architecture.md`）§2 + §4 | 6 階層のディレクトリツリー、依存方向、各関心事の所在 |
| 3 | 5 | [`pdb2reaction/cli/app.py`](https://github.com/t-0hmura/pdb2reaction/blob/main/pdb2reaction/cli/app.py) | Click root group、`_LAZY_SUBCOMMANDS` registry、絶対path解決 |
| 4 | 20 | [`pdb2reaction/workflows/all.py`](https://github.com/t-0hmura/pdb2reaction/blob/main/pdb2reaction/workflows/all.py)（流し読み） | 1 つのサブコマンドを上から下まで追う。`extract → MEP → tsopt → IRC → freq → dft` をたどる |
| 5 | 7 | [`CONTRIBUTING.md`](https://github.com/t-0hmura/pdb2reaction/blob/main/CONTRIBUTING.md) §3 + §4 | 機能を足す 5 つの手順（サブコマンド、MLIP バックエンド、出力形式、workflow の段、テスト）+ 触ってはいけない箇所の一覧 |

ステップ 5 の後は、§4 のファイル索引を辿ることで他のどのファイルでも読めます。本パッケージは **各階層内でフラット** です。`pdb2reaction/<layer>/` の下にネストしたパッケージは存在しないため、2 ディレクトリより深く辿る必要は決してありません。

---

## 4. ファイル索引 — 「この関心事はどこにあるか?」

### 4.1 CLI / エントリ（L1 `cli/`）

| 関心事 | ファイル |
|---|---|
| Click root group + サブコマンドディスパッチ | `pdb2reaction/cli/app.py` |
| サブコマンドリゾルバ（遅延 import） | `pdb2reaction/cli/default_group.py` |
| `python -m pdb2reaction` エントリ | `pdb2reaction/__main__.py` |
| YAML ソース解決 + 標準化された例外処理 | `pdb2reaction/cli/decorators.py` |
| `--help-advanced` ページャ | `pdb2reaction/cli/help_pages.py` |
| Bool flag 互換（`--flag` / `--no-flag` + value style） | `pdb2reaction/cli/bool_compat.py` |
| 共有 option-decorator ファクトリ（`--print-every`, `--irc-pos-def`, `--precision`, `--coord-type`, `--charge / --ligand-charge / --multiplicity`） | `pdb2reaction/cli/common_options.py` |

### 4.2 ワークフローステージランナー（L2 `workflows/`）

| 関心事 | ファイル |
|---|---|
| 全パイプラインオーケストレータ | `pdb2reaction/workflows/all.py` |
| 構造最適化（L-BFGS / RFO） | `pdb2reaction/workflows/opt.py` |
| Scanと2D/3D energy-landscape grid + 共有 | `pdb2reaction/workflows/scan{,2d,3d,_common}.py` |
| MEP 探索（GSM / DMF） | `pdb2reaction/workflows/path_search.py` |
| MEP optimizer コア（pysisyphus COS） | `pdb2reaction/workflows/path_opt.py` |
| TS 最適化（RS-P-RFO / RS-I-RFO / TRIM / Dimer + Bofill） | `pdb2reaction/workflows/tsopt.py` |
| 振動解析（backend-agnostic PHVA + active block） | `pdb2reaction/workflows/freq.py` |
| IRC 積分 | `pdb2reaction/workflows/irc.py` |
| 一点 DFT（PySCF / GPU4PySCF、同じプロセス内） | `pdb2reaction/workflows/dft.py` |
| 活性部位抽出（クラスターキャップ） | `pdb2reaction/workflows/extract.py` |
| 拘束ヘルパー | `pdb2reaction/workflows/restraints.py` |
| Kabsch / frozen-subset アラインメント | `pdb2reaction/workflows/align_freeze.py` |

### 4.3 化学ヘルパー（L3 `domain/`）

| 関心事 | ファイル |
|---|---|
| R↔P 結合変化検出 | `pdb2reaction/domain/bond_changes.py` |
| IRC 後の結合サマリ | `pdb2reaction/domain/bond_summary.py` |
| PDB 元素列正規化 | `pdb2reaction/domain/add_elem_info.py` |
| 残基・イオン電荷テーブル | `pdb2reaction/domain/residue_data.py` |

### 4.4 MLIP バックエンド（L4a `backends/`）

| 関心事 | ファイル |
|---|---|
| バックエンドディスパッチ + レジストリ | `pdb2reaction/backends/__init__.py` |
| `MLIPCalculator` プロトコル + base | `pdb2reaction/backends/base.py` |
| custom ASE calculator アダプタ | `pdb2reaction/backends/custom.py` |
| 決定論的 reduction shim | `pdb2reaction/backends/_determinism.py` |
| 任意の DFT calculator adapter | `pdb2reaction/backends/pyscf_dft.py` |
| バックエンドごとのアダプタ | `pdb2reaction/backends/{uma, orb, mace, aimnet2}.py` |

バックエンドを追加するときは、`CONTRIBUTING.md` の recipe 3.2「Add an MLIP backend」に従ってください。

### 4.5 I/O（L4b `io/`）

| 関心事 | ファイル |
|---|---|
| `summary.json` / `summary.log` ライタ | `pdb2reaction/io/summary.py` |
| Plotly エネルギーダイアグラム | `pdb2reaction/io/energy_diagram.py` |
| 軌跡 → PNG / JPEG / SVG / PDF / HTML / CSV | `pdb2reaction/io/trj2fig.py` |
| PDB altloc 解決 | `pdb2reaction/io/pdb_fix.py` |
| altloc identity と一貫した選択 | `pdb2reaction/io/altloc.py` |
| 残基を考慮した電荷解決 | `pdb2reaction/io/charge.py` |
| PDB/mmCIF 変換 + identifier 復元 | `pdb2reaction/io/structure_formats.py` |
| インメモリ Hessian キャッシュ（実行ごとの TTL） | `pdb2reaction/io/hessian_cache.py` |

### 4.6 Foundation（L5 `core/`）

| 関心事 | ファイル |
|---|---|
| **共有の数値の既定値（主な出典。コマンド固有の例外も確かめる）** | `pdb2reaction/core/defaults.py` |
| PDB / XYZ / plot ヘルパー | `pdb2reaction/core/utils.py` |
| `-v` / `--verbose LEVEL` logging 配線 | `pdb2reaction/core/logging.py` |
| output／result ownership helper | `pdb2reaction/core/output.py`、`pdb2reaction/core/result_commit.py` |
| energy 成分合成 | `pdb2reaction/core/pes_composition.py` |

### 4.7 repo 内部の同梱 fork

| ディレクトリ | 役割 | 主な分岐ファイル（完全な一覧は各ディレクトリの README を参照） |
|---|---|---|
| `pysisyphus/` | optimizer / TS / IRC エンジン | `irc/IRC.py`（オプトインの `require_pos_def_hessian` PSD 収束ガード）、`optimizers/hessian_updates.py`（GPU 常駐の in-place rank-two Bofill 更新と、明示的な `PYSIS_BOFILL_CPU_OFFLOAD=1` フォールバック）、`tsoptimizers/{RSIRFOptimizer,RSPRFOptimizer,TRIM,TSHessianOptimizer}.py`、`calculators/{Calculator,Dimer}.py`、`_array.py`（torch/numpy バックエンド shim） |
| `thermoanalysis/` | thermochemistry（ΔG, ZPE, 分配関数） | `QCData.py`（上流との branding 差分） |

---

## 5. 隠れた制約（パッチ前に必読）

### 5.1 化学ルール（grep レシピ）

正確性に直結する下の 3 つのルールは、smoke テストでは検出 **されません**。ここでの静かな乖離は反応経路の精度を壊します。インラインの `# CHEMISTRY-RULE:N` マーカーがルールを識別し、`.github/scripts/check_engineering_markers.py` が CI でマーカーの完全性を強制します。

編集前にすべての化学ルールを見つけるには:

```bash
# List all rule sites in the repo (host file + line)
grep -rnE '# CHEMISTRY-RULE:[0-9]+' pdb2reaction/

# List every # DOMAIN_PURE marker (on workflows/dft.py, tsopt.py, sp.py: modules that must not depend on a specific MLIP backend)
grep -rn '# DOMAIN_PURE' pdb2reaction/
```

CI が求めるマーカーは次の 3 つです:

| マーカー | ルール | 実装のファイル |
|---|---|---|
| 4 | gpu4pyscf `rks_lowmem` のclosed-shell/GPU/lowmem guard | `pdb2reaction/core/dft_settings.py` |
| 5 | def2 family auto-ECP injection | `pdb2reaction/workflows/dft.py` |
| 7 | Dimer の flatten loop での active Hessian block の Bofill 更新（`_bofill_update_active`）。advanced indexing は in-place の `+=` ではなく代入で書く。optimizer 自体の Bofill 更新は `pysisyphus/optimizers/hessian_updates.py` にある | `pdb2reaction/workflows/tsopt.py` |

これらのいずれかを編集するには、コミットメッセージの先頭の `[CHEMISTRY-RULE:N]`、メンテナの承認、該当する回帰テスト、記録を残した定期の数値ベンチマークが必要です。`CONTRIBUTING.md` §4.1 と §5 を参照してください。

### 5.2 VRAM 管理の不変条件（`del` チェーンをリファクタしない）

IRC / TSopt / Freq の各ステージは、CUDA メモリを解放するために、それぞれの段の境界で自分の GPU 常駐の `calc`・`geom`・`hess` と一時オブジェクトを明示的に `del` します。3 つの workflow が同じ `del` の並びを共有しているわけではありません。ステージ境界では `gc.collect()` に加え、CUDA allocation がある場合は `torch.cuda.empty_cache()` も実行し、`all` は子ステージの間でガベージコレクションを行います。**これらの解放処理をリファクタで取り除かないでください**。大きな活性部位モデルを伴う長時間の `all` ジョブは、これらがないと OOM（メモリ不足）で止まります。

### 5.3 同梱 fork: 上流を並べてインストールしない

同梱された `pysisyphus/` と `thermoanalysis/` パッケージは **fork** です。本パッケージと並べて `pip install pysisyphus` や `pip install thermoanalysis` で PyPI 版を入れ直すと、次が気づかないうちに動かなくなります:

- `pysisyphus/irc/IRC.py` — 初期変位のメモリ管理 + オプトインの `require_pos_def_hessian` kwarg
- `pysisyphus/optimizers/hessian_updates.py` — GPU 常駐の in-place rank-two Bofill 更新、オプトインの `PYSIS_BOFILL_CPU_OFFLOAD=1` フォールバック
- `pysisyphus/tsoptimizers/TSHessianOptimizer.py` — RSIRFO kwargs
- `pysisyphus/calculators/{Calculator,Dimer}.py` — GPU 対応バックエンドフック。fork にあるのはこの 2 つだけで、QM の calculator はありません
- `pysisyphus/_array.py` — `optimizers/hessian_updates.py` やほかのいくつかのホットパスファイルで使われる `get_xp` / `_outer` / `_dot` / `_eigh` shim
- `thermoanalysis/QCData.py` — 上流との branding / I/O 差分

### 5.4 パッケージ変更には隔離インストール検証が必要

`include` glob の `pdb2reaction*` は新しい層サブパッケージを自動探索するため、通常の内部ファイル追加で package discovery を変える必要はありません。package discovery や runtime dependency を変更する場合は、sdist と wheel を作成し、wheel の内容を検査し、クリーン環境でインストールした後に CLI smoke test を実行してください。これはついでのリファクタではなく、配布物の中身を変える変更として扱ってください。

### 5.5 `_LAZY_SUBCOMMANDS` レジストリは絶対パスを使う必要がある

`pdb2reaction/cli/app.py:_LAZY_SUBCOMMANDS` はすべてのサブコマンドを **絶対** モジュールパスで解決します。いずれかのエントリを相対 dotted import（`".all"` など）に変えると、`default_group.py` が移動した際にサブコマンド探索が気づかないうちに動かなくなります。リゾルバの `__package__` がパッケージルートから乖離するためです。

---

## 6. 同梱 fork（repo 内部）

`pdb2reaction` はリポジトリ最上位に **2 つ** の repo 内部モジュールを同梱しています:

| ディレクトリ | 上流の PyPI 版か | 用途 | 許される編集の範囲 |
|---|---|---|---|
| `pysisyphus/` | NO — fork、`pip install pysisyphus` を並べないこと | optimizer, TS, IRC, COS, calculators | logic edit には再現された不具合または承認済み機能、focused regression test、該当する numerical/GPU benchmark が必要 |
| `thermoanalysis/` | NO — fork（branding/I/O 差分） | ΔG, ZPE, 分配関数, `QCData` | logic edit には I/O または数値上の必要性と thermochemistry golden test が必要 |

各ディレクトリは分岐ファイルと変更の方針を列挙した独自の `README.md` を持ちます。

---

## 7. 推奨される深掘りの読み順

初めて読む人のための道筋（§3）の後は、この深さ優先の読み順に従ってください:

1. `pdb2reaction/core/defaults.py` — 共有の既定値の主な出典（§2.3）。
2. `pdb2reaction/workflows/extract.py` — 活性部位クラスターキャップ。
3. `pdb2reaction/backends/__init__.py` + `base.py` — MLIP ディスパッチャとバックエンドごとのアダプタ契約。
4. `pdb2reaction/workflows/tsopt.py` — TS 最適化の driver と、Dimer の flatten loop の Bofill 更新（CHEMISTRY-RULE:7）。
5. `pdb2reaction/workflows/freq.py` — クラスターモデル上での振動解析。
6. `pdb2reaction/workflows/irc.py` — VRAM 管理 + IRC 積分。
7. `pdb2reaction/workflows/dft.py` — PySCF / GPU4PySCF による一点 DFT（CHEMISTRY-RULE:5）。
8. `pdb2reaction/core/utils.py` — 共有 PDB / XYZ / plot ヘルパー。
