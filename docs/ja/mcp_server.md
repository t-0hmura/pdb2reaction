# pdb2reaction MCP サーバー

AI エージェントから MCP（Model Context Protocol）で pdb2reaction の 18 個のツールを呼ぶための、インストール、ツールの一覧、クライアントの設定をまとめたページです。サーバー `pdb2reaction-mcp`（別名 `p2r-mcp`）は stdio 上の JSON-RPC でやりとりするので、[MCP](https://modelcontextprotocol.io/) に対応したどのクライアントからも使えます。Claude Desktop、Claude Code、Cursor、Codeium のほか、公式の Python や TypeScript の MCP SDK で作ったエージェントからも使えます。

## インストール

```bash
pip install "pdb2reaction[mcp]"
```

これにより `mcp[cli]` 依存関係が追加され、`pdb2reaction-mcp` / `p2r-mcp` のコマンドが登録されます。

## ツール

18 個のツールがあり、それぞれが CLI サブコマンドに 1 対 1 で対応します。各ツールは次のフィールドを持つ構造化された dict を返します。

- `schema_version`: 結果の形式の版
- `execution_status`: `completed` | `failed`
- `scientific_status`: `success` | `partial` | `failed`
- `summary_status`: 中身のある `summary` が返るのは `ok` のときだけです
  - `ok`: この呼び出しの `summary` を読めた
  - `not_required`: 要約を書かない構造 / I/O ヘルパー
  - `summary_missing`: `out_dir` に `summary.json` が無い
  - `summary_parse_error`: `summary.json` を JSON の object として読めない
  - `summary_run_mismatch`: ファイルが別の実行のもの
- `exit_code`: CLI のプロセスの終了コード
- `out_dir`: ステージランナーとスキャン / 経路 / パイプラインのツールの出力ディレクトリ。構造 / I/O ヘルパーでは null
- `summary`: 読み込んだ `summary.json`。構造 / I/O ヘルパーでは空の object
- `stderr_tail` / `stdout_tail`: プロセス出力の末尾約 60 行
- `hint`: CLI のエラーメッセージの末尾の `; recover: <hint>` にある対処のヒント（ある場合）
- `argv`: 実行したコマンドライン全体
- `run_id`: この呼び出しの UUID

各ツールは CLI をサブプロセスで実行し、作業ディレクトリは呼び出した側のままなので、入力の相対パスの意味は変わりません。各ツールの必須の引数は下の表にあります。すべての引数とその型は、クライアントがツールの一覧と一緒に受け取る入力スキーマにあります。

- `input_pdb`・`ts_pdb`・`reactant_pdb` などの入力のパスはそのコマンドの `-i` に渡るので、XYZ など `-i` が読む形式を渡せます。
- 省略できる引数はそのコマンドの CLI オプションを指定します。例えば `charge` は `-q`、`max_cycles` は `--max-cycles` です。

### 構造化されたエラーエンベロープ

ステージランナーとスキャン / 経路 / パイプラインのツールが失敗すると、返された `summary` に次のエラーのフィールドが入ります。エージェントはテキストをパースせずに、エラーのクラスで場合分けできます。構造 / I/O ヘルパーと、`summary.json` を書く前に止まった実行（`summary_missing`）にはエラーのフィールドが無いので、`stderr_tail` と `hint` を読んでください。

- `error`: エラーメッセージ
- `error_type`: 例外クラス名
- `error_class_chain`: そのクラスと親クラスの名前を、具体的なものから順に並べたもの
- `error_module`: 例外クラスが定義されているモジュール
- `error_label`: 上位レベルの CLI ステージラベル

### ステージランナー

| MCP ツール | 必須の引数 | CLI サブコマンド | 目的 |
|---|---|---|---|
| `optimize_geometry` | `input_pdb` | `pdb2reaction opt` | 単一の分子構造を最適化 |
| `find_transition_state` | `ts_pdb` | `pdb2reaction tsopt` | TS 探索（RS-P-RFO / Dimer / TRIM / RS-I-RFO） |
| `run_irc` | `ts_pdb` | `pdb2reaction irc` | TS 構造からの IRC 積分 |
| `compute_frequencies` | `input_pdb` | `pdb2reaction freq` | 振動解析 + 熱化学 |
| `run_single_point` | `input_pdb` | `pdb2reaction sp` | 選んだ `backend` での一点エネルギー + 原子間力（+ オプションで Hessian） |

### スキャン / 経路 / パイプライン

| MCP ツール | 必須の引数 | CLI サブコマンド | 目的 |
|---|---|---|---|
| `scan_1d` / `scan_2d` / `scan_3d` | `input_pdb`, `scan_lists` | `pdb2reaction scan` / `pdb2reaction scan2d` / `pdb2reaction scan3d` | 拘束をかけた距離・角度・二面角のスキャン |
| `optimize_path` | `reactant_pdb`, `product_pdb` | `pdb2reaction path-opt` | 2 端点間の MEP 最適化 |
| `search_paths` | `input_pdb`, `product_pdb` | `pdb2reaction path-search` | 再帰的な反応経路探索 |
| `run_full_pipeline` | `reactant_pdb` | `pdb2reaction all` | エンドツーエンド: extract → MEP → TS → IRC → freq → DFT |
| `run_single_point_dft` | `input_pdb` | `pdb2reaction dft` | 一点 DFT のエネルギーと原子電荷（GPU4PySCF または PySCF） |

### 構造 / I/O ヘルパー

| MCP ツール | 必須の引数 | CLI サブコマンド | 目的 |
|---|---|---|---|
| `extract_active_site` | `complex_pdb`, `ligand_id`, `radius_angstrom`, `output_pdb` | `pdb2reaction extract` | 活性部位モデル: リガンドの近くの残基を切り出し、キャップ水素を付ける |
| `add_element_info` | `input_pdb`, `output_pdb` | `pdb2reaction add-elem-info` | PDB の元素欄を修復 |
| `fix_altloc` | `input_pdb`, `output_pdb` | `pdb2reaction fix-altloc` | PDB の代替位置（altloc）を解決 |
| `plot_trajectory` | `input_trj_xyz`, `output_png` | `pdb2reaction trj2fig` | エネルギープロファイル図（デフォルトは PNG。JPEG/SVG/PDF/HTML/CSV も可） |
| `plot_energy_diagram` | `energies`, `output_png` | `pdb2reaction energy-diagram` | 与えた値からの状態エネルギー図 |
| `detect_bond_changes` | `reactant_pdb`, `product_pdb` | `pdb2reaction bond-summary` | 2 つの構造（XYZ / PDB / mmCIF / GJF）の間の結合変化 |

### 電荷と順序付き入力

`charge` と `ligand_charge` を持つツールで、PDB の残基名から全電荷を求めたいときは、`charge` を省いて残基名ごとの `ligand_charge` を渡します。PDB の残基情報が無い XYZ の入力には、全電荷の `charge` を明示してください。有効な GJF の電荷・多重度の行があれば、上書きしない限り両方の値をそこから取ります。`charge` と `ligand_charge` の両方を渡すと、明示した `charge` が優先されます。

`search_paths` では、反応物の `input_pdb` と `product_pdb` の間に入る中間体を、順番に並べて `intermediate_pdbs` に渡します。`scan_1d` と `run_full_pipeline` で段階的にスキャンするときは、最初の段を `scan_lists` に、後の段を `additional_scan_stages` に入れます。CLI には `--scan-lists` が 1 回だけ渡され、その後にすべての段の値が続きます。

## IRC と TS 最適化の設定

IRC と TS の引数は同じ名前の CLI オプションです。意味とデフォルトは各コマンドのページにあります。

- `run_irc` — `step_size`・`never_stop`・`irc_pos_def`: [`irc`](irc.md)。`--irc-pos-def` は [自動生成のオプションの一覧（英語のみ）](../reference/commands/irc.md) にあります
- `find_transition_state` — `opt_mode`: [`tsopt`](tsopt.md) の `--opt-mode`。デフォルトは `hess`（RS-P-RFO）です。{ref}`コマンドごとの --opt-mode <ja-opt-mode-semantics>` も参照してください
- `run_full_pipeline` — `irc_step_size`・`irc_never_stop`・`flatten`・`refine_path`: [`all`](all.md)。`--irc-step-size` と `--irc-never-stop` は [自動生成のオプションの一覧（英語のみ）](../reference/commands/all.md) にあります

## クライアント設定

クライアントごとに設定スキーマは異なります。次のスニペットはトップレベルの `mcpServers` オブジェクトを受け付けるクライアント用です。

- Claude Desktop — `~/Library/Application Support/Claude/claude_desktop_config.json`（macOS） / `%APPDATA%\Claude\claude_desktop_config.json`（Windows）
- Cursor — `~/.cursor/mcp.json`
- Claude Code — ファイルは編集せず、`claude mcp add pdb2reaction -- pdb2reaction-mcp` を実行します。登録できると、`claude mcp list` でこのサーバーに `✔ Connected` が出ます
- その他のクライアント — 各クライアント自身の MCP サーバードキュメントを参照

クライアントがサーバーを起動すると、ツールの一覧に 18 個のツールが出ます。

```json
{
  "mcpServers": {
    "pdb2reaction": {
      "command": "pdb2reaction-mcp",
      "args": []
    }
  }
}
```

環境変数 PATH と CUDA_VISIBLE_DEVICES を設定する完全な例は [`examples/mcp_client_config.json`](https://github.com/t-0hmura/pdb2reaction/blob/main/examples/mcp_client_config.json) を参照してください。

VS Code は `.vscode/mcp.json` で[トップレベルの `servers` オブジェクト](https://code.visualstudio.com/docs/agents/reference/mcp-configuration)を使います。

```json
{
  "servers": {
    "pdb2reaction": {
      "command": "pdb2reaction-mcp",
      "args": []
    }
  }
}
```

### カスタム Python MCP クライアント

```python
import asyncio

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

async def main():
    server_params = StdioServerParameters(command="pdb2reaction-mcp")
    async with stdio_client(server_params) as (read, write):
        async with ClientSession(read, write) as session:
            await session.initialize()
            result = await session.call_tool(
                "optimize_geometry",
                arguments={
                    "input_pdb": "r.pdb",
                    "charge": -1,
                    "max_cycles": 50,
                },
            )
            print(result.content)

asyncio.run(main())
```

## サンドボックス / 安全性に関する注意

- サーバーは呼び出した側の PATH、conda 環境、CUDA の設定を引き継ぎます。opt・tsopt・irc のように時間のかかる呼び出しには `timeout_seconds` を設定し、止まらない計算を打ち切ってください（デフォルトは時間制限なし）。
- ステージランナーとスキャン / 経路 / パイプラインのツールの出力は `out_dir` の下に置かれます。指定しないときは呼び出しごとに別の一時ディレクトリ `p2r_mcp_<subcmd>_…` を使うので、同時の呼び出しがぶつかりません。
- 構造 / I/O ヘルパーは `out_dir` を持たず、指定した出力パスに書きます。`extra_args` で CLI のフラグを追加できますが、型付きの出力パス、`--out-dir`、`--out-json/--no-out-json` は上書きできません。コマンドに渡したパスは、返された `argv` ですべて確かめられます。
- サーバーは `~/.bashrc` やログイン環境を変えず、ソフトウェアのインストールやモデルの配布元へのログインもしません。MLIP の重みと入力の PDB は、前もってディスクに置いてください。

## 関連ドキュメント

* [JSON 出力の一覧](json-output.md) — ツールが返す状態の欄と `summary.json`
* [トラブルシューティング](troubleshooting.md) — 実行に失敗したときの対処
* [コマンドの一覧（英語のみ）](../reference/commands/index.md) — 各ツールの元の CLI オプション（`extra_args` 用）
