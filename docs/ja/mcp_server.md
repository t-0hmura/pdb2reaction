# pdb2reaction MCP サーバー

`pdb2reaction-mcp`（エイリアス `p2r-mcp`）は、[MCP](https://modelcontextprotocol.io/)
サーバーであり、MCP に対応した任意のエージェントが stdio 上の JSON-RPC を介して
すべての `pdb2reaction` CLI サブコマンドを駆動できるようにします。Claude Desktop /
Claude Code / Cursor / Codeium のほか、公式の Python または TypeScript MCP SDK で
構築したカスタムエージェントからも利用できます。

## インストール

```bash
pip install "pdb2reaction[mcp]"
```

これにより `mcp[cli]` 依存関係が追加され、2 つのコンソールスクリプト
`pdb2reaction-mcp` および `p2r-mcp`（エイリアス）が登録されます。

## ツール

18 個のツールがあり、それぞれが CLI サブコマンドに 1 対 1 で対応します。各ツールは次のフィールドを持つ構造化された dict を返します。

- `schema_version`: エンベロープのバージョン。実際の値は `pdb2reaction.mcp._runner.MCP_SUBCMD_RESULT_SCHEMA_VERSION`。値が上がるとフィールドセットや値の型の変更を意味するため、本ドキュメント中のリテラルではなく定数に対して固定してください。
- `status`: `ok` | `failed` | `summary_missing` | `summary_parse_error`
- `exit_code`: サブプロセスの終了コード
- `out_dir`: CLI が書き込んだ作業ディレクトリ
- `summary`: パース済みの `summary.json`（CLI 出力スキーマ。ステージごとの形式は [JSON 出力リファレンス](json-output.md) を参照）
- `stderr_tail` / `stdout_tail`: プロセス出力の末尾約 60 行
- `hint`: CLI エラーメッセージから抽出した `; recover: <hint>` サフィックス（存在する場合）
- `argv`: 実行された完全な argv（再現性のため）

### 構造化されたエラーエンベロープ

サブコマンドが失敗した場合、パース済みの `summary`（または同階層の `result.json`）に拡張エラーエンベロープが含まれます。これにより、エージェントはテキストをパースせずに例外クラスの階層をパターンマッチできます。

- `error`: 元の例外の `str(exc)`
- `error_type`: 例外クラス名（例: `"OptimizationError"`）
- `error_class_chain`: MRO のクラス名（例: `["OptimizationError", "RuntimeError", "Exception", "BaseException"]`）
- `error_module`: 例外クラスが定義されているモジュール
- `error_label`: 上位レベルの CLI ステージラベル

### ステージランナー

| MCP ツール | CLI サブコマンド | 目的 |
|---|---|---|
| `optimize_geometry` | `pdb2reaction opt` | 単一の分子構造を最適化 |
| `find_transition_state` | `pdb2reaction tsopt` | TS 探索（RS-I-RFO / Dimer / TRIM / RS-P-RFO） |
| `run_irc` | `pdb2reaction irc` | TS 構造からの IRC 積分 |
| `compute_frequencies` | `pdb2reaction freq` | 振動解析 + 熱化学 |
| `run_single_point` | `pdb2reaction sp` | MLIP の一点エネルギー + 原子間力（+ オプションで Hessian） |

### スキャン / 経路 / パイプライン

| MCP ツール | CLI サブコマンド | 目的 |
|---|---|---|
| `scan_1d` / `scan_2d` / `scan_3d` | `pdb2reaction scan{,2d,3d}` | 拘束駆動の距離スキャン |
| `optimize_path` | `pdb2reaction path-opt` | 2 端点間の MEP 最適化 |
| `search_paths` | `pdb2reaction path-search` | 再帰的な反応経路探索 |
| `run_full_pipeline` | `pdb2reaction`（`all` サブコマンド） | エンドツーエンド: extract → MEP → TS → IRC → freq → DFT |
| `run_single_point_dft` | `pdb2reaction dft` | gpu4pyscf による一点 DFT |

### 構造 / I/O ヘルパー

| MCP ツール | CLI サブコマンド | 目的 |
|---|---|---|
| `extract_active_site` | `pdb2reaction extract` | リガンド周辺の球を切り出す |
| `add_element_info` | `pdb2reaction add-elem-info` | PDB の元素列を修復 |
| `fix_altloc` | `pdb2reaction fix-altloc` | PDB の代替位置（altloc）を解決 |
| `plot_trajectory` | `pdb2reaction trj2fig` | エネルギープロファイル図（デフォルトは PNG。SVG/PDF/HTML/CSV も可） |
| `plot_energy_diagram` | `pdb2reaction energy-diagram` | カテゴリ別エネルギーダイアグラム |
| `detect_bond_changes` | `pdb2reaction bond-summary` | 2 つの PDB 間の結合変化の差分 |

## オプトインの IRC 収束ガード

`run_irc`（CLI: `pdb2reaction irc`）は `irc_pos_def: bool` を受け付けます。これにより IRC
収束には、質量重み付き Hessian が正定値であることが追加で要求され、rms のみの
基準が局所極小に到達する前に成功と判定してしまう IRC の「ショルダー」での偽収束を
防ぎます。デフォルトは `None`（rms のみ、従来動作）です。

`find_transition_state`（CLI: `pdb2reaction tsopt`）は `--opt-mode` により
代替の TS オプティマイザも公開しています。

- `opt_mode="trim"` — Helgaker (1991) の trust-region image-minimization TS opt
- `opt_mode="rsprfo"` — Banerjee (1985) の restricted-step P-RFO TS opt

## クライアント設定

すべての MCP クライアントは同じ `mcpServers` スキーマを取ります。次のスニペットを
クライアントの MCP 設定ファイルに記述してください。

- Claude Desktop — `~/Library/Application Support/Claude/claude_desktop_config.json`（macOS） / `%APPDATA%\Claude\claude_desktop_config.json`（Windows）
- Cursor — `~/.cursor/mcp.json`
- その他のクライアント — 各クライアント自身の MCP サーバードキュメントを参照

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

明示的な環境変数の上書き（PATH / CUDA_VISIBLE_DEVICES）を含む完全な例は
[`examples/mcp_client_config.json`](../../examples/mcp_client_config.json)
を参照してください。

### カスタム Python MCP クライアント

```python
import subprocess
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

server_params = StdioServerParameters(command="pdb2reaction-mcp")
async with stdio_client(server_params) as (read, write):
    async with ClientSession(read, write) as session:
        await session.initialize()
        tools = await session.list_tools()
        result = await session.call_tool(
            "optimize_geometry",
            arguments={
                "input_pdb": "r.pdb",
                "charge": -1,
                "max_cycles": 50,
            },
        )
        print(result.content)
```

## サンドボックス / 安全性に関する注意

- MCP サーバーは呼び出し元の環境の PATH、conda 環境、CUDA セットアップを継承します。
  長時間実行されるツール（opt / tsopt / irc）は `pdb2reaction` CLI をサブプロセスで
  起動するため、エージェントは各ツール呼び出しで `timeout_seconds` を設定し、
  暴走する計算を制限してください。
- 出力ファイルは `out_dir` キーワード引数の配下に置かれます（デフォルトは一意の
  `tempfile.mkdtemp(prefix="p2r_mcp_<subcmd>_…")` で、同時並行のエージェント呼び出しが
  衝突しないようになっています）。
- サーバーは `~/.bashrc` / ログイン環境を変更したり、ソフトウェアをインストールしたり、
  `out_dir` の外に書き込んだりはしません。すべての MLIP の重みや PDB 入力は
  あらかじめディスク上に存在している必要があります。
