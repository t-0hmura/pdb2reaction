# `sp`

`pdb2reaction sp` は、単一の構造に対して MLIP のエネルギー + 原子間力（オプションで完全な Hessian）を評価します。用途は次のとおりです。

- 最適化の実行前に、構造のエネルギー / 原子間力 / Hessian を手早くサニティチェックする
- バックエンド同士を直接比較する
- optimizer のループ外で参照値 / Hessian を生成する

## 実行例

コマンド形式:

```bash
pdb2reaction sp -i FILE -q INT -m INT [-b uma|orb|mace|aimnet2] [--hess] [options]
```

エネルギー + 原子間力（UMA バックエンド、中性閉殻）:

```bash
# energy + forces (UMA backend, neutral closed-shell)
pdb2reaction sp -i structure.pdb -q 0 -m 1
```

完全な Hessian も計算（自動モードでは UMA は解析的、その他は有限差分）:

```bash
# also compute the full Hessian (automatic backend selection)
pdb2reaction sp -i structure.pdb -q 0 -m 1 --hess
```

## 出力

`sp` はデフォルトで出力を `result_sp/` 以下に書き出します。エネルギー（スカラー値）と `|force|_max` は stdout に出力され、`forces.npy`（`--hess` 指定時は `hessian.npy` も）は常にそこへ書き出されます。

| ファイル | 内容 | 書き出し |
|---|---|---|
| _stdout_ | エネルギー（スカラー値）(a.u.) と `|force|_max`。`[sp] energy = …` の形式で出力 | 常に |
| `forces.npy` | 原子単位 (Hartree / Bohr) の `(N, 3)` 力配列 | 常に |
| `hessian.npy` | 質量重みなしの `(3N, 3N)` Hessian (Hartree / Bohr²) | `--hess` 指定時のみ |
| `result.json` / `summary.json` | 機械可読なエネルギー (a.u.)、バックエンド、電荷/スピン、npy 出力へのパス、経過時間 | `--out-json` 指定時のみ |

`sp` は人間可読な `summary.log` を書き出しません。

### Hessian バックエンド

`--hess` を設定すると、バックエンドの選択が Hessian の計算戦略を決めます。

- `--backend uma`（デフォルト）→ UMA の torch autograd 経路による `Analytical` Hessian
- `--backend orb` / `mace` / `aimnet2` → 自動モードでは `FiniteDifference`

UMA、ORB、MACE、AIMNet2 はすべて解析 Hessian を実装しています。明示的に使う場合は `--hessian-calc-mode Analytical`、数値的なクロスチェックには `FiniteDifference` を指定します。UMA では `workers > 1` と明示的な解析 Hessian を併用できずエラーになるため、`workers = 1` または有限差分を使用してください。

## CLI オプション

フラグの完全な一覧は自動生成された [コマンドリファレンス](../reference/commands/index.md) にあります。下表では説明が必要なオプションを扱います。

| フラグ | デフォルト | 意味 |
|---|---|---|
| `-i, --input FILE` | — | PDB / XYZ / GJF の構造ファイル（必須） |
| `-q, --charge INT` | — | 系の総電荷（非 GJF では必須；GJF はテンプレートから継承） |
| `-l, --ligand-charge TEXT` | — | 残基別の電荷マッピング（例: `SAM:1,GPP:-3`）。`-q` の自動導出に使用 |
| `-m, --multiplicity INT` | `1` | スピン多重度、2S+1（任意；省略時は 1。GJF はテンプレートから継承） |
| `-b, --backend [uma\|orb\|mace\|aimnet2]` | `uma` | MLIP バックエンドの選択 |
| `--hess / --no-hess` | `--no-hess` | `hessian.npy` も計算して書き出す |
| `--hessian-calc-mode [Analytical\|FiniteDifference]` | auto | 特定の Hessian モードを強制（`--hess` 指定時のみ有効） |
| `-o, --out-dir PATH` | `./result_sp/` | 出力ディレクトリ |
| `--precision [fp32\|fp64]` | `fp32` | バックエンドに渡す数値精度 |
| `--config PATH` | — | `calc.*`, `geom.*` のデフォルトを与える YAML 設定 |
| `--out-json / --no-out-json` | `--no-out-json` | 機械可読な `result.json`（`summary.json` にもミラー）を出力ディレクトリに書き出す |
| `--show-config / --dry-run` | off | 実効的なマージ済み設定を出力 / 実行せずに検証 |

完全な一覧（workers、溶媒補正など）を見るには `pdb2reaction sp --help-advanced` を実行してください。

## 注記

- 一点 DFT（gpu4pyscf / PySCF）のベンチマークには、代わりに [`dft`](dft.md) を使用してください。

## 関連項目

- [`opt`](opt.md) — 構造を最適化する
- [`tsopt`](tsopt.md) — TS 候補を精密化する
- [`freq`](freq.md) — 熱化学を含む振動解析
- [`dft`](dft.md) — 一点 DFT の対応コマンド（PySCF / gpu4pyscf を使用）
