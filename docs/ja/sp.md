# `sp`（一点計算）

## 概要

`sp` サブコマンドは、選んだバックエンドで 1 つの構造の**エネルギーと原子に働く力**を計算し、`--hess` を付けると **Hessian** も計算します。構造最適化は行わず、入力の構造のまま評価します。

### 主な用途

* **最適化の前の確認**: 電荷と多重度が受け付けられ、バックエンドが有限のエネルギーと力を返すかを確かめる
* **バックエンドの比較**: 同じ構造を、[MLIP](backends.md)（機械学習原子間ポテンシャル）の UMA・ORB・MACE・AIMNet2 か、DFT（`-b dft`）で評価する
* **参照値の作成**: 力と Hessian を `.npy` ファイルとして、エネルギーを端末か `result.json` から得て、自分の解析に使う

---

## 基本的な実行例

### 1. エネルギーと力

デフォルトのバックエンド（UMA）で、中性の一重項を評価します。

```bash
pdb2reaction sp -i structure.pdb -q 0 -m 1 --out-json
```

端末に `[sp] energy = … a.u.  |force|_max = … a.u./bohr` が出て、`result_sp/` に `forces.npy` と、`energy_au` を持つ `result.json` があれば成功です。

### 2. Hessian も計算する

`--hess` を付けると Hessian も計算します。

```bash
pdb2reaction sp -i structure.pdb -q 0 -m 1 --hess
```

---

## 処理の仕組みと計算仕様

1. **構造の読み込み**:
PDB・mmCIF・XYZ・GJF を読み込みます。電荷は `-q`、`-l`（PDB/mmCIF 入力）、YAML の `calc.charge`、`.gjf` のヘッダーのいずれかから決まります。`--freeze-atoms` で指定した原子は凍結します。
2. **エネルギーと力**:
入力の構造でバックエンドを 1 回呼び、エネルギーと力の最大成分を端末に表示して、力を `forces.npy` に保存します。
3. **Hessian（`--hess` 指定時）**:
`--hessian-calc-mode FiniteDifference` は力を数値微分し、`Analytical` は UMA・ORB・MACE・AIMNet2・DFT の解析 Hessian を使います。UMA で `--uma-workers` を 2 以上にすると `Analytical` は{ref}`使えません <ja-workers-analytical-error>`。

---

## 主な出力ファイル

`--out-dir` に以下のファイルを書き出します。

| ファイル | 内容 | 書き出す条件 |
| --- | --- | --- |
| `forces.npy` | 力の `(N, 3)` 配列（Hartree/bohr） | 常に |
| `hessian.npy` | 質量重み付けなしの Cartesian Hessian（Hartree/bohr²）。`(3N, 3N)`、凍結原子があるときは動ける M 原子（入力の順）の `(3M, 3M)` | `--hess` 指定時 |
| `result.json` | エネルギー（`energy_au`）、バックエンド、モデル、電荷、多重度、原子数、`.npy` ファイルのパス、経過時間 | `--out-json` 指定時 |
| `summary.json` | `result.json` の写し。`result.json` を読む | `--out-json` 指定時 |

---

## 主な CLI オプション

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス | （必須） | 入力構造ファイル（`.pdb`, `.cif`, `.xyz`, `.gjf` 等） |
| `-q, --charge` | 整数 | `None` | 系全体の総電荷。`-l`、YAML の `calc.charge`、`.gjf` 入力のどれも無ければ必須 |
| `-m, --multiplicity` | 整数 | `1` | スピン多重度（2S+1）。`.gjf` 入力ではファイルの値を使用 |
| `-l, --ligand-charge` | 文字列 | `None` | 残基ごとの形式電荷（例: `'SAM:1,GPP:-3'`）またはリガンドの総電荷。PDB/mmCIF 入力が必要 |
| `-b, --backend` | 文字列 | `uma` | 計算バックエンド（`uma`, `orb`, `mace`, `aimnet2`, `dft`）。`-b dft` の設定は [MLIP の TS を DFT で確かめる](dft-backend.md) を参照 |
| `--hess/--no-hess` | フラグ | `False` | Hessian も計算して `hessian.npy` に書き出す |
| `--hessian-calc-mode` | `FiniteDifference` / `Analytical` | `FiniteDifference` | Hessian の計算法（有限差分 / 解析的）。`--hess` と併用 |
| `--freeze-atoms` | 文字列 | `None` | 凍結する原子インデックス（1 始まり、カンマ区切り: 例 `'1,3,5'`） |
| `-o, --out-dir` | パス | `./result_sp/` | 出力先ディレクトリ |
| `--out-json/--no-out-json` | フラグ | `False` | `result.json` と `summary.json` を出力 |

全オプションの一覧は [自動生成 CLI リファレンス](../reference/commands/sp.md) を参照してください。

> **補足:** YAML（`--config`）では、`calc` でバックエンドを設定し、`geom.freeze_atoms`（1 始まり）で `--freeze-atoms` に凍結原子を追加できます。

---

## 使用上の注意点

* **エネルギーがおかしいとき**: {ref}`電荷と多重度 <ja-charge-spin-problems>`を見直してください。
* **凍結原子**に働く力は 0 になります。
* **キャップ水素**: `sp` は `extract` が付けたキャップ水素の親原子を自動では凍結しません。固定したい場合は `--freeze-atoms` に指定してください。
* **原子電荷**: `sp -b dft` が出すのは DFT のエネルギーと力だけです。Mulliken・meta-Löwdin・IAO の電荷が必要なときは [`dft`](dft.md) を使ってください。
* **失敗したとき**: 1 行の `Error: …` か、トレースバック付きの `Unhandled error during single-point calculation:` が出て、0 以外の終了コードで終わります。[エラー処理](json-output.md#エラー処理)を参照してください。
* **終了コード**: {ref}`終了コード <ja-exit-codes>`を参照してください。

---

## 関連ドキュメント

* [opt](opt.md) — 構造最適化
* [tsopt](tsopt.md) — 遷移状態（TS）候補の構造最適化
* [freq](freq.md) — 振動解析と熱化学
* [dft](dft.md) — 原子電荷も出す DFT 一点計算
* [MLIP バックエンド](backends.md) — バックエンドの選び方と、UMA に要る Hugging Face へのログイン
* [MLIP の TS を DFT で確かめる](dft-backend.md) — `-b dft` の設定と GPU メモリ
* [トラブルシューティング](troubleshooting.md) — 実行に失敗したときの対処
