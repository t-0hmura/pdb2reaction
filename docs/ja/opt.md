# `opt`（構造最適化）

`opt` サブコマンドは、1 つの構造を局所極小点へ最適化します。

---

## 主な用途

* **R・P・中間体の準備**: 経路探索や振動解析の前に、反応物・生成物・中間体の構造を緩和し、[`freq`](freq.md) で極小点（n_imag = 0）であることを確かめる
* **距離を保った緩和**: 選んだ原子の組の距離を保ったまま、ほかの自由度を緩和する
* **IRC の端点から R と P へ**: [`irc`](irc.md) の端点を、それぞれがつながる極小点まで最適化する

計算バックエンドにはデフォルトの **UMA**（Meta）のほか、`-b/--backend` で **ORB**、**MACE**、**AIMNet2**、**DFT** も選べます。

---

## 基本的な実行例

### 1. 標準の最小化

電荷とスピン多重度を明示し、`--out-json` で結果の要約も書き出します。

```bash
pdb2reaction opt -i input.pdb -q 0 -m 1 --out-json --out-dir ./result_opt
```

端末に `[opt] Converged!` が出て、`result_opt/result.json` の `"optimization_status"` が `"converged"` であれば収束しています。

### 2. 厳しい収束条件と軌跡の保存

収束条件を `gau_tight` にし、最適化の軌跡を残します。

```bash
pdb2reaction opt -i input.pdb -q 0 -m 1 --thresh gau_tight --dump \
    --out-dir ./result_opt_tight
```

### 3. 距離拘束

弱い調和拘束（20 eV·Å⁻²）で、原子 1 と 5 の距離を 2.0 Å へ近づけます。

```bash
pdb2reaction opt -i input.pdb -q 0 -m 1 \
    --distance-restraint '[(1,5,2.0)]' --restraint-k 20.0 --out-dir ./result_opt_rest
```

### 4. RFO

`--opt-mode hess` で、厳密な Hessian から始める RFO に切り替えます。

```bash
pdb2reaction opt -i input.pdb -q 0 -m 1 --opt-mode hess --out-dir ./result_opt_hess
```

---

## 処理の仕組みと計算仕様

1. **構造の読み込みと境界の凍結**: {ref}`電荷 <ja-charge-specification>`は `-q` または `-l` から決まります。`--freeze-links`（デフォルト有効）では、切り出したクラスターの{ref}`キャップ水素 <ja-link-hydrogen-and-frozen-atoms>`の親原子を凍結します。`--freeze-atoms` でほかの原子も凍結できます。
2. **最適化法の選択**（`--opt-mode`）: `grad`（別名 `lbfgs`）は勾配だけを使う **L-BFGS** を実行します。`hess`（別名 `rfo`）は **RFO** を実行し、厳密な Hessian から始めて [TS-BFGS](glossary.md#最適化アルゴリズム) 式で更新し、500 サイクルごとに計算し直します。更新式は YAML の `rfo.hessian_update` で変えられます。`tsopt` では同じ指定が{ref}`別の方法を選びます <ja-opt-mode-semantics>`。
3. **距離拘束の追加**: `--distance-restraint` の `(i, j, target)` のそれぞれが、力の定数 `--restraint-k`（eV·Å⁻²）の調和項を加え、原子 i と j の距離を `target`（Å）へ引き寄せます。`(i, j)` は最初の距離を保ちます。番号は 1 始まりで、`--zero-based` を付けると 0 始まりになります。
4. **最小化**: 収束条件を満たすか `--max-cycles` に達するまで構造を動かします。デフォルトの `--thresh gau` は、力の最大値が 4.5 × 10⁻⁴、RMS が 3.0 × 10⁻⁴ hartree/bohr 未満、ステップの最大値が 1.8 × 10⁻³、RMS が 1.2 × 10⁻³ bohr 未満を求め、Gaussian の既定と同じ条件です。
5. **虚振動の除去（`--flatten`）**: 最適化の後に Hessian を計算し、すべての虚振動モード（ν < −5.00 cm⁻¹）に沿って構造を 0.10 Å ずらして最適化し直します。虚振動が無くなるか 50 回に達するまで繰り返します。`--flatten` では、各回の後に端末の `[Imaginary modes] n=…` の行に n_imag が出て、最後の回の後にも虚振動が残ると `[flatten] WARNING: Remaining imaginary modes after the flatten loop: N` が出ます。

---

## 収束の判定

実行の終わり方は、端末と `result.json`（`--out-json`）に出ます。

| 終わり方 | `optimization_status` | 端末の行 | `scientific_status` / 終了コード |
| --- | --- | --- | --- |
| 収束 | `converged` | `[opt] Converged!` | `success` / 0 |
| `--max-cycles` に達して未収束 | `not_converged` | `[opt] Reached max cycles (N/M).` | `failed` / 1 |
| エネルギーが変わらなくなって停止（`--stop-plateau`） | `stalled` | `[opt] Stalled (energy plateau; not converged)` | `failed` / 1 |

どの行の後にも `[opt] Total cycles: N` が出ます。`not_converged` や `stalled` のときに何を変えるかは、{ref}`max_cycles とプラトー停止 <ja-troubleshooting-max-cycles>` を参照してください。

収束して得られるのは停留点で、極小点とは限りません。`opt` は `--flatten` のとき以外は最後の Hessian を計算しないので、final geometry に [`freq`](freq.md) を実行し、n_imag = 0 を確かめてください。

---

## 主な出力ファイル

実行が終わると、`--out-dir` に次のファイルができます。

```text
result_opt/
├─ final_geometry.xyz      # final geometry（常に出力）
├─ final_geometry.pdb      # 同じ構造の PDB（PDB/mmCIF 入力。Gaussian 入力では .gjf）
├─ optimization_trj.xyz    # 最適化の軌跡（--dump）
├─ optimization.pdb        # 同じ軌跡の PDB（--dump、PDB/mmCIF 入力）
├─ restart_NNN.yaml        # オプティマイザの状態（--dump と YAML の opt.dump_restart）
└─ result.json             # 結果の要約（--out-json）
```

{ref}`mmCIF の入力 <ja-mmcif-input>`と、PDB の欄に入りきらない大きな PDB の入力では、元の識別子を保った `.cif` も書きます。

* **final geometry**: `final_geometry.*` が最適化した構造です。[`freq`](freq.md) や経路探索に渡してください。
* **要約**: `--out-json` を付けると、[`result.json`](json-output.md) に `optimization_status`、最後のエネルギー `energy_hartree`（拘束のエネルギーを除いた値）、サイクル数 `n_opt_cycles` が記録されます。
* **端末**: サイクルごとの表と実行時間が出ます。`-v 3` では、実際に使った `geom`・`calc`・`opt`・`lbfgs` / `rfo` の設定も出ます。

---

## 主な CLI オプション

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス | （必須） | 入力構造ファイル（`.pdb`, `.cif`, `.mmcif`, `.xyz`, `.gjf`）。軌跡は 1 フレームを `.xyz` に切り出してから指定（{ref}`軌跡から 1 フレームを取り出す <ja-trajectory-one-frame>` を参照） |
| `-q, --charge` | 整数 | `None` | 系全体の総電荷。`-l` を使う場合と `.gjf` 入力のほかは必須 |
| `-l, --ligand-charge` | 文字列 | `None` | リガンドの総電荷（例: `-1`）または残基名ごとの電荷（例: `'GPP:-3,SAM:1'`）。`-q` を省いたときに使用（PDB/mmCIF 入力または `--ref-pdb`） |
| `-m, --multiplicity` | 整数 | `1` | スピン多重度（2S+1） |
| `--ref-pdb` | パス | `None` | `.xyz` / `.gjf` 入力に使う PDB/mmCIF のトポロジー。座標は `-i` から取る（例: IRC の端点。[irc](irc.md) を参照） |
| `-b, --backend` | 文字列 | `uma` | バックエンド（`uma`, `orb`, `mace`, `aimnet2`, `dft`） |
| `--opt-mode` | `grad` / `hess` | `grad` | 最適化法: L-BFGS / RFO（`lbfgs` と `rfo` は別名） |
| `--coord-type` | `cart` / `redund` / `dlc` / `tric` | `cart` | 最適化に使う座標系：デカルト座標 / 冗長内部座標 / 非局在化内部座標（DLC）/ 並進・回転を含む内部座標（TRIC） |
| `--thresh` | プリセット | `gau` | 収束条件（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） |
| `--max-cycles` | 整数 | `100000` | 最適化サイクルの上限。`--flatten` の各回と共有 |
| `--dump/--no-dump` | フラグ | `False` | 最適化の軌跡 `optimization_trj.xyz` を書き出す |
| `--distance-restraint` | 文字列 | `None` | 調和の距離拘束。直接書く（`'[(i,j,target_Å),...]'`）か、同じ項目を `constraints:` に並べた YAML/JSON ファイルで指定。`(i,j)` は最初の距離を保つ。原子は `'SAM,320,CS1'` のような[原子セレクタ](cli-conventions.md#原子セレクタ)でも書ける |
| `--restraint-k` | 実数 | `300` | 距離拘束の力の定数（eV·Å⁻²） |
| `--one-based/--zero-based` | フラグ | `--one-based` | `--distance-restraint` の番号を 1 から数えるか 0 から数えるか |
| `--freeze-links/--no-freeze-links` | フラグ | `True` | キャップ水素の親原子を凍結（PDB/mmCIF 入力または `--ref-pdb`） |
| `--freeze-atoms` | 文字列 | `None` | 凍結する原子（1 始まり、カンマ区切り: 例 `'1,3,5'`） |
| `--flatten/--no-flatten` | フラグ | `False` | 最適化の後に虚振動を除く |
| `--reject-uphill/--no-reject-uphill` | フラグ | `False` | `hess` で、エネルギーが 1e-4 hartree を超えて上がる RFO のステップを捨て、信頼半径を縮める |
| `--stop-plateau/--no-stop-plateau` | フラグ | `False` | エネルギーが変わらなくなったら（直近 50 サイクルの幅が 1e-4 hartree 未満）止め、`stalled` と報告 |
| `-o, --out-dir` | パス | `./result_opt/` | 出力先ディレクトリ |

全オプションの一覧は [自動生成 CLI リファレンス](../reference/commands/opt.md) を参照してください。

> **補足:** YAML（`--config`）では、`geom.freeze_atoms` で凍結する原子（1 始まり）を足せます。足した原子は `--freeze-links` と `--freeze-atoms` の原子と合わせて凍結されます。キーの一覧は YAML リファレンスの [`geom`](yaml-reference.md#geom)、[`opt`](yaml-reference.md#opt)、[`lbfgs`](yaml-reference.md#lbfgs)、[`rfo`](yaml-reference.md#rfo) にあります。

---

## 使用上の注意点

* **プラトーでの停止**: `--stop-plateau` は、力のノイズで力の収束条件に届かないときにサイクルを節約できますが、エネルギーが平坦であることは停留点の証拠になりません。実質的な上限は `--max-cycles` です。エネルギーの幅とサイクル数は `--stop-plateau-thresh` と `--stop-plateau-window` で指定できます。
* **凍結原子があるときの剛体運動**: Cartesian 座標での RFO の曲率の確認と `--flatten` は、剛体運動を [`freq`](freq.md#凍結境界での剛体モード) と同じように扱います。L-BFGS には影響しません。
* **凍結原子と拘束の全体**: クラスターモデルで凍結する原子や拘束の選び方は、{ref}`原子の固定と距離の拘束 <ja-freeze-atoms-and-restraints>` を参照してください。
* **オプティマイザの状態の書き出し**: `--dump` を付け、YAML の `opt.dump_restart` に正の整数 N を指定すると、N サイクルごとに `restart_NNN.yaml` を書きます。pdb2reaction はこのファイルを読み戻さないので、止まった計算は final geometry から `opt` をやり直してください。
* **モデルと精度**: `--backend-model` でバックエンドのモデルを、`--precision` で精度（`fp32`・`fp64`）を選べます。詳しくは自動生成 CLI リファレンスを参照してください。

---

## 関連ドキュメント

* [freq](freq.md) — 最適化した構造が極小点（n_imag = 0）かの確認
* [tsopt](tsopt.md) — 極小点ではなく TS（鞍点）の最適化
* [irc](irc.md) — TS から反応経路をたどり、最適化する端点を得る
* [extract](extract.md) — 最適化の前に活性部位モデルを切り出す
* [all](all.md) — IRC の端点の最適化まで含む一連のワークフロー
* [トラブルシューティング](troubleshooting.md) — 実行が失敗したときの切り分け
* [YAML リファレンス](yaml-reference.md) — `opt`、`lbfgs`、`rfo` のすべての設定
* [用語集](glossary.md) — L-BFGS、RFO などの用語
* {ref}`終了コード <ja-exit-codes>` — 終了ステータスの意味
