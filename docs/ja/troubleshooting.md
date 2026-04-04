# トラブルシューティング

このページでは、`pdb2reaction` でよく遭遇するエラーと対処法をまとめます。 コンソールに出てきたメッセージをそのまま検索（ページ内検索）すると見つけやすいように書いています。
症状から先に当たりを付けたい場合は、先に [典型エラー別レシピ](recipes-common-errors.md) を見てからこのページに戻ってください。

---

## 実行前チェックリスト

長い計算を回す前に、最低限次を確認してください。

- `pdb2reaction -h` でヘルプが表示される
- UMA のモデルがダウンロードできる（Hugging Face のログイン/トークンが利用可能）
- 酵素系ワークフローでは、入力 PDB に **水素** と **元素記号（element column）** が入っている
- 複数の PDB を与える場合、**同じ原子が同じ順序** で並んでいる（座標だけが異なる）

---

## 入力 / 抽出の問題

### 「Element symbols are missing … add-elem-info を実行してください」
典型的なメッセージ:

```text
Element symbols are missing in '...'.
Please run `pdb2reaction add-elem-info -i...` to populate element columns before running extract.
```

対処:
- 次を実行して element 列（元素記号列）を補完します。

 ```bash
 pdb2reaction add-elem-info -i input.pdb -o input_with_elem.pdb
 ```

- その後、`extract` / `all` を補完後の PDB で再実行します。

原因:
- PDB の element 列が空だったり不統一だったりすることが多く、`extract` は原子タイプ判定のために元素記号を必要とします。

---

### 「[multi] Atom count mismatch …」「[multi] Atom order mismatch …」
典型的なメッセージ:

```text
[multi] Atom count mismatch between input #1 and input #2:...
[multi] Atom order mismatch between input #1 and input #2.
```

対処:
- **すべて** の構造を同じ前処理ワークフロー（同じプロトン化ツール、同じ設定）で作り直します。
- 水素付加を行う場合、全フレームで同一順序になる手順を選びます。

ヒント:
- MD 由来なら、同一のトポロジー/軌跡からフレーム抽出する方が安全です（異なるツールで生成した PDB を混ぜると順序がズレやすい）。

---

### 「活性部位モデル（バインディングポケット）が空っぽ / 必要な残基が落ちる」
症状:
- 抽出された活性部位モデルが想定より小さい
- 触媒残基が含まれない

対処の例:
- `--radius` を増やしてください（例: 2.6 → 3.5 Å）
- `--selected-resn` で残基を強制包含してください（例: `--selected-resn 'A:123,B:456'`）
- 主鎖削除が強すぎる場合は `--no-exclude-backbone` を試してください

---

### 非標準残基が正しく切断されない

抽出された活性部位モデルに非標準の3文字コードを持つ修飾アミノ酸残基（リン酸化セリン、メチル化リシンなど）が含まれている場合、デフォルトでは主鎖切断やリンク水素付加が適用されません。`--modified-residue` で登録してください:

```bash
pdb2reaction extract -i complex.pdb -c PRE --modified-residue "SEP,TPO,MLY" -o pocket.pdb
```

`--modified-residue` で対応できない場合（残基の主鎖トポロジーが特殊な場合など）は、活性部位モデルを手動で構築し、下流のコマンド（`opt`、`tsopt`、`path-opt` など）に直接渡してください。

---

## 電荷 / スピンの問題

### 「電荷が必須」系のエラー（非 GJF 入力）
`.gjf` でない入力では、複数ステージで総電荷が必要になります。`-q/--charge` を省略した場合、PDB なら `--ligand-charge` を使って推定しようとしますが、推定できないと停止します。

対処:
- 電荷と多重度を明示する:

 ```bash
 pdb2reaction path-search -i R.pdb P.pdb -q 0 -m 1
 ```

- あるいは（抽出ありの場合）残基名ごとの電荷マッピングを与える:

 ```bash
 pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3'
 ```

---

## インストール / 環境の問題

### UMA のダウンロード/認証エラー（Hugging Face）
症状:
- モデルをダウンロードできない、認証が必要、といったエラー。

対処:
- 環境/マシンごとに一度ログインします。

 ```bash
 huggingface-cli login
 ```

- HPC では、計算ノードから HF キャッシュ（ホームディレクトリ等）に書き込み可能か確認してください。

---

### CUDA / PyTorch の不整合
症状:
- GPU があるのに `torch.cuda.is_available()` が False
- import 時に CUDA runtime error が出る

対処:
- クラスターの CUDA と整合する PyTorch を入れます。
- GPU が見えているか確認します:

 ```bash
 nvidia-smi
 python -c "import torch; print(torch.version.cuda, torch.cuda.is_available())"
 ```

---

### DMF モードが動かない（cyipopt がない）
DMF（`--mep-mode dmf`）を使うときに IPOPT/cyipopt の import エラーが出る場合:

対処:
- `pdb2reaction` を入れる前に conda-forge から `cyipopt` を入れるのが簡単です。

 ```bash
 conda install -c conda-forge cyipopt
 ```

---

### 図のエクスポートが失敗する（Chrome がない）
Plotly/Chrome 系のエラーで静的画像が出ない場合:

対処:
- headless Chrome を一度入れます。

 ```bash
 plotly_get_chrome -y
 ```

---

## 計算 / 収束の問題

### TS 最適化が収束しない

症状:
- TS 最適化が多くのサイクルを回しても収束しない
- 最適化後もヘシアン行列に複数の負の固有値が残る（虚振動数が 2 本以上）

対処の例:
- オプティマイザーモードを切り替えてください: `--opt-mode grad`（Dimer 法）または `--opt-mode hess`（RS-I-RFO 法、デフォルト）
- 余分な虚振動数モードのフラット化を有効にしてください: `--flatten`
- 最大サイクル数を増やしてください: `--max-cycles 20000`（単独の `tsopt` の場合）、`--tsopt-max-cycles 20000`（`all` の場合）
- より厳しい収束閾値を使ってください: `--thresh baker` または `--thresh gau_tight`

---

### IRC が正常に終了しない

症状:
- IRC が明確な極小構造（停留点）に到達する前に停止する
- エネルギーが振動したり勾配ノルムが大きいままになる

対処の例:
- ステップサイズを減らしてください: `--step-size 0.05`（デフォルト: 0.10 Bohr sqrt(amu)）
- 最大サイクル数を増やしてください: `--max-cycles 200`
- IRC 実行前に TS 候補の虚振動数が 1 本（|ν| >= 100 cm⁻¹）だけであることを確認してください

---

### MEP 探索（GSM/DMF）が失敗または予期しない結果

症状:
- 経路探索が有効な MEP なしで終了
- 結合変化が正しく検出されない

対処の例:
- `--max-nodes` を増やしてください（複雑な反応には 15 や 20 など）
- 端点の事前最適化を有効にしてください: `--preopt`
- 別の MEP 手法を試してください: `--mep-mode dmf`（GSM が失敗した場合）またはその逆
- YAML で結合検出パラメータを調整してください（`bond.bond_factor`、`bond.delta_fraction`）

---

## パフォーマンス / 安定性のヒント

- **VRAM 不足**: `--radius` の値を減らして活性部位モデルを小さくする、`--max-nodes` を減らす、軽い最適化設定にする（`--opt-mode grad`）
- **解析ヘシアン（Analytical Hessian）が遅いまたは OOM**: デフォルトの `FiniteDifference` を維持してください。`--hessian-calc-mode Analytical` は十分な VRAM がある場合のみ使用してください（500 原子以上では 16 GB 以上推奨）
- **workers > 1**: HPC で UMA のスループットは改善しますが、解析ヘシアンは無効になります
- **大規模系（1000 原子以上）**: より小さな活性部位モデル（`--radius 2.5` Å）を抽出するか、マルチ GPU セットアップでの実行を検討してください

---

## バックエンド選択ガイド

ベンチマーク: LBFGS 構造最適化、29-177 原子クラスターモデル、NVIDIA RTX 5080 (16 GB VRAM)。

| バックエンド | 精度 | 速度 (中央値 s/step) | VRAM | 備考 |
|------------|------|---------------------|------|------|
| **UMA-s1p1** | 良好 | 0.03 s | ~2 GB | デフォルト。高速、探索向き |
| **UMA-s1p2** | より高精度 | 0.08 s | ~4 GB | 2-3 倍遅い |
| **UMA-m1p1** | より高精度 | 0.22 s | ~8 GB | 中規模モデル、VRAM 大 |
| **MACE** | 最高精度 | 0.37 s | ~4 GB | 最高精度だが別環境が必要（e3nn 競合） |
| **ORB** | 変動あり | 0.02 s | ~2 GB | 最速だが複雑な反応で失敗率高い |

**推奨:** UMA-s1p1 で高速スクリーニング → MACE/UMA-s1p2 で主要結果を検証。

## GPU メモリ (VRAM) 目安

| 原子数 | LBFGS 最適化 | ヘシアン（解析的） | ヘシアン（有限差分） |
|-------|------------|-----------------|------------------|
| 50 | ~2 GB | ~3 GB | ~2 GB |
| 100 | ~3 GB | ~6 GB | ~3 GB |
| 200 | ~4 GB | ~12 GB | ~4 GB |
| 500 | ~6 GB | 16 GB で OOM | ~6 GB |

`torch.cuda.OutOfMemoryError` の場合: `--hessian-calc-mode FiniteDifference`、`--radius` の縮小、小さいモデル（`uma-s-1p1`）を検討してください。

## 不具合報告のときに添えると助かる情報

- 実行したコマンド（コピペ可能な形）
- `summary.log`（またはコンソール出力）
- 再現する最小入力（可能なら）
- OS / Python / CUDA / PyTorch バージョン
