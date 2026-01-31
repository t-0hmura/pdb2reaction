# トラブルシューティング

このページでは、`pdb2reaction` でよく遭遇するエラーと対処法をまとめます。  
コンソールに出てきたメッセージをそのまま検索（ページ内検索）すると見つけやすいように書いています。

---

## 実行前チェックリスト

長い計算を回す前に、最低限つぎを確認してください。

- `pdb2reaction -h` でヘルプが表示される
- UMA のモデルがダウンロードできる（Hugging Face のログイン/トークンが利用可能）
- 酵素系ワークフローでは、入力 PDB に **水素** と **元素記号（element column）** が入っている
- 複数の PDB を与える場合、**同じ原子が同じ順序** で並んでいる（座標だけが異なる）

---

## 入力 / 抽出の問題

### 「Element symbols are missing … add-elem-info を実行して下さい」
典型的なメッセージ:

```text
Element symbols are missing in '...'.
Please run `pdb2reaction add-elem-info -i ...` to populate element columns before running extract.
```

対処:
- 次を実行して element 列を補完します。

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
[multi] Atom count mismatch between input #1 and input #2: ...
[multi] Atom order mismatch between input #1 and input #2.
```

対処:
- **すべて** の構造を同じ前処理ワークフロー（同じプロトン化ツール、同じ設定）で作り直します。
- 水素付加を行う場合、全フレームで同一順序になる手順を選びます。

ヒント:
- MD 由来なら、同一のトポロジ/軌跡からフレーム抽出する方が安全です（異なるツールで生成した PDB を混ぜると順序がズレやすい）。

---

### 「ポケットが空っぽ / 必要な残基が落ちる」
症状:
- 抽出ポケットが想定より小さい
- 触媒残基が含まれない

対処の例:
- `--radius` を増やす（例: 2.6 → 3.5 Å）
- `--selected-resn` で残基を強制包含する（例: `--selected-resn 'A:123,B:456'`）
- 主鎖削除が強すぎる場合は `--exclude-backbone False` を試す

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
  pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
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

- HPC では、compute node から HF キャッシュ（ホームディレクトリ等）が書き込み可能か確認してください。

---

### CUDA / PyTorch の不整合
症状:
- GPU があるのに `torch.cuda.is_available()` が False
- import 時に CUDA runtime error が出る

対処:
- クラスタの CUDA と整合する PyTorch を入れます。
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

### 図の export が失敗する（Chrome がない）
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
- 最適化後も複数の虚振動数が残る

対処の例:
- オプティマイザモードを切り替える: `--opt-mode light` (Dimer) または `--opt-mode heavy` (RS-I-RFO)
- 余分な虚モードのフラット化を有効にする: `--flatten-imag-mode True`
- 最大サイクル数を増やす: `--tsopt-max-cycles 20000`
- より厳しい収束条件を使う: `--thresh baker` または `--thresh gau_tight`

---

### IRCが正常に終了しない

症状:
- IRCが明確な極小に到達する前に停止
- エネルギーが振動したり勾配が大きいまま

対処の例:
- ステップサイズを減らす: `--step-size 0.05`（デフォルトは0.10）
- 最大サイクル数を増やす: `--max-cycles 200`
- IRC実行前にTSに虚振動数が1つだけあることを確認

---

### MEP 探索（GSM/DMF）が失敗または予期しない結果

症状:
- 経路探索が有効なMEPなしで終了
- 結合変化が正しく検出されない

対処の例:
- `--max-nodes` を増やす（複雑な反応には15や20など）
- 端点の事前最適化を有効にする: `--preopt True`
- 別のMEP手法を試す: `--mep-mode dmf`（GSMが失敗した場合）またはその逆
- YAMLで結合検出パラメータを調整（`bond.bond_factor`、`bond.delta_fraction`）

---

## パフォーマンス / 安定性のヒント

- **VRAM不足**: `--radius` を減らす、`--max-nodes` を減らす、軽い最適化設定にする（`--opt-mode light`）
- **解析ヘシアンが遅いまたはOOM**: デフォルトの `FiniteDifference` を維持。`--hessian-calc-mode Analytical` は十分なVRAMがある場合のみ使用（500原子以上には16 GB以上推奨）
- **workers > 1**: HPC で UMA スループットは改善しますが、解析ヘシアンは無効になります
- **大規模系（1000原子以上）**: より小さなポケット（`--radius 2.5`）を抽出するか、マルチ GPUセットアップで実行を検討

---

## 不具合報告のときに添えると助かる情報

- 実行したコマンド（コピペ可能な形）
- `summary.log`（またはコンソール出力）
- 再現する最小入力（可能なら）
- OS / Python / CUDA / PyTorch バージョン
