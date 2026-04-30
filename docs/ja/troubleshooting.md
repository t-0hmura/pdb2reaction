# トラブルシューティング

このページでは、`pdb2reaction` でよく遭遇するエラーと対処法をまとめます。
症状から先に当たりを付けたい場合は、先に [典型エラー別レシピ](recipes-common-errors.md) を見てからこのページに戻ってください。

## 実行前チェックリスト

長い計算を回す前に、最低限次を確認してください。

- UMA のモデルがダウンロードできる（Hugging Face のログイン/トークンが利用可能）
- 酵素系ワークフローでは、入力 PDB に **水素** と **元素記号（element column）** が入っている
- 複数の PDB を与える場合、**同じ原子が同じ順序** で並んでいる（座標だけが異なる）

---

(ts-input-extraction)=
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
- `--selected-resn` で残基を強制包含してください（例: `--selected-resn 'A:123,B:456'`）。残基 ID 仕様の詳細は CLI 規約の {ref}`ja-selected-resn-takes-ids` を参照。
- 以前に `--exclude-backbone` を明示的に有効化していて、その削除が強すぎる場合は当該フラグを外す（または `--no-exclude-backbone` を渡す）ことで主鎖原子を保持できます

---

### エネルギー・障壁の計算値が不正確

症状:
- 計算されたエネルギーや反応障壁が不合理に見える
- モデルサイズを大きくすると結果が大幅に変わる

対処:
- 抽出された活性部位モデルが小さすぎると、エネルギーや障壁の計算値が不正確になることがあります。抽出半径を大きくする（例: `-r 4.0` 以上）ことで、タンパク質環境をより多く含めて精度を改善できます:

  ```bash
  pdb2reaction extract -i complex.pdb -c 'SUB' -o model.pdb -r 4.0
  ```

---

### 非標準残基が正しく切断されない

抽出された活性部位モデルに非標準の3文字コードを持つ修飾アミノ酸残基（リン酸化セリン、メチル化リシンなど）が含まれている場合、デフォルトでは主鎖切断やリンク水素付加が適用されません。`--modified-residue` で登録してください:

```bash
pdb2reaction extract -i complex.pdb -c PRE --modified-residue "SEP,TPO,MLY" -o pocket.pdb
```

`--modified-residue` で対応できない場合（残基の主鎖トポロジーが特殊な場合など）は、活性部位モデルを手動で構築し、下流のコマンド（`opt`、`tsopt`、`path-opt` など）に直接渡してください。

---

(ts-charge-spin)=
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
  pdb2reaction extract -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3'
  ```

---

(ts-install-env)=
## インストール / 環境の問題

### UMA のダウンロード/認証エラー（Hugging Face）
症状:
- モデルをダウンロードできない、認証が必要、といったエラー。

対処:
- 環境/マシンごとに一度ログインします。

  ```bash
  hf auth login
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

(ts-calc-conv)=
## 計算 / 収束の問題

### 最適化が `max_cycles` に達し、`max(force)` が閾値をわずかに超える

症状:
- オプティマイザーが `max_cycles` 上限まで回りきり、最終サマリで `max(force)` や `rms(force)` が目標をわずかに上回る（例: `baker` の 3×10⁻⁴ au に対して 4×10⁻⁴ au）。
- 一方で、エネルギー自体は明らかにフラット化しており、10⁻⁵–10⁻⁴ au レベルで振動している。

原因:
- MLIP の勾配（力）計算には小さな確率的ノイズフロアがあり（UMA 系で典型的に ~4×10⁻⁴ au）、これが力ベースの収束閾値（`baker` = 3×10⁻⁴ au）を上回る場合があります。その結果、構造はすでに収束しているにもかかわらず、力閾値を満たすことができません。

対処:
- **エネルギープラトーによるフォールバック収束**（v0.3.5 新機能）が自動でこの状況を処理します: `opt.energy_plateau: true` のとき、直近 `opt.energy_plateau_window`（デフォルト 50）ステップのエネルギーレンジが `opt.energy_plateau_thresh`（デフォルト `1×10⁻⁴ au ≈ 0.06 kcal/mol`）を下回ると収束と判定されます。多くの場合、ユーザー側での対応は不要です。
- 自動フォールバックを上書きしたい場合は、力の閾値を手動で緩めてください: `--thresh gau`（`opt` のデフォルト）または `--thresh gau_loose`。
- `opt.energy_plateau_thresh` / `opt.energy_plateau_window` は YAML からチューニングでき、`opt.energy_plateau: false` で無効化できます。
- 注意: このプラトーフォールバックは **chain-of-states オプティマイザー**（`path-opt`、`path-search` の string/GSM/DMF 段階）では**スキップ**されます（単一のスカラーエネルギー履歴ではなく、イメージごとのエネルギー配列を保持しているため）。

---

### TS 最適化が収束しない

症状:
- TS 最適化が多くのサイクルを回しても収束しない
- 最適化後もヘシアン行列に複数の負の固有値が残る（虚振動数が 2 本以上）

対処の例（CLI フラグと YAML キーは補完的、必要に応じて併用してください）:
- オプティマイザーモードを切り替えてください: `--opt-mode grad`（Dimer 法）または `--opt-mode hess`（RS-I-RFO 法、デフォルト）
- 余分な虚振動数モードのフラット化を有効にしてください: `--flatten`（単独の `tsopt`、`opt`、および `pdb2reaction all` で利用可能。デフォルトは無効）
- 最大サイクル数を増やしてください: `--max-cycles 20000`（単独の `tsopt` の場合）、`--tsopt-max-cycles 20000`（`all` の場合）
- より厳しい収束閾値を使ってください: `--thresh baker` または `--thresh gau_tight`
- YAML でステップサイズ / 信頼半径を縮小してください — LBFGS/Dimer: `lbfgs.max_step` / `hessian_dimer.lbfgs.max_step`、RFO/RS-I-RFO: `rfo.trust_radius` / `rfo.trust_min` / `rfo.trust_max`（および `rsirfo` セクション）。セクション構成は [YAML リファレンス](yaml-reference.md) を参照

---

### IRC が正常に終了しない

症状:
- IRC が明確な極小構造（停留点）に到達する前に停止する
- エネルギーが振動したり勾配ノルムが大きいままになる

対処の例:
- ステップサイズを減らしてください: `--step-size 0.05`（デフォルト: 0.10 bohr、質量重み付けなしのデカルト座標）
- 最大サイクル数を増やしてください: `--max-cycles 200`
- IRC 実行前に TS 候補の虚振動数が 1 本（|ν| >= 100 cm⁻¹）だけであることを確認してください。5 cm⁻¹ 検出閾値と 100 cm⁻¹ 品質ゲートの違いは用語集 {ref}`ja-imaginary-mode-thresholds` を参照。

---

### MEP 探索（GSM/DMF）が失敗または予期しない結果

症状:
- 経路探索が有効な MEP なしで終了
- 結合変化が正しく検出されない

対処の例:
- `--max-nodes` をデフォルトの 20 より増やしてください（複雑な反応には 30 や 40 など）
- 端点の事前最適化を有効にしてください: `--preopt`
- 別の MEP 手法を試してください: `--mep-mode dmf`（GSM が失敗した場合）またはその逆
- YAML で結合検出パラメータを調整してください（`bond.bond_factor`、`bond.delta_fraction`）

## パフォーマンス / 安定性のヒント

- **VRAM 不足**: `--radius` の値を減らして活性部位モデルを小さくする、`--max-nodes` を減らす、軽い最適化設定にする（`--opt-mode grad`）
- **解析ヘシアン（Analytical Hessian）が遅いまたは OOM**: デフォルトの `FiniteDifference` を維持してください。`--hessian-calc-mode Analytical` は十分な VRAM がある場合のみ使用してください（500 原子以上では 16 GB 以上推奨）
- **workers > 1**: HPC で UMA の処理速度は改善しますが、解析ヘシアンは無効になります
- **大規模系（1000 原子以上）**: より小さな活性部位モデル（`--radius 2.5` Å）を抽出するか、マルチ GPU セットアップでの実行を検討してください

## バックエンド選択ガイド

以下は NVIDIA RTX 5080 (16 GB VRAM) で小〜中規模クラスターモデルを LBFGS 1 ステップ動かした場合の目安値です。バックエンド選択の参考に使う想定で、厳密なベンチマークではありません。

| バックエンド | 速度 (中央値 s/step) | VRAM | 備考 |
|------------|---------------------|------|------|
| **UMA-s-1.1** (デフォルト) | 0.03 s | ~2 GB | 高速、探索向き |
| **UMA-m-1.1** | 0.22 s | ~8 GB | 中規模モデル、VRAM 大 |
| **MACE-OMOL-0** | 0.37 s | ~4 GB | 別 conda 環境が必要（`e3nn` 競合） |
| **Orb-v3-omol** | 0.02 s | ~2 GB | 最速。注意点は下を参照 |

**推奨:**
- まず UMA-s-1.1 で高速スクリーニングし、その後 MACE または UMA-m-1.1 で主要結果を相互確認してください。
- SAM 依存 S~N~2 / メチル転移反応では、MACE と UMA-s-1.1 は相補的に働く傾向があるので、片方の TS が疑わしいときはもう片方で検証するのが有効です。
- Orb-v3-omol は反応座標自体は正しく捉えるものの、TS が虚振動を複数持つ状態に収束しやすく、clean な単一虚振動 TS は保証されません。Orb は機構回収の first-pass には向きますが、定量的な反応速度論や振動解析には UMA / MACE / DFT で再評価するか backend を切り替えてください。

## GPU メモリ (VRAM) 目安

| 原子数 | LBFGS 最適化 | ヘシアン（解析的） | ヘシアン（有限差分） |
|-------|------------|-----------------|------------------|
| 50 | ~2 GB | ~3 GB | ~2 GB |
| 100 | ~3 GB | ~6 GB | ~3 GB |
| 200 | ~4 GB | ~12 GB | ~4 GB |
| 500 | ~6 GB | 16 GB で OOM | ~6 GB |

`torch.cuda.OutOfMemoryError` の場合: `--hessian-calc-mode FiniteDifference`、`--radius` の縮小、小さいモデル（YAML設定で `calc.model: uma-s-1p1`）を検討してください。

## 不具合報告のときに添えると助かる情報

- 実行したコマンド（コピペ可能な形）
- `summary.log`（またはコンソール出力）
- 再現する最小入力（可能なら）
- OS / Python / CUDA / PyTorch バージョン
