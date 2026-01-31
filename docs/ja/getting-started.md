# はじめに

## 概要

`pdb2reaction` は、機械学習原子間ポテンシャル (MLIP) を使用して **PDB構造** を **酵素反応経路** に変換するPython CLIツールキットです。

基本的には、以下のような**単一コマンド**で酵素反応経路をモデリングできます:
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
```

---
さらに、`--tsopt True --thermo True --dft True` を追加すると、**MEP探索 → TS精密化 → IRC → 熱化学解析 → DFT一点計算** を一度に実行できます:
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' --tsopt True --thermo True --dft True
```
---

**(i) 2つ以上の完全なタンパク質-リガンドPDB** `.pdb`（R → … → P）、**または (ii) `--scan-lists` 付きの1つのPDB**、**または (iii) `--tsopt True` 付きの1つのTS候補** が与えられると、`pdb2reaction` は自動的に:

- ユーザー定義の基質周辺の**活性部位**を抽出して**クラスターモデル**を構築
- Growing String Method (GSM) や Direct Max Flux (DMF) などの経路最適化手法で**最小エネルギー経路 (MEP)** を探索
- _オプションで_**遷移状態**を精密化し、**振動解析**、**IRC計算**、**DFT一点計算**を実行

MetaのUMAモデルを使用した**機械学習原子間ポテンシャル (MLIP)** で計算を行います。

これらすべてはコマンドラインインターフェース (CLI) を通じて提供され、**多段階酵素反応メカニズム**を最小限の手動介入で生成できるように設計されています。もちろん、このツールキットは小分子系も扱えます。全構造でワークフローを実行する場合（つまり `--center/-c` と `--ligand-charge` を省略）、`.xyz` または `.gjf` 形式の入力構造も使用できます。

**HPCクラスターやマルチGPUワークステーション**では、`pdb2reaction` はノード間でUMA推論自体を並列化することで、大規模なクラスターモデル（およびオプションで**完全なタンパク質-リガンド複合体**）を処理できます。並列計算を有効にするには `workers` と `workers_per_node` を設定してください。設定の詳細については [UMA計算機](uma_pysis.md) を参照してください。

```{important}
- 入力PDBファイルには**水素原子**が含まれている必要があります。
- 複数のPDBを提供する場合、**同じ原子が同じ順序**で含まれている必要があります（座標のみ異なる可能性があります）。そうでない場合はエラーが発生します。
- ブール値CLIオプションは `True`/`False` として明示的に渡されます（例: `--tsopt True`）。
```

```{tip}
初めて使う場合は、まず [概念とワークフロー](concepts.md) を読むと全体像が掴みやすいです。  
セットアップや実行でエラーに遭遇したら [トラブルシューティング](troubleshooting.md) も参照してください。
```


### 水素原子付与の推奨ツール

PDBに水素原子がない場合、pdb2reactionを実行する前に以下のツールを使用してください：

| ツール | コマンド例 | 備考 |
|--------|------------|------|
| **reduce** (Richardson Lab) | `reduce input.pdb > output.pdb` | 高速、結晶構造に広く使用 |
| **pdb2pqr** | `pdb2pqr --ff=AMBER input.pdb output.pqr` | 水素を追加し部分電荷を割り当て |
| **Open Babel** | `obabel input.pdb -O output.pdb -h` | 汎用ケモインフォマティクスツールキット |

複数のPDB入力で同一の原子順序を確保するには、すべての構造に同じ水素付与ツールを一貫した設定で適用してください。

```{warning}
このソフトウェアはまだ開発中です。自己責任でご使用ください。
```

---

## インストール

`pdb2reaction` は、CUDA対応GPUを備えたLinux環境（ローカルワークステーションまたはHPCクラスター）向けに設計されています。特に **PyTorch**、**fairchem-core (UMA)**、**gpu4pyscf-cuda12x** などの依存関係は、動作するCUDAインストールを前提としています。

詳細は上流プロジェクトを参照してください:

- fairchem / UMA: <https://github.com/facebookresearch/fairchem>, <https://huggingface.co/facebook/UMA>
- Hugging Faceトークンとセキュリティ: <https://huggingface.co/docs/hub/security-tokens>

### クイックスタート

以下は多くのCUDA 12.9クラスターで動作する最小限のセットアップ例です。モジュール名とバージョンはシステムに合わせて調整してください。以下はデフォルトのGSM MEPモード（DMFなし）を想定しています。DMFを使用する場合は、最初にcondaでcyipoptをインストールしてください。

```bash
# 1) CUDA対応のPyTorchビルドをインストール
# 2) GitHubからpdb2reactionをインストール
# 3) Plotly図表エクスポート用のヘッドレスChromeをインストール

pip install torch --index-url https://download.pytorch.org/whl/cu129
pip install git+https://github.com/t-0hmura/pdb2reaction.git
plotly_get_chrome -y
```

最後に、UMAモデルをダウンロードできるように **Hugging Face Hub** にログインします:

```bash
# Hugging Face CLI
hf auth login --token '<YOUR_ACCESS_TOKEN>' --add-to-git-credential
```

または

```bash
# クラシックCLI
huggingface-cli login
```

これはマシン/環境ごとに1回だけ行う必要があります。

- MEP探索でDirect Max Flux法を使用する場合は、conda環境を作成し、インストール前にcyipoptをインストールしてください:
  ```bash
  # 専用のconda環境を作成してアクティブ化
  conda create -n pdb2reaction python=3.11 -y
  conda activate pdb2reaction

  # cyipoptをインストール（MEP探索のDMF法に必要）
  conda install -c conda-forge cyipopt -y
  ```

- *環境モジュール*を使用するHPCクラスターでは、PyTorchをインストールする**前に**CUDAをロードしてください:
  ```bash
  module load cuda/12.9
  ```


### ステップバイステップインストール

環境を段階的に構築する場合:

1. **CUDAをロード（HPCで環境モジュールを使用する場合）**

   ```bash
   module load cuda/12.9
   ```

2. **conda環境を作成してアクティブ化**

   ```bash
   conda create -n pdb2reaction python=3.11 -y
   conda activate pdb2reaction
   ```

3. **cyipoptをインストール**
   MEP探索でDMF法を使用する場合に必要です。

   ```bash
   conda install -c conda-forge cyipopt -y
   ```

4. **適切なCUDAビルドのPyTorchをインストール**

   CUDA 12.9の場合:

   ```bash
   pip install torch --index-url https://download.pytorch.org/whl/cu129
   ```

   （クラスターが推奨する場合は別の互換バージョンを使用できます。）

5. **`pdb2reaction` 本体と可視化用Chromeをインストール**

   ```bash
   pip install git+https://github.com/t-0hmura/pdb2reaction.git
   plotly_get_chrome -y
   ```

6. **Hugging Face Hub (UMAモデル) にログイン**

   ```bash
   huggingface-cli login
   ```

   参照:

   - <https://github.com/facebookresearch/fairchem>
   - <https://huggingface.co/facebook/UMA>
   - <https://huggingface.co/docs/hub/security-tokens>

---

## コマンドラインの基本

メインのエントリーポイントは `pip` でインストールされる `pdb2reaction` コマンドです。内部的には **Click** ライブラリを使用しており、デフォルトのサブコマンドは `all` です。

つまり:

```bash
pdb2reaction [OPTIONS] ...
# は以下と同等
pdb2reaction all [OPTIONS] ...
```

`all` ワークフローは**オーケストレーター**です: クラスター抽出、MEP探索、TS最適化、振動解析、オプションのDFT一点計算を単一コマンドに連鎖させます。

クラスター抽出を行う場合、すべての高レベルワークフローは2つの重要なオプションを共有します:

- `-i/--input`: 1つ以上の**完全構造**（反応物、中間体、生成物）
- `-c/--center`: **基質/抽出中心**の定義方法（例: 残基名または残基ID）

`--center/-c` を省略すると、クラスター抽出はスキップされ、**完全な入力構造**が直接使用されます。

---

## メインワークフローモード

### 複数構造MEPパイプライン（反応物 → 生成物）

推定反応座標に沿った複数の完全PDB構造（例: R → I1 → I2 → P）がすでにある場合に使用します。

**最小例**

```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
```

**詳細例**

```bash
pdb2reaction -i R.pdb I1.pdb I2.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' --out-dir ./result_all --tsopt True --thermo True --dft True
```

動作:

- 反応順序で2つ以上の**完全系**を受け取る
- 各構造の触媒クラスターモデルを抽出
- デフォルトで `path_search` による**再帰的MEP探索**を実行
- `--refine-path False` でオプションで**シングルパス** `path-opt` 実行に切り替え
- PDBテンプレートが利用可能な場合、クラスターモデルMEPを**完全系**にマージ
- オプションで各セグメントに対してTS最適化、振動解析、DFT一点計算を実行

これは適度に間隔を置いた中間体（ドッキング、MD、または手動モデリングなど）を生成できる場合に推奨されるモードです。

```{important}
`pdb2reaction` は複数の入力PDBが**まったく同じ原子を同じ順序**で含むことを前提としています（座標のみ異なる可能性があります）。入力間で座標以外のフィールドが異なる場合、エラーが発生します。入力PDBファイルには**水素原子**も含まれている必要があります。
```

---

### 単一構造 + 段階的スキャン（MEP精密化への入力）

**1つのPDB構造**しかないが、反応に沿ってどの原子間距離が変化するかを知っている場合に使用します。

`--scan-lists` と一緒に単一の `-i` を指定します:

**最小例**

```bash
pdb2reaction -i R.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' --scan-lists '[("TYR 285 CA","MMT 309 C10",2.20),("TYR 285 CB","MMT 309 C11",1.80)]' '[("TYR 285 CB","MMT 309 C11",1.20)]'
```

**詳細例**

```bash
pdb2reaction -i SINGLE.pdb -c 'SAM,GPP' --scan-lists '[("TYR 285 CA","MMT 309 C10",2.20),("TYR 285 CB","MMT 309 C11",1.80)]' '[("TYR 285 CB","MMT 309 C11",1.20)]' --mult 1 --out-dir ./result_scan_all --tsopt True --thermo True --dft True
```

ポイント:

- `--scan-lists` は抽出されたクラスターモデルでの**段階的距離スキャン**を記述
- 各タプル `(i, j, target_Å)` は:
  - `'TYR,285,CA'` のようなPDB原子セレクター文字列（**区切り文字**: スペース/カンマ/スラッシュ/バッククォート/バックスラッシュ ` ` `,` `/` `` ` `` `\`）**または**1始まりの原子インデックス
  - クラスターモデルのインデックスに自動的に再マッピング
- 1つの `--scan-lists` リテラルは単一スキャンステージを実行; 複数のリテラルは順次ステージを実行。単一フラグの後に複数のリテラルを渡します（フラグの繰り返しは不可）
- 各ステージは `stage_XX/result.pdb` を書き出し、候補中間体または生成物として扱われる
- デフォルトの `all` ワークフローは連結されたステージを再帰的 `path_search` で精密化
- `--refine-path False` を使用すると、シングルパス `path-opt` チェーンを実行し、再帰的精密化をスキップ（マージされた `mep_w_ref*.pdb` なし）

このモードは単一構造から反応経路を構築するのに便利です。

---

### 単一構造TSOPTのみモード

すでに**遷移状態候補**があり、それを精密化してIRC計算を行いたい場合に使用します。

正確に1つのPDBを指定し、`--tsopt` を有効にします:

**最小例**

```bash
pdb2reaction -i TS_CANDIDATE.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' --tsopt True
```

**詳細例**

```bash
pdb2reaction -i TS_CANDIDATE.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' --tsopt True --thermo True --dft True --out-dir ./result_tsopt_only
```

動作:

- MEP/経路探索を完全にスキップ
- TS最適化で**クラスターモデルTS**を精密化
- 両方向で**IRC**を実行し、両端点を最適化してRおよびP極小に緩和
- R/TS/Pに対して `freq` と `dft` を実行可能
- UMA、Gibbs、DFT//UMAエネルギーダイアグラムを生成

`energy_diagram_*_all.png` や `irc_plot_all.png` などの出力は、トップレベルの `--out-dir` の下にミラーされます。

```{important}
単一入力実行には **`--scan-lists`**（段階的スキャン → GSM）**または** `--tsopt True`（TSOPTのみ）のいずれかが必要です。これらのいずれかなしで単一の `-i` のみを指定しても、完全なワークフローはトリガーされません。
```

---

## 重要なCLIオプションと動作

以下はワークフロー全体で最もよく使用されるオプションです。

| オプション | 説明 |
|----------|------|
| `-i, --input PATH...` | 入力構造。**2つ以上のPDB** → MEP探索; **1つのPDB + `--scan-lists`** → 段階的スキャン → GSM; **1つのPDB + `--tsopt True`** → TSOPTのみモード |
| `-c, --center TEXT` | 基質/抽出中心を定義。残基名（`'SAM,GPP'`）、残基ID（`A:123,B:456`）、またはPDBパスをサポート |
| `--ligand-charge TEXT` | 電荷情報: マッピング（`'SAM:1,GPP:-3'`）または単一整数 |
| `-q, --charge INT` | 総電荷のハードオーバーライド |
| `-m, --mult INT` | スピン多重度（例: シングレットは `1`）。注: `all` 以外のサブコマンドでは `--multiplicity` を使用してください。 |
| `--scan-lists TEXT...` | 単一入力実行用の段階的距離スキャン |
| `--out-dir PATH` | トップレベル出力ディレクトリ |
| `--tsopt {True\|False}` | TS最適化とIRCを有効化 |
| `--thermo {True\|False}` | 振動解析と熱化学を実行 |
| `--dft {True\|False}` | DFT一点計算を実行 |
| `--refine-path {True\|False}` | 再帰的MEP精密化（デフォルト）vs シングルパス |
| `--opt-mode light\|heavy` | 最適化手法: Light (LBFGS/Dimer) または Heavy (RFO/RS-I-RFO) |
| `--mep-mode gsm\|dmf` | MEP手法: Growing String Method または Direct Max Flux |
| `--hessian-calc-mode Analytical\|FiniteDifference` | ヘシアン計算モード。**VRAMが十分な場合はAnalytical推奨** |

すべてのオプションとYAMLスキーマについては [all](all.md) および [YAMLリファレンス](yaml-reference.md) を参照してください。

---

## 実行サマリー

すべての `pdb2reaction all` 実行は以下を書き出します:

- `summary.log` – クイック検査用のフォーマット済みサマリー
- `summary.yaml` – YAMLバージョンサマリー

通常含まれる内容:

- 呼び出された正確なCLIコマンド
- グローバルMEP統計（例: 最大障壁、経路長）
- セグメントごとの障壁高さと主要な結合変化
- UMA、熱化学、DFT後処理からのエネルギー（有効な場合）

各 `path_search` セグメントディレクトリにも独自の `summary.log` と `summary.yaml` があり、ローカルな精密化を個別に検査できます。

---

## CLIサブコマンド

ほとんどのユーザーは主に `pdb2reaction all` を呼び出しますが、CLIは `pdb2reaction opt` のようなサブコマンドも公開しています。各サブコマンドは `-h/--help` をサポートしています。

| サブコマンド | 役割 | ドキュメント |
|------------|------|------------|
| `all` | エンドツーエンドワークフロー | [all](all.md) |
| `extract` | 活性部位クラスター抽出 | [extract](extract.md) |
| `opt` | 構造最適化 | [opt](opt.md) |
| `tsopt` | 遷移状態最適化 | [tsopt](tsopt.md) |
| `path-opt` | MEP最適化 (GSM/DMF) | [path_opt](path_opt.md) |
| `path-search` | 再帰的MEP探索 | [path_search](path_search.md) |
| `scan` | 1D結合長スキャン | [scan](scan.md) |
| `scan2d` | 2D距離スキャン | [scan2d](scan2d.md) |
| `scan3d` | 3D距離スキャン | [scan3d](scan3d.md) |
| `irc` | IRC計算 | [irc](irc.md) |
| `freq` | 振動解析 | [freq](freq.md) |
| `dft` | DFT一点計算 | [dft](dft.md) |
| `trj2fig` | エネルギープロファイルプロット | [trj2fig](trj2fig.md) |
| `add-elem-info` | PDB元素カラム修復 | [add_elem_info](add_elem_info.md) |

```{important}
コマンド（**`all` 以外**すべて）は、`extract` で生成された**クラスターモデル**を入力として想定しています。これらのクラスターモデルでは、Link-Hキャップに最も近い原子が自動的に**凍結**されます。クラスターモデルを自分で構築する場合は、Link-H残基名を `LKH`、原子名を `HL` に設定するか、`--args-yaml` → `geom.freeze_atoms` で凍結する原子を指定してください。
```

```{tip}
`all`、`tsopt`、`freq`、`irc` では、VRAMが十分にある場合は **`--hessian-calc-mode Analytical`** を設定することを強く推奨します。
```

---

## ヘルプ

任意のサブコマンドについて:

```bash
pdb2reaction <subcommand> --help
```

これは利用可能なオプション、デフォルト、および短い説明を表示します。詳細なUMA計算機オプションについては [UMA計算機](uma_pysis.md) を参照してください。

問題が発生した場合は、[GitHubリポジトリ](https://github.com/t-0hmura/pdb2reaction) でIssueを開いてください。
