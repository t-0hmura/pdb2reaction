# インストール

`pdb2reaction` は Linux 環境（ローカルワークステーションまたは HPC クラスター）向けで、本番計算では通常 CUDA 対応 GPU を使用します。prebuilt の **PyTorch** wheel は CUDA runtime library を同梱するため、必要なのは互換 NVIDIA driver であり、local CUDA toolkit ではありません。toolkit が必要なのは CUDA extension や GPU package を source build する場合です。

詳細は上流プロジェクトを参照してください:

- fairchem / UMA: <https://github.com/facebookresearch/fairchem>, <https://huggingface.co/facebook/UMA>
- Hugging Face トークンとセキュリティ: <https://huggingface.co/docs/hub/security-tokens>

## クイックスタート

以下はデフォルトの GSM による MEP 探索（`--mep-mode gsm`）を前提とした最小セットアップです。DMF（`--mep-mode dmf`）を使用する場合は、先に conda で cyipopt をインストールしてください。PyTorch 2.8.0 は `cu126`、`cu128`、`cu129` wheel を配布しているため、site で検証済みの driver・GPU architecture に合う index を選びます。

### 必須

```bash
# 1) CUDA 対応の PyTorchビルドをインストール
# 2) pdb2reactionをインストール
# 3) Plotly 静的画像 (PNG) エクスポート用のヘッドレス Chrome をインストール
#    Chromium binaryをdownload（internet接続が必要）

TORCH_INDEX=cu126  # GPU/site stack が必要とする場合は cu128/cu129
pip install 'torch==2.8.0' --index-url "https://download.pytorch.org/whl/${TORCH_INDEX}"
pip install pdb2reaction
plotly_get_chrome -y
```

最後に、UMA モデルをダウンロードできるように **Hugging Face Hub** にログインします（無料の HF アカウントと読み取り専用トークンが必要。<https://huggingface.co/facebook/UMA> でモデルライセンスの承認が必要な場合あり）:

```bash
hf auth login
# またはスクリプト内でトークン指定する場合:
hf auth login --token '<YOUR_ACCESS_TOKEN>' --add-to-git-credential
```

（新しい `huggingface_hub` は `hf` CLI を同梱しています。旧バージョンでは `huggingface-cli login` が引き続き利用できますが、これは非推奨化が進められています。）

これはマシン/環境ごとに 1 回だけ行う必要があります。

### 任意

- MEP 探索で Direct Max Flux（DMF）法を使用する場合は、conda 環境を作成し、pdb2reaction のインストール前に cyipopt をインストールしてください。

  ```bash
  # 専用のconda環境を作成してアクティブ化
  conda create -n <your-env> python=3.11 -y
  conda activate <your-env>

  # cyipoptをインストール（MEP 探索のDMF法に必要）
  conda install -c conda-forge cyipopt -y
  ```

- HPC site が local build extension 用に CUDA module を要求する場合のみ、site 指定の module をロードします。prebuilt PyTorch wheel のインストールだけのために別の CUDA runtime を追加せず、まず NVIDIA driver と wheel の組合せを検証してください。

  ```bash
  module load cuda/<your-version>   # 例: cuda/12.6 または cuda/12.9
  ```

> **ヒント:** UMA がデフォルトの MLIP バックエンドです。ORB や AIMNet2 を使用するには、対応する extra をインストール（例: `pip install "pdb2reaction[orb]"`）し、コマンドに `-b/--backend orb` を渡してください。下の手順 7 を参照。

```{warning}
**MACE:** `mace-torch` は `e3nn==0.4.4` を要求し、`fairchem-core` の `e3nn>=0.5` pin（UMA）と競合します。両者は共存できないため、MACE には専用の conda env が必要です。標準 recipe はその env で `pip uninstall -y fairchem-core && pip install mace-torch` です。
```


(ja-step-by-step-installation)=
## 詳細なインストール手順

環境を段階的に構築する場合:

1. **site/build が必要とする場合のみ CUDA toolkit をロード**

    prebuilt PyTorch wheel に `nvcc` は不要です。依存 package を source
    build する場合は `module avail cuda` を確認し、cluster が指定する
    compiler/toolkit の組合せをロードしてください:

    ```bash
    module load cuda/<your-version>
    ```

2. **conda 環境を作成してアクティブ化**

    ```bash
    conda create -n <your-env> python=3.11 -y
    conda activate <your-env>
    ```

3. **cyipopt をインストール**
    MEP 探索で DMF 法（`--mep-mode dmf`）を使用する場合に必要です。GSM のみを使用する場合はスキップできます。

    ```bash
    conda install -c conda-forge cyipopt -y
    ```

4. **適切な CUDA ビルドの PyTorch をインストール**

    Blackwell より前の GPU に対する保守的な例:

    ```bash
    pip install 'torch==2.8.0' --index-url https://download.pytorch.org/whl/cu126
    ```

    公式 2.8.0 matrix には `cu128`、`cu129`、`cpu` もあります。driver と
    GPU architecture で選び、`torch.cuda.is_available()` で検証します。
    `nvidia-smi` の "CUDA Version" 表示を local toolkit の選択値として扱わないでください。

5. **`pdb2reaction` 本体と可視化用 Chrome をインストール**

    ```bash
    pip install pdb2reaction
    plotly_get_chrome -y
    ```

6. **Hugging Face Hub (UMA モデル) にログイン**

    ```bash
    hf auth login
    ```

    参照:

    - <https://github.com/facebookresearch/fairchem>
    - <https://huggingface.co/facebook/UMA>
    - <https://huggingface.co/docs/hub/security-tokens>

7. **（任意）追加の MLIP バックエンドをインストール**

    pdb2reaction はデフォルトで UMA を使用します。他のバックエンドを使用する場合は、対応するオプション依存関係をインストールしてください:

    ```bash
    # ORB バックエンド
    pip install "pdb2reaction[orb]"

    # AIMNet2 バックエンド
    pip install "pdb2reaction[aimnet]"

    # MACE バックエンド（mace-torch が要求する e3nn==0.4.4 が UMA の
    # fairchem-core と衝突するため、別 conda 環境で実施してください）
    conda create -n <mace-env> python=3.11 -y && conda activate <mace-env> \
        && pip install pdb2reaction \
        && pip uninstall -y fairchem-core \
        && pip install mace-torch

    # DFT 一点計算の後処理（`--dft` / `pdb2reaction dft`）
    # gpu4pyscf-cuda12x、PySCF、および関連依存をインストールします。
    # 注: gpu4pyscf-cuda12x は PyPI で x86_64 wheel を配布。aarch64 では
    # ソースからビルドしてください (https://github.com/pyscf/gpu4pyscf)。
    pip install "pdb2reaction[dft]"
    ```

8. **インストールの確認**

    ```bash
    pdb2reaction --version
    ```

    インストールされたバージョンが表示されます。GPU アクセスの確認:

    ```bash
    python -c "import torch; print('CUDA:', torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'N/A')"
    ```

    `CUDA: False` の場合、version を変える前に installed wheel、scheduler
    の GPU visibility、driver、environment library を確認してください:

    ```bash
    python -m torch.utils.collect_env
    python -m pip check
    ```

## システム要件

**GPU / CUDA / VRAM:** PyTorch 2.8.0 の公式 CUDA wheel（`cu126`、`cu128`、`cu129`）から、driver と GPU architecture の両方に対応するものを選びます。新しい GPU architecture では新しい wheel が必要なことがありますが、prebuilt wheel と同じ番号の local toolkit は不要です。必要VRAMはbackend/model、原子数、Hessian mode、precision、active自由度に依存します。代表的なproduction stageをpilot実行してpeak allocationを測定してください。smoke suiteはcorrectness checkでありproduction memoryの見積りではありません。

**RAM:** 代表runで必要量を測定します。dense Hessian、model load、並行worker/process stageが支配的になり得ます。

**ディスク:** 選択したenvironment、backend weight cache、生成trajectory/Hessian、`plotly_get_chrome` が任意installするChromiumを含めて見積もります。production前に対象filesystem上の実sizeを確認してください。

CPU のみでも実行できますが、通常は大幅に遅くなります。backend/model ごとに測定し、固定の GPU/CPU 比を仮定しないでください。

## 次のステップ

- [はじめに](getting-started.md) — プロジェクト概要、パイプラインの各ステージ、ワークフローモード
- [クイックスタート: `pdb2reaction all`](quickstart-all.md) — 2 つの PDB から end-to-end 実行
- [クイックスタート: 単一構造スキャン](quickstart-scan.md) — `--scan-lists` で 1 つの PDB から MEP
- [クイックスタート: TS のみモード](quickstart-tsopt-freq.md) — TS 候補を end-to-end で検証
- [CLI 規約](cli-conventions.md) — フラグの優先順位、原子/残基セレクタ、共通オプション
- [トラブルシューティング](troubleshooting.md) と [典型エラー別レシピ](recipes-common-errors.md)
