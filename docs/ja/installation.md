# インストール

`pdb2reaction` は Linux 環境（ローカルワークステーションまたは HPC クラスター）向けで、本番計算では通常 CUDA 対応 GPU を使用します。ビルド済みの **PyTorch** wheel は CUDA ランタイムを同梱するため、通常は互換性のある NVIDIA ドライバーだけが必要です。CUDA 拡張や GPU パッケージをソースからビルドする場合は、CUDA toolkit も必要です。

詳細は上流プロジェクトを参照してください:

- fairchem / UMA: <https://github.com/facebookresearch/fairchem>, <https://huggingface.co/facebook/UMA>
- Hugging Face トークンとセキュリティ: <https://huggingface.co/docs/hub/security-tokens>

## クイックスタート

以下はデフォルトの GSM による MEP 探索（`--mep-mode gsm`）を前提とした最小セットアップです。DMF（`--mep-mode dmf`）を使用する場合は、先に conda で cyipopt をインストールしてください。PyTorch 2.13.0 は `cu126`、`cu130`、`cu132` wheel を配布しているため、site で検証済みの driver・GPU architecture に合う index を選びます。

### 必須

```bash
# 1) CUDA 対応の PyTorchビルドをインストール
# 2) pdb2reactionをインストール
# 3) Plotly 静的画像 (PNG) エクスポート用のヘッドレス Chrome をインストール
#    Chromium binaryをdownload（internet接続が必要）

TORCH_INDEX=cu130  # GPU/site stack が必要とする場合は cu126/cu132
pip install 'torch==2.13.0' --index-url "https://download.pytorch.org/whl/${TORCH_INDEX}"
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

DMF を使う場合は、必須コマンドの実行前に環境を作成し、cyipopt を導入します（下記の手順 2〜3）。CUDA module はソースビルドで必要な場合のみ使用します（手順 1）。ORB・AIMNet2・MACE・DFT の導入は手順 7 を参照してください。MACE の `e3nn==0.4.4` は UMA の `e3nn>=0.5` と競合するため、専用環境が必要です。

(ja-step-by-step-installation)=
## 詳細なインストール手順

環境を段階的に構築する場合:

1. **site/build が必要とする場合のみ CUDA toolkit をロード**

    prebuilt PyTorch wheel に `nvcc` は不要です。依存 package を source
    build する場合は `module avail cuda` を確認し、cluster が指定する
    compiler/toolkit の組合せをロードしてください:

    ```bash
    module load cuda/<your-version>   # 例: cuda/12.6 または cuda/12.9
    ```

2. **conda 環境を作成してアクティブ化**

    ```bash
    conda create -n <your-env> python=3.12 -y
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
    pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cu126
    ```

    公式 2.13.0 matrix には `cu130`、`cu132`、`cpu` もあります。driver と
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

    利用許諾と非対話型ログインは、上の「必須」を参照してください。

7. **（任意）追加の MLIP バックエンドをインストール**

    pdb2reaction はデフォルトで UMA を使用します。他のバックエンドは対応する extra を導入し、`-b/--backend`（例: `-b orb`）で選択します:

    ```bash
    # ORB バックエンド（Python 3.11／3.12 が必要、3.12 推奨）
    pip install --only-binary=dm-tree "pdb2reaction[orb]"

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

**GPU / CUDA / VRAM:** PyTorch 2.13.0 の公式 CUDA wheel（`cu126`、`cu130`、`cu132`）から、ドライバーと GPU アーキテクチャの両方に対応するものを選びます。新しい GPU では新しい wheel が必要なことがありますが、同じ番号の CUDA toolkit は不要です。必要な VRAM はバックエンド・モデル、原子数、Hessian 計算モード、精度、Active DOF に依存します。本番計算で使う代表的な処理を試行し、最大メモリ使用量を測定してください。スモークテストは動作確認用で、本番計算のメモリ見積もりには使えません。

**RAM:** 代表的な計算で必要量を測定します。密な Hessian、モデルの読み込み、ワーカーやプロセスの並列実行が主なメモリ消費要因になる場合があります。

**ディスク:** 使用環境、バックエンドのモデル重みキャッシュ、生成する軌跡・Hessian、任意で `plotly_get_chrome` がインストールする Chromium を含めて見積もります。本番計算前に、保存先のファイルシステムで実際の使用量を確認してください。

CPU のみでも実行できますが、通常は大幅に遅くなります。バックエンド・モデルごとに測定し、固定の GPU/CPU 比を仮定しないでください。

## 次のステップ

- [はじめに](getting-started.md) — プロジェクト概要、パイプラインの各ステージ、ワークフローモード
- [クイックスタート: `pdb2reaction all`](quickstart-all.md) — 2 つの PDB から end-to-end 実行
- [クイックスタート: 単一構造スキャン](quickstart-scan.md) — `--scan-lists` で 1 つの PDB から MEP
- [クイックスタート: TS のみモード](quickstart-tsopt-freq.md) — TS 候補を end-to-end で検証
- [CLI 規約](cli-conventions.md) — フラグの優先順位、原子/残基セレクタ、共通オプション
- [トラブルシューティング](troubleshooting.md) と [典型エラー別レシピ](recipes-common-errors.md)
