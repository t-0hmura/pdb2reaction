# インストール

`pdb2reaction` は、CUDA 対応 GPU を備えた Linux 環境（ローカルワークステーションまたは HPC クラスター）での動作を前提としています。特に **PyTorch**、**fairchem-core (UMA)**、**gpu4pyscf-cuda12x** などの依存関係は、動作する CUDA インストールを必要とします。

詳細は上流プロジェクトを参照してください:

- fairchem / UMA: <https://github.com/facebookresearch/fairchem>, <https://huggingface.co/facebook/UMA>
- Hugging Face トークンとセキュリティ: <https://huggingface.co/docs/hub/security-tokens>

## クイックスタート

以下は多くの CUDA 12.9 クラスターで動作する最小限のセットアップ例です。モジュール名やバージョンはお使いの環境に合わせて調整してください。この例はデフォルトの GSM による MEP 探索（`--mep-mode gsm`）を前提としています。DMF（`--mep-mode dmf`）を使用する場合は、先に conda で cyipopt をインストールしてください。

### 必須

```bash
# 1) CUDA 対応の PyTorchビルドをインストール
# 2) pdb2reactionをインストール
# 3) Plotly 静的画像 (PNG) エクスポート用のヘッドレス Chrome をインストール
#    ~150 MB の Chromium バイナリをダウンロード（インターネット接続必要）

pip install torch --index-url https://download.pytorch.org/whl/cu129
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

- *環境モジュール*を使用する HPC クラスターでは、PyTorch をインストールする**前に** CUDA をロードしてください。`module avail cuda` で利用可能なバージョンを確認し、ターゲットの PyTorch wheel に合うバージョン（例: `cu126` ↔ CUDA 12.6、`cu129` ↔ CUDA 12.9）をロードしてください:

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

1. **CUDA をロード（HPC で環境モジュールを使用する場合）**

    `module avail cuda` で利用可能なバージョンを確認し、ターゲットの
    PyTorch wheel に合うバージョン（例: `cu126` は CUDA 12.6、`cu129`
    は CUDA 12.9）をロードしてください:

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

    CUDA 12.9 の場合:

    ```bash
    pip install torch --index-url https://download.pytorch.org/whl/cu129
    ```

    PyTorch は CUDA ドライバーバージョンに合わせたビルドが必要です。互換性は [PyTorch Get Started](https://pytorch.org/get-started/locally/) で確認してください。CPU のみの実行もサポートされますが、大幅に遅くなります（10-100 倍）。

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
    # ORB バックエンド。orb-models が引く torch_scatter は prebuilt wheel が PyPI でなく
    # PyG の index にあるため、pip はソースビルドし、build isolation で
    # "No module named 'torch'" と失敗することがあります。torch+CUDA タグに合わせて PyG index を追加:
    #   pip install "pdb2reaction[orb]" -f https://data.pyg.org/whl/torch-2.8.0+cu129.html
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

    暗黙溶媒補正を使用するには、[xTB](https://github.com/grimme-lab/xtb) をインストールし、`xtb` コマンドが `PATH` 上で利用可能であることを確認してください。

    #### xTB のインストール

    **ALPB 溶媒和モデルの場合**（推奨の出発点）:

    ```bash
    conda install -c conda-forge xtb
    ```

    **CPCM-X 溶媒和モデルの場合** conda-forge の xtb は CPCM-X を含まないため、ソースからのビルドが必要です。手順は {ref}`ja-recipe-cpcmx-build` を参照。

    カスタム xTB バイナリを使用するには、YAML 設定で `xtb_cmd` キーを設定するか、Python で `calc.xtb_cmd` を使用してください。

8. **インストールの確認**

    ```bash
    pdb2reaction --version
    ```

    インストールされたバージョンが表示されます。GPU アクセスの確認:

    ```bash
    python -c "import torch; print('CUDA:', torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'N/A')"
    ```

    `CUDA: False` の場合、CUDA モジュールのロードと PyTorch ビルドの CUDA バージョンを確認してください。

## システム要件

**GPU / CUDA / VRAM:** 実行環境の CUDA に合った CUDA タグの PyTorch wheel をインストールしてください（CUDA 12.6 なら `cu126`、CUDA 12.9 なら `cu129`。RTX 50 系では 12.9 が必須です）。VRAM は 8 GB 以上を推奨します。大きな ML 領域の解析的Hessianでは GPU 容量が大きいほど有利です。`tests/smoke/` のスイートはデフォルトの `uma-s-1p1` モデルでピーク約 0.9 GB のため、本番の TS / IRC / Hessianのワークフローが収まらない小容量 GPU でも実行できます。

**RAM:** 16 GB 以上を推奨します（GPU 計算と並行して大規模な活性部位モデルを扱う場合は余裕があるほど安心です）。

**ディスク:** 空き容量約 20 GB を見込んでください。内訳は conda 環境（約 8 GB）、UMA の Hugging Face モデルキャッシュ（約 1〜4 GB）、静的 PNG エクスポート用に Plotly が取得するヘッドレス Chromium（約 150 MB、`plotly_get_chrome` 経由）です。

CPU のみでの実行は可能ですが 10〜100 倍遅く、TS / IRC / Hessianの一連のワークフローには推奨されません。

## 次のステップ

- [はじめに](getting-started.md) — プロジェクト概要、パイプラインの各ステージ、ワークフローモード
- [クイックスタート: `pdb2reaction all`](quickstart-all.md) — 2 つの PDB から end-to-end 実行
- [クイックスタート: 単一構造スキャン](quickstart-scan.md) — `--scan-lists` で 1 つの PDB から MEP
- [クイックスタート: TS のみモード](quickstart-tsopt-freq.md) — TS 候補を end-to-end で検証
- [CLI 規約](cli-conventions.md) — フラグの優先順位、原子/残基セレクタ、共通オプション
- [トラブルシューティング](troubleshooting.md) と [典型エラー別レシピ](recipes-common-errors.md)

