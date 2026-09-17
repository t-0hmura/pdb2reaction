# インストール

標準手順は Linux、Python 3.12、NVIDIA GPU 向けです。以下は PyTorch の CUDA 13 wheel を使うため、対応する NVIDIA ドライバーが必要です。別のドライバー・GPU の組合せでは、[PyTorch の対応表](https://pytorch.org/get-started/previous-versions/)から wheel を選んでください。CUDA ランタイムは wheel に含まれます。

(ja-step-by-step-installation)=
## クイックスタート

### 必須

[UMA のモデルライセンス](https://huggingface.co/facebook/UMA)に同意した後、以下をターミナルへコピーしてください。`hf auth login` で Hugging Face の認証情報を入力します。

```bash
conda create -n pdb2reaction python=3.12 pip -y
conda activate pdb2reaction
pip install torch==2.13.0 --index-url https://download.pytorch.org/whl/cu130
pip install pdb2reaction
plotly_get_chrome -y
hf auth login
pdb2reaction --version
```

デフォルトの UMA バックエンドと、PNG 出力用の Chrome が入ります。新しいターミナルでは `conda activate pdb2reaction` で環境を有効にしてください。

### 任意

使う機能だけ、同じ環境を有効にした状態で追加します。

| 機能 | インストール |
|---|---|
| ORB（`-b orb`） | `pip install --only-binary=dm-tree "pdb2reaction[orb]"` |
| AIMNet2（`-b aimnet2`） | `pip install "pdb2reaction[aimnet]"` |
| DFT（`--dft` / `pdb2reaction dft`） | `pip install "pdb2reaction[dft]"` |
| MCP サーバー | `pip install "pdb2reaction[mcp]"` |
| DMF 経路探索（`--mep-mode dmf`） | `conda install -c conda-forge cyipopt "numpy>=2,<2.5" -y` |

ORB は Python 3.11／3.12 が必要です。上記手順では 3.12 を使用します。ORB と AIMNet2 は UMA の認証を必要としません。PyDMF は本体の依存として導入され、DMF の利用には cyipopt が追加で必要です。

DFT extra は x86_64 では GPU4PySCF を導入します。aarch64 では `dft` コマンドに `--engine cpu` を指定して CPU PySCF を使用してください。[DFT](dft.md)も参照してください。

MACE は UMA と e3nn の依存条件が競合するため、専用環境を使用します。[MACE の導入手順](https://github.com/t-0hmura/pdb2reaction/blob/main/skills/pdb2reaction-install-backends/mace.md)を参照してください。

## GPU の確認

```bash
python -c "import torch; print('CUDA:', torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'N/A')"
python -m pip check
```

クラスターでは GPU が割り当てられたジョブ内で確認してください。CUDA が利用できない場合は、ドライバー、PyTorch wheel、ジョブの GPU 割当を確認します。詳細は `python -m torch.utils.collect_env` で取得できます。

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
