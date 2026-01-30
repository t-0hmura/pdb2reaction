# `tsopt` サブコマンド

## 概要
`pdb2reaction tsopt` は2つの補完的なワークフローを使用して遷移状態を最適化します:

- **light** モード: 定期的な正確ヘシアン更新を伴うHessian Dimer探索、余分な虚数モードを除去するためのオプションのメモリ効率的なフラットンループ（デフォルトで無効）、活性自由度のPHVA対応ヘシアン更新
- **heavy** モード: 設定可能な信頼領域セーフガードを持つRS-I-RFOオプティマイザー、収束後に余分な虚数モードが残る場合のオプションの後最適化フラットンループ

両モードはエネルギー/勾配/ヘシアンにUMA計算機を使用し、YAMLから `geom`/`calc`/`opt` 設定を継承し、最終的な虚数モードを常に `.trj` に書き込みます。デフォルトの `--opt-mode` は **heavy**（RS-I-RFO）です; Hessian Dimerワークフローを実行するには `--opt-mode light` に切り替えてください。

## 使用法
```bash
pdb2reaction tsopt -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [-m 2S+1] \
                    [--opt-mode light|heavy] [--flatten-imag-mode {True\|False}] \
                    [--freeze-links {True\|False}] [--max-cycles N] [--thresh PRESET] \
                    [--dump {True\|False}] [--out-dir DIR] [--args-yaml FILE] \
                    [--hessian-calc-mode Analytical|FiniteDifference] \
                    [--convert-files {True\|False}] [--ref-pdb FILE]
```

### 例
```bash
# 推奨ベースライン: 電荷/多重度を指定しlightワークフローを選択
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --opt-mode light --out-dir ./result_tsopt/

# YAMLオーバーライド、有限差分ヘシアン、freeze-links処理を伴うlightモード
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --freeze-links True \
    --opt-mode light --max-cycles 10000 --dump False \
    --out-dir ./result_tsopt/ --args-yaml ./args.yaml \
    --hessian-calc-mode FiniteDifference

# YAMLのみで駆動されるheavyモード（RS-I-RFO）
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --opt-mode heavy \
    --args-yaml ./args.yaml --out-dir ./result_tsopt/
```

## ワークフロー
- **電荷/スピン解決**: 入力が `.gjf` の場合、電荷と多重度はテンプレート値を継承。`-q` が省略されているが `--ligand-charge` が提供されている場合、構造は酵素-基質複合体として扱われ、PDB入力（または `--ref-pdb` 付きXYZ/GJF）で総電荷が導出される
- **構造ロード & freeze-links**: 構造は `pysisyphus.helpers.geom_loader` を介して読み込まれる。PDB入力では、`--freeze-links True` がリンク水素を見つけてその親原子を凍結
- **UMAヘシアン**: `--hessian-calc-mode` は解析的評価と有限差分評価を切り替え。VRAMが十分な場合は `--hessian-calc-mode` を `Analytical` に設定することを強く推奨

## CLIオプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル | 必須 |
| `-q, --charge INT` | 総電荷。`.gjf` テンプレートまたは `--ligand-charge` が提供しない限り必須 | テンプレート/導出が適用されない限り必須 |
| `--ligand-charge TEXT` | `-q` が省略された場合に使用される総電荷または残基名ごとのマッピング | _None_ |
| `-m, --multiplicity INT` | スピン多重度（2S+1） | `.gjf` テンプレート値または `1` |
| `--freeze-links {True\|False}` | PDBのみ。リンク水素の親を凍結 | `True` |
| `--max-cycles INT` | `opt.max_cycles` に転送されるマクロサイクル上限 | `10000` |
| `--opt-mode TEXT` | 上記のLight/Heavyエイリアス | `heavy` |
| `--dump {True\|False}` | 軌跡をダンプ | `False` |
| `--out-dir TEXT` | 出力ディレクトリ | `./result_tsopt/` |
| `--thresh TEXT` | 収束プリセットのオーバーライド | `baker` |
| `--flatten-imag-mode {True\|False}` | 余分な虚数モードフラットンループを有効化 | `False` |
| `--hessian-calc-mode CHOICE` | UMAヘシアンモード（`Analytical` または `FiniteDifference`） | `FiniteDifference` |
| `--convert-files {True\|False}` | PDBまたはGaussian入力用のXYZ/TRJ → PDB/GJFコンパニオンをトグル | `True` |
| `--ref-pdb FILE` | 入力がXYZ/GJFの場合に使用する参照PDBトポロジー | _None_ |
| `--args-yaml FILE` | YAMLオーバーライド（`geom`、`calc`、`opt`、`hessian_dimer`、`rsirfo`） | _None_ |

## 出力（& ディレクトリレイアウト）
```
out_dir/ (デフォルト: ./result_tsopt/)
├─ final_geometry.xyz            # 常に書き込み
├─ final_geometry.pdb            # 入力がPDBの場合（変換有効時）
├─ final_geometry.gjf            # 入力がGaussianの場合（変換有効時）
├─ optimization_all.trj          # --dumpがTrueのときのLightモードダンプ
├─ optimization.trj              # --dumpがTrueのときのHeavyモード軌跡
├─ vib/
│  ├─ final_imag_mode_±XXXX.Xcm-1.trj
│  └─ final_imag_mode_±XXXX.Xcm-1.pdb
└─ .dimer_mode.dat               # Lightモード方向シード
```

## 注意事項
- `--opt-mode` エイリアスは上記のワークフローに正確にマップされる; YAMLキーを手動で調整するよりも意図したアルゴリズム用に1つを選択（デフォルト: `heavy`）
- 虚数モード検出は〜5 cm⁻¹がデフォルト（`hessian_dimer.neg_freq_thresh_cm` で設定可能）
- `--hessian-calc-mode` はYAMLマージ後に `calc.hessian_calc_mode` をオーバーライド
