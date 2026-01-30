# `path-search` サブコマンド

## 概要
反応座標に沿って順序付けられた**2つ以上**の構造にわたって連続的な最小エネルギー経路（MEP）を構築します。`path-search` はGSM**またはDMF**セグメントを連鎖させ、共有結合変化のある領域のみを選択的に精密化し、（オプションで）PDBポケットをフルサイズテンプレートにマージします。**GSMがデフォルト**です。

## 使用法
```bash
pdb2reaction path-search -i R.pdb [I.pdb ...] P.pdb [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [--multiplicity 2S+1]
                         [--workers N] [--workers-per-node N]
                         [--mep-mode {gsm|dmf}] [--freeze-links {True\|False}] [--thresh PRESET]
                         [--refine-mode {peak|minima}]
                         [--max-nodes N] [--max-cycles N] [--climb {True\|False}]
                         [--opt-mode light|heavy] [--dump {True\|False}]
                         [--out-dir DIR] [--preopt {True\|False}]
                         [--align {True\|False}] [--ref-full-pdb FILE ...] [--ref-pdb FILE ...]
                         [--convert-files {True\|False}]
                         [--args-yaml FILE]
```

### 例
- **ポケットのみ**の2つのエンドポイント間のMEP:
  ```bash
  pdb2reaction path-search -i reactant.pdb product.pdb -q 0
  ```
- YAMLオーバーライドとマージされた全系出力を使用した**マルチステップ**探索:
  ```bash
  pdb2reaction path-search \
      -i R.pdb IM1.pdb IM2.pdb P.pdb -q -1 \
      --args-yaml params.yaml --ref-full-pdb holo_template.pdb --out-dir ./run_ps
  ```

## CLIオプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH...` | 反応順序の2つ以上の構造（反応物 → 生成物） | 必須 |
| `-q, --charge INT` | 総電荷 | テンプレート/導出が適用されない限り必須 |
| `--ligand-charge TEXT` | `-q` が省略された場合に使用される総電荷または残基名ごとのマッピング | _None_ |
| `-m, --multiplicity INT` | スピン多重度（2S+1） | `.gjf` テンプレート値または `1` |
| `--freeze-links {True\|False}` | PDBポケットロード時、リンク水素の親原子を凍結 | `True` |
| `--max-nodes INT` | MEPセグメントごとの内部ノード | `10` |
| `--max-cycles INT` | 最大MEP最適化サイクル | `300` |
| `--climb {True\|False}` | GSMセグメントのクライミングイメージを有効化 | `True` |
| `--opt-mode TEXT` | HEI±1/kinkノード用の単一構造オプティマイザー | `light` |
| `--mep-mode {gsm\|dmf}` | セグメント生成器: GSMまたはDMF | `gsm` |
| `--refine-mode {peak\|minima}` | 精密化のシード | _Auto_ |
| `--dump {True\|False}` | MEP（GSM/DMF）と単一構造軌跡/リスタートをダンプ | `False` |
| `--convert-files {True\|False}` | PDBまたはGaussian入力用のXYZ/TRJ → PDB/GJFコンパニオンをトグル | `True` |
| `--out-dir TEXT` | 出力ディレクトリ | `./result_path_search/` |
| `--thresh TEXT` | GSMおよびイメージごとの最適化の収束プリセットをオーバーライド | `gau` |
| `--args-yaml FILE` | YAMLオーバーライド | _None_ |
| `--preopt {True\|False}` | MEP探索前に各エンドポイントを事前最適化（推奨） | `True` |
| `--align {True\|False}` | 探索前にすべての入力を最初の構造にアライメント | `True` |
| `--ref-full-pdb PATH...` | フルサイズテンプレートPDB | _None_ |
| `--ref-pdb PATH...` | 入力がXYZ/GJFの場合のポケット参照PDB | _None_ |

## ワークフロー
1. **ペアごとの初期セグメント（GSM/DMF）** – 各隣接入力（A→B）間で粗いMEPを取得しHEIを特定
2. **HEI周辺の局所緩和** – 選択した単一構造オプティマイザーでHEI±1または最寄りの局所極小を精密化
3. **kink vs. 精密化の決定** – `End1` と `End2` 間に共有結合変化がなければkinkとして扱い、そうでなければ精密化セグメントを起動
4. **選択的再帰** – サブインターバルに共有結合更新が含まれる場合のみ再帰
5. **スティッチング & ブリッジング** – 解決されたサブパスを連結
6. **アライメント & マージング（オプション）** – `--ref-full-pdb` でポケット軌跡をフルサイズPDBテンプレートにマージ

## 出力
```
out_dir/ (デフォルト: ./result_path_search/)
├─ mep.trj                  # プライマリMEP軌跡
├─ mep.pdb                  # 入力がPDBテンプレートで変換が有効な場合のPDBコンパニオン
├─ mep_w_ref.pdb            # マージされた全系MEP（参照PDB/テンプレートが必要）
├─ mep_w_ref_seg_XX.pdb     # 共有結合変化がある場合のマージされたセグメントごとのパス
├─ summary.yaml             # すべての再帰セグメントの障壁と分類サマリー
├─ mep_plot.png             # ΔEプロファイル（kcal/mol、反応物基準）
├─ energy_diagram_MEP.png   # MEP状態エネルギーダイアグラムの静的エクスポート
└─ seg_000_*/               # セグメントごとのGSM/DMFダンプ、HEIスナップショット
```
