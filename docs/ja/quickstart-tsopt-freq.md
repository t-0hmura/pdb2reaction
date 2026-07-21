# クイックスタート: `pdb2reaction all --tsopt`（TS-only モード）

## 目的

既に手元にある TS 候補構造に対して、`tsopt → irc → freq →（任意で dft）` の検証・熱力学計算の一連の処理だけを `pdb2reaction all` で一括実行します。上流の `extract` / `path-opt` はスキップされます。

## 事前に必要なもの

- pdb2reaction がインストール済み（[インストール](installation.md)を参照）
- TS 候補構造 1 つ: `.pdb`（残基／電荷情報を持つため推奨）または `.xyz`
- 電荷指定（次のいずれか **ひとつ**を必ず指定）:
  - `-q INT` — 全電荷を整数で直接指定
  - `-l 'RES:Q,...'` — 残基ごとのリガンド電荷（タンパク質–リガンド複合体 PDB の場合）
  - `.gjf` ヘッダから自動取得（Gaussian 入力を渡した場合のみ）
- `-m/--multiplicity` — デフォルトは `1`（一重項）。ラジカル種では明示が必要です
- `.xyz` 入力では `-q` と `-m` が必須です。`-l` 形式で電荷を渡したい場合は `--ref-pdb cluster.pdb` を併用してください
- TS のみモードに入る条件は **3 つすべて成立**: (1) `-i` 入力がちょうど 1 つ、(2) `--scan-lists` が無い、(3) `--tsopt` が指定されている。そうでない場合 CLI は入力ゲートで `BadParameter` を送出します（`Provide at least two structures with -i/--input in reaction order, or use a single structure with --scan-lists, or a single structure with --tsopt.`）。

## 最小コマンド

```bash
pdb2reaction all -i ts_candidate.pdb -l 'SAM:1,GPP:-3' \
    --tsopt --thermo -o ./result_ts_only
```

XYZ 入力では `-q` と `-m` を明示します:

```bash
pdb2reaction all -i ts_candidate.xyz -q -1 -m 1 -b uma \
    --tsopt --thermo -o ./result_ts_only
```

`--tsopt` で検証チェーンを起動し、`--thermo` で freq から ZPE / Gibbs 補正を追加します。両ステージとも同一バックエンド（デフォルト UMA）で実行されます。

### （任意）DFT 一点計算を追加

```bash
pdb2reaction all -i ts_candidate.pdb -l 'SAM:1,GPP:-3' \
    --tsopt --thermo --dft --dft-func-basis 'wb97m-v/def2-tzvpd' \
    -o ./result_ts_only
```

> **VRAM 注意:** `--dft` は GPU4PySCF で一点計算を実行します。必要メモリは
> 構造・基底・汎関数・精度・software stackに依存するため、対象nodeで代表構造を
> pilot実行し、peak memoryを測定してください。OOM時は `--dft` を外し、より小さい
> 基底または縮小clusterで `pdb2reaction dft` を単独実行するか、より大きいnodeへ
> 移してください。`[dft]` extraのinstallも必要です
> （[インストール](installation.md) を参照）。

## 期待される出力

成功時の出力ツリー:

```text
result_ts_only/
├── summary.log                                # 人間可読サマリー
├── summary.json                               # status: success | partial | failed
└── segments/
    └── seg_01/                                # TS のみモードの成果物
        ├── reactant.pdb                       # 正準 R/TS/P は seg_01/ 直下（TS のみモード）
        ├── ts.pdb
        ├── product.pdb
        ├── ts/
        │   ├── final_geometry.{xyz,pdb}
        │   └── vib/imag_*_trj.xyz             # 虚振動モードアニメーション
        ├── irc/
        │   └── {forward,backward,finished}_irc_trj.xyz
        ├── freq/{R,TS,P}/
        │   ├── frequencies_cm-1.txt
        │   └── thermoanalysis.yaml
        └── dft/{R,TS,P}/                      # --dft 時のみ
            └── result.yaml                    # 常に出力（--dft 時）
```

**確認ポイント（実行順）:**

1. `summary.json` の `status` が `"success"` であること。これは要求した結果が揃い、TS虚振動validatorも通過したことを表します。`"partial"` なら利用可能なpathはありますが、要求stageの失敗・欠損またはvalidator不通過があるため `status_reasons` を確認します。`"failed"` は利用可能なpath結果が無い状態です。`segments[0].barrier_kcal` と `segments[0].delta_kcal` も確認します。
2. `post_segments[0].ts_imag.n_imag == 1` — 一次鞍点の必要条件です。振動数の大きさだけで反応性を判定せず、modeを可視化してIRC接続を確認します。系ごとのnoise評価で柔らかい非反応modeを特定した場合だけ、YAMLの opt-in filter `irc.imag_below` を既定の `0.0` より負側へ設定します。IRCが受理するのは `ν <= imag_below` のmodeです。
3. `segments/seg_01/irc/{forward,backward}_irc_trj.xyz` を PyMOL で開き、R 端・P 端まで到達していることを確認。
4. `segments[0].bond_changes` が空でなく、想定どおりの結合切断・形成が記録されていること。
5. `segments/seg_01/freq/{R,TS,P}/frequencies_cm-1.txt` で R / P に虚振動が無く、TS で 1 つだけあること。

**トラブルシュート:**

| 症状 | 原因 | 対処 |
|---|---|---|
| `post_segments[0].ts_imag.n_imag == 0` | TS 候補が極小に落ちてしまう | `path-search` で TS 候補を取り直します。`all` は MEP 接線を内部で渡して有界回復を行いますが、経路情報のない通常の TS のみモードでは目的の隣接鞍点を特定できません |
| `n_imag >= 2` | 縮退した負モードあり | 既定 `baker` で現れる**近ゼロ**の余剰モード（数 cm⁻¹）は多くが収束アーティファクトです。`baker` は大量計算にコスパが良い既定ですが、`n_imag >= 2` が出たら freq/tsopt を厳しい `--thresh-post`（`gau_tight` 以上）で再実行してください — 通常は `n_imag = 1` に解消します。それでも残る場合は `--flatten` で余剰モードを除去（[tsopt](tsopt.md) の `hessian_dimer.flatten_max_iter` を参照）。 |
| `segments[0].bond_changes` が空（`""` または `(no covalent changes detected)`）、または IRC が想定と違う終点に到達 | 虚振動が反応座標方向と一致していない、または TS が同じ井戸同士を結んでいる（反応物側と生成物側が同一極小） | `segments/seg_01/ts/vib/imag_*_trj.xyz` を PyMOL で可視化し、虚振動が想定の反応方向か確認。違う場合は TS 候補を取り直す |
| `freq/{R,P}/frequencies_cm-1.txt` に虚振動が残る | IRC の終点が真の極小に達していない | 収束をより厳しくする（`--thresh-post gau_tight` など; `baker` はデフォルト値なので指定しても変化なし）か YAML で IRC max cycles を伸ばす。[freq](freq.md) を参照 |

## 補足

- `tsopt` 単独でのパラメータ調整（`--opt-mode`、`--max-cycles`、Hessian オプション）は [tsopt](tsopt.md) を参照。
- デフォルトの `FiniteDifference` を基準にし、選択した backend/model と対象系で速度・メモリ・結果を検証できた場合にだけ `--hessian-calc-mode Analytical` を明示してください。
- 全オプションを確認するには `pdb2reaction all --help-advanced`。

## 次のステップ

- 複数構造からの MEP 経路: [クイックスタート: `pdb2reaction all`](quickstart-all.md)
- 単一構造からのスキャン駆動: [クイックスタート: `pdb2reaction all --scan-lists`](quickstart-scan.md)
- 全オプションリファレンス: [all](all.md) / [tsopt](tsopt.md) / [irc](irc.md) / [freq](freq.md) / [dft](dft.md)
