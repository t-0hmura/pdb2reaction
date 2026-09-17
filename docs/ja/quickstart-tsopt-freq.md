# クイックスタート: `pdb2reaction all --tsopt`（TS-only モード）

## 目的

手元の TS 候補に対して、MEP を省き、`pdb2reaction all --tsopt` で `tsopt → irc` を実行します。`--thermo` で `freq`、`--dft` で DFT 一点計算を追加できます。PDB/mmCIF の抽出は `-c` 指定時のみ行います。

## 事前に必要なもの

- pdb2reaction がインストール済み（[インストール](installation.md)を参照）
- TS 候補構造 1 つ: PDB/mmCIF（残基／電荷情報を持つため推奨）、XYZ、または GJF
- 電荷は `-q`、`-l`、GJF ヘッダー、または設定ファイルで指定します。優先順位は [電荷の指定](cli-conventions.md#電荷の指定)を参照してください。
- 多重度は `-m` → YAML `calc.spin` → GJF ヘッダー → `1` の順で決まります。
- TS-only モードは、入力が 1 つ、`--scan-lists` なし、`--tsopt` 指定で選択されます。2 構造以上なら MEP、単一構造と `--scan-lists` ならスキャンを実行します。

## 最小コマンド

```bash
pdb2reaction all -i ts_candidate.pdb -l 'SAM:1,GPP:-3' \
    --tsopt --thermo -o ./result_ts_only
```

XYZ 一重項では `-q` を明示し、多重度は省略できます:

```bash
pdb2reaction all -i ts_candidate.xyz -q -1 -b uma \
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
├── summary.log                                # 実行要約
├── summary.json                               # status: success | partial | failed
└── segments/
    └── seg_01/                                # TS のみモードの成果物
        ├── reactant.pdb                       # 正準 R/TS/P は seg_01/ 直下（TS のみモード）
        ├── ts.pdb
        ├── product.pdb
        ├── ts/
        │   ├── final_geometry.{xyz,pdb}
        │   └── vib/imag_*_trj.xyz             # 虚振動モード軌跡
        ├── irc/
        │   └── {forward,backward,finished}_irc_trj.xyz
        ├── freq/{R,TS,P}/
        │   ├── frequencies_cm-1.txt
        │   └── thermoanalysis.yaml
        └── dft/{R,TS,P}/                      # --dft 時のみ
            └── result.yaml                    # 常に出力（--dft 時）
```

## 結果の確認

1. **完了状況:** `summary.json` の `scientific_status` と `scientific_status_reasons` を確認します。必要な計算・数値最適化の完了を表し、虚振動数や化学的接続性は別に確認します。
2. **TSのモード:** `post_segments[0].ts_imag.n_imag` が記録された分類基準で `1` か確認します。波数は `nu_imag_max_cm` です。`segments/seg_01/ts/vib/imag_*_trj.xyz` を可視化し、想定する結合の変位を確認してください。大きさだけでは反応性を判定できません。`irc.imag_below`（既定 `0.0` cm⁻¹、`ν <= imag_below` を受理）は、系ごとのノイズ評価後にだけ負側へ変更します。
3. **接続性と構造:** `segments/seg_01/irc/finished_irc_trj.xyz`、同じセグメントの `reactant.pdb`・`ts.pdb`・`product.pdb`、`segments[0].bond_changes` を確認します。原子対応、構造、想定する結合変化を点検してください。TS-onlyモードは高エネルギーのIRC端点をRと呼ぶため、化学的な方向を判断する前に `endpoint_assignment` を確認します。詳細は [all](all.md) を参照してください。
4. **端点の振動:** `segments/seg_01/freq/{R,TS,P}/frequencies_cm-1.txt` を確認します。完全な符号付き振動数を保持し、`n_negative_modes` はすべての負符号を数えます。R/Pに虚振動が残っても熱化学は計算されますが、極小点と判断する前にモードを確認します。
5. **エネルギー:** `rate_limiting_step.barrier_kcal` と `segments[0].delta_kcal` はΔE‡とΔE、`post_segments[0].gibbs_mlip.barrier_kcal` / `.delta_kcal` はΔG‡とΔGです。各状態の `thermoanalysis.yaml` に `electronic_energy_ha`、`zpe_correction_ha`、`sum_EE_and_ZPE_ha`、`sum_EE_and_thermal_free_energy_ha` と温度・圧力（既定298.15 K、1 atm）を記録します。TSまたはPから選択したRの値を引きます。ラベルだけで化学的なR/Pを判断しないでください。

エネルギー差の単位はkcal/mol、`thermoanalysis.yaml` の状態エネルギーはhartreeです。

| 結果 | 次の確認 |
|---|---|
| `n_imag == 0` | TS候補またはMEPを改善します。経路情報のないTS-only計算は目的の隣接鞍点を特定できず、saddle recoveryの既定上限は0です。 |
| `n_imag >= 2` | 各虚振動を確認し、`all --thresh-post gau_tight` または `tsopt --thresh gau_tight` で再最適化します。`--flatten` は余分なモードに対する明示的な選択肢です。一次鞍点の分類には選択した基準で1本が必要です。 |
| `bond_changes` が空、または想定外の端点 | TSモードとIRCを確認します。意図した極小点同士を結んでいない可能性があります。 |
| R/Pに虚振動が残る | 端点構造とモードを確認し、必要なら端点最適化を厳しくするかIRCを延長します。[freq](freq.md) を参照してください。 |

## 補足

- `tsopt` 単独でのパラメータ調整（`--opt-mode`、`--max-cycles`、Hessian オプション）は [tsopt](tsopt.md) を参照。
- デフォルトの `FiniteDifference` を基準にし、選択した backend/model と対象系で速度・メモリ・結果を検証できた場合にだけ `--hessian-calc-mode Analytical` を明示してください。
- 全オプションを確認するには `pdb2reaction all --help-advanced`。

## 次のステップ

- 複数構造からの MEP 経路: [クイックスタート: `pdb2reaction all`](quickstart-all.md)
- 単一構造からのスキャン駆動: [クイックスタート: `pdb2reaction all --scan-lists`](quickstart-scan.md)
- 全オプションリファレンス: [all](all.md) / [tsopt](tsopt.md) / [irc](irc.md) / [freq](freq.md) / [dft](dft.md)
