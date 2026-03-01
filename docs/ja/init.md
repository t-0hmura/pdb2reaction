# `init`

## 概要

`pdb2reaction init` は `pdb2reaction all` 用のスターター YAML 設定ファイルを生成します。

設定ファイル中心で再現性の高い実行を始めるときに使用します:

```bash
pdb2reaction init --out pdb2reaction_all.config.yaml
pdb2reaction all --config pdb2reaction_all.config.yaml --dry-run
```

## 使い方

```bash
pdb2reaction init [--out PATH] [--force]
```

## オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-o, --out PATH` | 生成先 YAML パス。 | `pdb2reaction_all.config.yaml` |
| `--force/--no-force` | 既存ファイルを上書き。 | `False` |

## 生成されるテンプレート

生成されるファイルには次のスターター設定が含まれます:

```yaml
# Starter config for `pdb2reaction all`

extract:
  radius: 2.6
  radius_het2het: 0.0

calc:
  workers: 1
  workers_per_node: 1

path_search:
  mep_mode: gsm
  max_nodes: 10
  max_cycles: 300

scan:
  max_step_size: 0.2
  bias_k: 300.0
  relax_max_cycles: 10000

tsopt:
  max_cycles: 10000

freq:
  max_write: 10
  amplitude_ang: 0.8
  n_frames: 20
  sort: value
  temperature: 298.15
  pressure: 1.0

dft:
  func_basis: wb97m-v/def2-tzvpd
  max_cycle: 100
  conv_tol: 1.0e-9
  grid_level: 3
```

## 次のステップ

1. 生成された YAML を系に合わせて編集します（電荷、多重度、基質など）。
2. ドライランで検証します: `pdb2reaction all --config pdb2reaction_all.config.yaml --dry-run -i R.pdb P.pdb -c "SUBSTRATE"`。
3. 問題なければ `--dry-run` を外して本実行します。
4. 詳細な調整は [YAML リファレンス](yaml_reference.md) を参照してください。

## 補足

- 生成ファイルは出発点であり、完全なスキーマ定義ではありません。

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
- [YAML リファレンス](yaml_reference.md) -- 設定キーの全体像
