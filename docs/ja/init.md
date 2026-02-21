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

## 補足

- 生成される YAML は最小テンプレートです（完全スキーマではありません）。
- 優先順位は `defaults < --config < CLI < --override-yaml` です。

---

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
- [all](all.md) -- `--config` / `--override-yaml` を使った実行
- [YAML リファレンス](yaml-reference.md) -- 設定キーの全体像
