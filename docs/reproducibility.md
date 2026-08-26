# 再現方法

## 1. 変更ソースを適用

`modified_source/` の相対パスを、対応するLQCD-DWF-HMCソースツリーへ配置します。適用前に元ファイルを保存し、差分を確認してください。

## 2. 2種類をビルド

- 融合なし: `build/Makefile_fair_timing_no_fusion`
- 融合あり: `build/Makefile_fair_timing_all_fusion`

両方を同一モジュール、同一コンパイラ、同一ソーススナップショットからビルドします。

## 3. ペア実験を作成

`scripts/run_pair_template.pjm` のプレースホルダーを置換し、同じノード割当内で融合なしと融合ありを順番に実行します。奇数・偶数ペアで順序を反転してください。

## 4. ログを確認

```bash
grep 'CG_PROFILE' hmc.log
grep 'DAGD_FUSION_PROFILE' hmc.log
grep 'DAGD_DETAIL_PROFILE' hmc.log
grep 'H(diff)' hmc.log
```

次を確認します。

- `converged=1`
- `timing_scheme=single_outer_block`
- 反復回数が比較条件間で一致
- `H(diff)` が一致
- 実行ファイルと入力ファイルのSHA256を保存

## 5. 集約

ペアディレクトリとログを `scripts/aggregate_fair_combined.py` と同じ配置に置き、実行します。

```bash
python3 scripts/aggregate_fair_combined.py
```

元のスクリプトは10個の `pair01`–`pair10` を想定しています。

## 注意

- `Main`、trajectory、CG、局所融合区間は包含範囲が異なります。
- 内側の削減時間を、そのまま外側のMain時間差として説明しないでください。
- 入力配位はこのリポジトリに含まれません。
