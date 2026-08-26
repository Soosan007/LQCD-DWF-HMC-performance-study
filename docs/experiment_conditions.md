# 実験条件

## 富岳: 最終公平比較

| 項目 | 条件 |
|---|---|
| ノード数 | 1ノード/ペア |
| ペア数 | 10 |
| 実行数 | 融合なし10、融合あり10 |
| MPI | 4プロセス |
| OpenMP | 12スレッド/MPIプロセス |
| プロセス格子 | 1×1×2×2 |
| 格子 | 8×8×8×8 |
| 第5次元 | Ls=12 |
| CPU周波数 | 2200 MHz |
| eco_state | 2 |
| ソルバー | CG |
| 前回解再利用 | false |
| 比較順序 | 奇数ペア: no_fusion→fusion、偶数ペア: fusion→no_fusion |
| 計測方式 | `single_outer_block` |

共通の入力配位、YAML、合理近似パラメータを使用し、各入力と実行ファイルのSHA256を実行前に記録しました。入力配位そのものは容量と配布範囲の都合からGitHubには含めません。

### 比較するマクロ

融合なし:

```text
USE_DAGD_DETAIL_PROFILE
```

融合あり:

```text
USE_CG_FUSED_AXPY_NORM2
USE_DAGD_FUSED_COPY_AXPY
USE_DAGD_DETAIL_PROFILE
```

## Genoa: OpenMPスケーリング

| 項目 | 条件 |
|---|---|
| CPU | AMD EPYC 9684X, 96 cores |
| ノード数 | 1 |
| MPI | 1プロセス |
| OpenMP | 1, 2, 4, 8, 16, 32, 64, 96 |
| 実装 | x86 QXS general |
| 格子 | 8×8×8×8 |
| 第5次元 | Ls=12 |
| trajectory | 1 |

Genoaの詳細結果は `data/genoa/` のワークブックに保存しています。
