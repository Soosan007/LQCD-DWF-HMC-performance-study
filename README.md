# LQCD-DWF-HMC Performance Study

富岳およびGenoa上で、Domain-Wall Fermion HMC（DWF-HMC）の並列性能と、CG/DdagD内部のメモリ走査融合を評価した研究記録です。

このリポジトリでは、次を再現可能な形で保存します。

- CGソルバーの呼出回数・反復回数・内部時間の計測
- 前回解を初期解として再利用する実装
- CG内の `axpy + norm2` 融合
- DdagD内の `copy + axpy` 融合
- 融合あり／なしを同一区間で測る公平計測
- 富岳の対応あり10組比較とGenoa OpenMPスケーリングの要約

## 主要結果

富岳の最終公平計測では、1ノード、4 MPI × 12 OpenMP、格子 `8×8×8×8` の条件で、融合なし10回と融合あり10回を同一ノード上で対応させて比較しました。

| 指標 | 融合なし | 融合あり | 短縮率 |
|---|---:|---:|---:|
| CG残差更新ブロック | 0.104544 s | 0.090146 s | 13.77% |
| CG全体 | 5.715252 s | 5.679530 s | 0.63% |
| DdagD combine | 0.134884 s | 0.093181 s | 30.92% |
| DdagD全体 | 7.987897 s | 7.955877 s | 0.40% |
| HMC Main | 23.232 s | 23.131 s | 0.43% |

局所的な融合対象では短縮を確認できました。一方、trajectory、Main、外部wall clockの対応差は測定変動より小さく、この10組だけからHMC全体が高速化したとは断定しません。

数値の妥当性確認では、両条件でCG反復回数が5549、`H(diff)=0.00083995` と一致しています。

詳細は [docs/results.md](docs/results.md) を参照してください。

## リポジトリ構成

```text
.
├── README.md
├── LICENSE
├── NOTICE.md
├── docs/
│   ├── research_overview.md
│   ├── code_changes.md
│   ├── experiment_conditions.md
│   ├── results.md
│   └── reproducibility.md
├── modified_source/
│   ├── lib_alt/Solver/
│   └── lib_alt_QXS/{Field,Fopr}/
├── build/
├── scripts/
├── presentation/
│   ├── research_presentation_ja.pptx
│   └── research_presentation_ja.pdf
└── data/
    ├── fugaku/
    ├── genoa/
    └── summary/              # 発表用の集計Excel（01～10）
```

## 集計データと発表資料

発表に使用する主要な集計結果は [`data/summary`](data/summary/) に保存しています。Genoa、富岳のACLE/General比較、MPI/OpenMP構成、DdagD・CG内部、関数融合、通信オーバーラップ、初期解、GPU、およびDWF Multigridの結果を、番号順のExcelファイルとして収録しています。

発表スライドとPDF版は [`presentation`](presentation/) に保存しています。

## 変更の有効化

主要なコンパイルマクロは次のとおりです。

| マクロ | 役割 |
|---|---|
| `USE_CG_FUSED_AXPY_NORM2` | CGの残差更新とノルム計算を融合 |
| `USE_DAGD_FUSED_COPY_AXPY` | DdagDのcopy+axpyを融合 |
| `USE_DAGD_DETAIL_PROFILE` | DdagD内部15区間を計測 |
| `USE_SOLVER_MULT_TRACE` | ソルバーの行列演算呼出情報を出力 |

融合なし／ありの完全な切替は [build](build/) に保存しています。

## データ方針

GitHubには、変更ソース、実験条件、集約スクリプト、要約CSV・Excel、および発表資料を置きます。以下は含めません。

- 入力配位ファイル
- 生ログ一式
- 実行ファイル、静的ライブラリ
- tar/zipバックアップ
- 個人PCや計算機アカウントの絶対パス

## ライセンス

変更ソースはLQCD-DWF-HMC/Bridge++系コードを基にしており、元コードと同様にGNU GPL v3として扱います。詳細は [LICENSE](LICENSE) と [NOTICE.md](NOTICE.md) を参照してください。

このリポジトリは個人の研究記録であり、理化学研究所や富岳運用組織の公式リポジトリではありません。
