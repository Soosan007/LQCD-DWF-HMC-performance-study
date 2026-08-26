# コード変更一覧

行番号はこのリポジトリに保存したファイルを基準にしています。今後編集するとずれるため、関数名・マクロ名・ログ識別子でも検索してください。

## CGソルバー

### `modified_source/lib_alt/Solver/asolver_CG.h`

- 37–41行付近: 前回解キャッシュと有効性・mode管理
- 45–66行付近: CGの呼出回数と内部時間を保存する `m_profile_*`
- 115行付近: 初期解と初期解modeを受け取る `solve_CG_init`

### `modified_source/lib_alt/Solver/asolver_CG-tmpl.h`

- `set_parameters`: YAMLの `reuse_previous_solution` を取得
- `solve`: CG呼出回数、反復回数、収束、最終残差、各処理時間を記録
- `solve_CG_init`: `RHS`、`GIVEN`、`ZERO` の初期解を処理
- `solve_CG_step`: `mult`、`dot`、`axpy`、残差更新、`aypx` を計測
- `USE_CG_FUSED_AXPY_NORM2`: `axpy + norm2` と `axpy_norm2` を切替
- `step_residual_block_sec`: 融合あり／なしを同じ外側タイマーで計測

主なログ:

```text
CG_PROFILE
SOLVER_MULT_TRACE
previous_guess_used=
iterations=
final_relative_residual=
step_residual_block_sec=
```

## QXS Field

### `modified_source/lib_alt_QXS/Field/afield.h`

- `axpy_norm2`: ベクトル更新と更新後ノルムを融合
- `copy_axpy`: `z = x + a*y` を1走査で実行

### `modified_source/lib_alt_QXS/Field/afield-tmpl.h`

- `AField::copy_axpy`: copyとaxpyを一つのACLE/SVEループに統合
- `AField::axpy_norm2`: axpy、保存、ノルム加算を一つのループに統合

### `modified_source/lib_alt_QXS/Field/afield-inc.h`

- ソルバー／演算子側から呼ぶ `axpy_norm2` と `copy_axpy` のラッパー

## DdagD

### `modified_source/lib_alt_QXS/Fopr/afopr_Domainwall_5din_eo.h`

- `m_dagd_profile_calls`
- `m_dagd_profile_sec[15]`

### `modified_source/lib_alt_QXS/Fopr/afopr_Domainwall_5din_eo-tmpl.h`

- DdagD全体と15内部区間の時間を記録
- `USE_DAGD_FUSED_COPY_AXPY` でcombine処理を切替
- `timing_scheme=single_outer_block` をログへ追加

主なログ:

```text
DAGD_DETAIL_PROFILE
DAGD_FUSION_PROFILE
fused_copy_axpy=
dagd_total_sec=
combine_fraction_percent=
```

## ビルド切替

- `build/Makefile_fair_timing_no_fusion`: 詳細計測のみ
- `build/Makefile_fair_timing_all_fusion`: CG融合、DdagD融合、詳細計測

両版は同じソーススナップショットとコンパイラ設定からビルドし、マクロだけを変える設計です。
