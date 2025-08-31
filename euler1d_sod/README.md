# 1次元オイラー方程式（Sod衝撃管, MUSCL+HLLC）

理想気体の1Dオイラー方程式を，MUSCL 勾配再構築 + HLLC（必要時 HLLE フォールバック）
で空間離散化し，陽解法で時間発展させます。Sod 衝撃管問題を解き，厳密解と
数値解を比較します。

- 方程式: 連続，運動量，エネルギー（比熱比 `gamma`）
- 再構築: MUSCL（minmod スロープ）/ OFF 切替
- リーマン解法: HLLC（安全装置として HLLE にフォールバック）
- 出力: 数値解 `solution.csv`，厳密解 `exact.csv`（t = 0.2）
- 可視化: `sod_compare.gp` → `sod.png`

## ビルド・実行

```bash
make -C euler1d_sod run
# 実行後、gnuplot で重ね描き
gnuplot euler1d_sod/sod_compare.gp
```

`make run` は `main` をビルドして実行し，`solution.csv` と `exact.csv` を生成します。
`gnuplot` スクリプトで密度の重ね描きを `sod.png` として保存します。

## パラメータ
`main.cpp` 内の `Params` でセル数 `N`，CFL 数，最終時刻 `t_end`，比熱比 `gamma` などを
変更できます。`use_muscl=false` にすると一階精度の区分一定からの Riemann 入力になります。

## 出力ファイル
- `solution.csv`: 列は `x, rho, u, p, E`
- `exact.csv`: 厳密解（Sod）`x, rho, u, p`
- `sod.png`: gnuplot による密度プロファイルの比較図
