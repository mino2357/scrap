# ---- plot.gp ----
# conc.dat: x cA cB cC
# temp.dat: x Tf Ts

if (!exists("delay")) delay = 0.05   # 秒。必要に応じて変更

set term qt size 1100,500
set grid
set key outside right

# 2パネル描画をマクロにまとめる（毎回ファイルを読み直す）
DRAW = \
"set multiplot layout 1,2 title 'Porous reactor (live file polling)';" \
. "set xlabel 'x [m]'; set ylabel 'Concentration [mol/m^3]';" \
. "plot 'conc.dat' using 1:2 w l t 'cA', '' using 1:3 w l t 'cB', '' using 1:4 w l t 'cC';" \
. "set xlabel 'x [m]'; set ylabel 'Temperature [K]';" \
. "plot 'temp.dat' using 1:2 w l t 'Tf', '' using 1:3 w l t 'Ts';" \
. "unset multiplot;"

# 初回描画
eval DRAW

# live ループ（ファイルを毎回読み直して再描画）
do for [i=1:1e9] {
    pause delay
    eval DRAW
}

# ウィンドウを閉じるまで続けたい場合は下の一行でもOK
# pause mouse close