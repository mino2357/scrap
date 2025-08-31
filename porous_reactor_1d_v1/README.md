# 1次元 多孔質触媒反応器（v1: α形式 + 体積基準熱交換）

この v1 は、熱拡散率（α = k/(ρCp)）で記述し、固体–流体の熱交換係数 `h_sf` を
体積基準 [W/(m³·K)] として直接用いるバリアントです。

このコードは、**1次元の多孔質触媒反応器**における流体相と固体相の**化学種輸送と熱輸送**を連成して計算します。  
支配方程式には以下を含みます：

- **移流**（流速 $u$ 一定、正方向：左 → 右）
- **有効拡散**（化学種と熱の拡散）
- **不均一触媒反応** $A + B \to C$
- **発熱反応／吸熱反応による熱生成**
- **流体相と固体相間の熱交換**

---

## 支配方程式

記号の定義：
- $c_A, c_B, c_C$: 流体相における化学種A, B, Cのモル濃度 [mol/m³]
- $T_f$: 流体温度 [K]
- $T_s$: 固体温度 [K]
- $\varepsilon$: 空隙率 $(-)$
- $a_s$: 反応器単位体積あたりの触媒表面積 [m²/m³]
- $u$: 見かけ流速 [m/s]
- $D_A, D_B, D_C$: 有効拡散係数 [m²/s]
- $k_f$: 流体相の有効熱拡散率（$k_f/(\rho_f C_{p,f})$）[m²/s]
- $k_s^{\mathrm{eff}}$: 固体相の有効熱拡散率 [m²/s]
- $h_{sf}$: 固体と流体間の体積基準熱伝達係数 [W/(m³·K)]
- $k_0, E_a, R_g$: アレニウス式パラメータ
- $\Delta H$: 反応エンタルピー変化 [J/mol]（発熱反応なら負）
- $\gamma$: 反応熱のうち流体相に分配される割合
- $r_s$: 表面反応速度 [mol/(m²·s)]
- $R_{\mathrm{vol}}$: 体積基準反応速度 [mol/(m³·s)]

---

### 化学種輸送式（流体相, $A$, $B$, $C$）

$$
\varepsilon \frac{\partial c_A}{\partial t} + u \frac{\partial c_A}{\partial x}
= D_A \frac{\partial^2 c_A}{\partial x^2} - R_{\mathrm{vol}}
$$

$$
\varepsilon \frac{\partial c_B}{\partial t} + u \frac{\partial c_B}{\partial x}
= D_B \frac{\partial^2 c_B}{\partial x^2} - R_{\mathrm{vol}}
$$

$$
\varepsilon \frac{\partial c_C}{\partial t} + u \frac{\partial c_C}{\partial x}
= D_C \frac{\partial^2 c_C}{\partial x^2} + R_{\mathrm{vol}}
$$

---

### 反応速度
**アレニウス式**と温度因子 $f_T(T_s)$ を使用：

$$
r_s(T_s, c_A, c_B) = k_0 \exp\left(-\frac{E_a}{R_g \, T_s}\right) \cdot c_A \cdot c_B \cdot f_T(T_s)
$$

**温度因子**（300K付近でS字型変化）：

$$
f_T(T) = \frac12 \left[ 1 + \frac{T - 300}{\sqrt{10000 + (T - 300)^2}} \right]
$$

**体積基準反応速度**：

$$
R_{\mathrm{vol}} = a_s \, r_s
$$

---

### エネルギー輸送式

**流体相:**

$$
\varepsilon \rho_f C_{p,f} \left( \frac{\partial T_f}{\partial t} + u \frac{\partial T_f}{\partial x} \right)
= \varepsilon k_f \frac{\partial^2 T_f}{\partial x^2}
+ h_{sf}(T_s - T_f)
+ \gamma \, q_{rx}
$$

**固体相:**

$$
(1-\varepsilon)\rho_s C_{p,s} \frac{\partial T_s}{\partial t}
= (1-\varepsilon) k_s^{\mathrm{eff}} \frac{\partial^2 T_s}{\partial x^2}
+ h_{sf}(T_f - T_s)
+ (1-\gamma) \, q_{rx}
$$

ここで

$$
q_{rx} = a_s \, (-\Delta H) \, r_s \quad [\mathrm{W/m^3}]
$$
$\gamma$ は反応熱の流体相への分配割合、$(1-\gamma)$ は固体相への分配割合。

---

## 境界条件

- **入口 ($x = 0$)**:
  - 流体濃度・温度: **Dirichlet条件**
    $$
    c_A(0,t) = c_{A,\mathrm{in}}, \quad
    c_B(0,t) = c_{B,\mathrm{in}}, \quad
    c_C(0,t) = c_{C,\mathrm{in}}, \quad
    T_f(0,t) = T_{f,\mathrm{in}}
    $$
  - 固体温度: **Neumannゼロ勾配**（断熱条件）
    $$
    \left.\frac{\partial T_s}{\partial x}\right|_{x=0} = 0
    $$

- **出口 ($x = L$)**:
  - 全変数: **Neumannゼロ勾配**
    $$
    \frac{\partial}{\partial x} (\cdot) = 0
    $$

---

## 離散化

- **移流**: 正の $u$ に対して一次精度風上差分
- **拡散**: 中心差分
- **時間積分**: 陽的オイラー法
- **格子**: 一様格子、区間 $[0,L]$ を $N_x$ 分割

---

## コード構成

- `Params<T>` : モデルパラメータ（形状、輸送係数、反応条件、入口BCなど）
- `State<T>` : シミュレーション状態（x座標、濃度、温度）
- `make_state` : 初期条件および入口条件を設定して状態を初期化
- `compute_reaction` : $r_s$ および $R_{\mathrm{vol}}$ を計算
- `update_species_explicit` : $c_A$, $c_B$, $c_C$ の陽的更新
- `update_energy_explicit` : $T_f$, $T_s$ の陽的更新
- `apply_bc` : 入口・出口の境界条件を適用

---

## `make_state` による初期化

```cpp
template<class T>
State<T> make_state(
    Params<T>& P,
    // 初期条件
    T cA0 = T(1), T cB0 = T(1), T cC0 = T(0),
    T Tf0 = T(300), T Ts0 = T(300),
    // 入口条件（Dirichlet）
    T cA_in = T(1), T cB_in = T(0.4), T cC_in = T(0), T Tf_in = T(300)
);

```

---

## ビルド・実行と可視化

```bash
make -C porous_reactor_1d_v1 run   # 実行（conc.dat / temp.dat を出力）
gnuplot -persist porous_reactor_1d_v1/plot_live.gp  # ライブ描画
```

- `conc.dat`: 列は `x cA cB cC`
- `temp.dat`: 列は `x Tf Ts`
- `plot_live.gp`: 2 画面（濃度・温度）を定期リロードで更新描画
