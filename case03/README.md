# 1D Porous Reactor Model with Heat and Reaction

## 概要
このコードは，多孔質触媒反応における1次元の移流・拡散・反応・熱移動をモデル化したものである．
流体相と固相を別温度場として扱う **二温度モデル** を採用し，触媒表面での発熱反応による温度上昇と熱伝達を計算します．

反応系：

```math
A + B \longrightarrow C
```

反応速度はアレニウス型で，温度依存の活性化関数 $f_T(T)$ を組み込みます．

---

## 支配方程式（有次元）

### 種輸送（流体相）

```math
\varepsilon\frac{\partial c_A}{\partial t} + u\frac{\partial c_A}{\partial x} = D_A\frac{\partial^2 c_A}{\partial x^2} - R_\mathrm{vol}
```

```math
\varepsilon\frac{\partial c_B}{\partial t} + u\frac{\partial c_B}{\partial x} = D_B\frac{\partial^2 c_B}{\partial x^2} - R_\mathrm{vol}
```

```math
\varepsilon\frac{\partial c_C}{\partial t} + u\frac{\partial c_C}{\partial x} = D_C\frac{\partial^2 c_C}{\partial x^2} + R_\mathrm{vol}
```

### エネルギー

流体：

```math
\varepsilon\rho_f C_{p,f} \frac{\partial T_f}{\partial t} + \rho_f C_{p,f} u \frac{\partial T_f}{\partial x} = k_f\frac{\partial^2 T_f}{\partial x^2} + h_{sf}(T_s-T_f) + \gamma\,a_s(-\Delta H)r_s
```

固体：

```math
(1-\varepsilon)\rho_s C_{p,s} \frac{\partial T_s}{\partial t} = k_s^{\mathrm{eff}}\frac{\partial^2 T_s}{\partial x^2} + h_{sf}(T_f-T_s) + (1-\gamma)\,a_s(-\Delta H)r_s
```

### 反応速度

```math
r_s = k_0 \exp\!\left(-\frac{E_a}{R_g T_s}\right) c_A c_B\, f_T(T_s)
```

ここで

```math
f_T(T_s) = \frac12\left[1 + \frac{T_s - 300}{\sqrt{10000 + (T_s-300)^2}}\right]
```

$R_\mathrm{vol} = a_s r_s$

---

## 無次元化

### スケールの選択
- 長さ： $L$
- 速度： $u>0$
- 時間： $t_\mathrm{ref}=L/u$（移流時間）
- 濃度スケール： $C_0$（例：入口 $c_{A,\mathrm{in}}$）
- 温度：基準 $T_0=300\ \mathrm{K}$，温度差 $\Delta T=100\ \mathrm{K}$

### 無次元変数

```math
\hat{x}=\frac{x}{L},\quad
\hat{t}=\frac{u\,t}{L},\quad
\hat{c}_i=\frac{c_i}{C_0},\quad
\hat{T}_f=\frac{T_f-T_0}{\Delta T},\quad
\hat{T}_s=\frac{T_s-T_0}{\Delta T}.
```

### 無次元パラメータ

- ペクレ数：

```math
\mathrm{Pe}_A=\frac{uL}{D_A},\quad
\mathrm{Pe}_B=\frac{uL}{D_B},\quad
\mathrm{Pe}_C=\frac{uL}{D_C}.
```

- ダムケラー数：

```math
\mathrm{Da}=\frac{a_s\,k_0\,C_0\,L}{\varepsilon\,u}.
```

- 熱ペクレ数：

```math
\mathrm{Pe}_{Tf}=\frac{uL}{k_f},\quad
\mathrm{Pe}_{Ts}=\frac{uL}{k_s^{\mathrm{eff}}}.
```

- 無次元熱交換数：

```math
H_f=\frac{h_{sf}L}{\varepsilon\,\rho_f C_{p,f}\,u},\quad
H_s=\frac{h_{sf}L}{(1-\varepsilon)\,\rho_s C_{p,s}\,u}.
```

- 無次元反応熱：

```math
\chi_f=\frac{(-\Delta H)C_0}{\rho_f C_{p,f}\,\Delta T},\quad
\chi_s=\frac{(-\Delta H)C_0}{\rho_s C_{p,s}\,\Delta T}.
```

### 無次元方程式

#### 種輸送

```math
\frac{\partial \hat{c}_A}{\partial \hat{t}}+\frac{\partial \hat{c}_A}{\partial \hat{x}} =\frac{1}{\mathrm{Pe}_A}\frac{\partial^2 \hat{c}_A}{\partial \hat{x}^2} -\mathrm{Da}\,\mathrm{e}^{-\Theta/(1+\delta\hat{T}_s)}\,\hat{c}_A\hat{c}_B\,f_T(\hat{T}_s)
```

```math
\frac{\partial \hat{c}_B}{\partial \hat{t}}+\frac{\partial \hat{c}_B}{\partial \hat{x}} =\frac{1}{\mathrm{Pe}_B}\frac{\partial^2 \hat{c}_B}{\partial \hat{x}^2} -\mathrm{Da}\,\mathrm{e}^{-\Theta/(1+\delta\hat{T}_s)}\,\hat{c}_B\hat{c}_A\,f_T(\hat{T}_s)
```

```math
\frac{\partial \hat{c}_C}{\partial \hat{t}}+\frac{\partial \hat{c}_C}{\partial \hat{x}} =\frac{1}{\mathrm{Pe}_C}\frac{\partial^2 \hat{c}_C}{\partial \hat{x}^2} +\mathrm{Da}\,\mathrm{e}^{-\Theta/(1+\delta\hat{T}_s)}\,\hat{c}_A\hat{c}_B\,f_T(\hat{T}_s)
```

#### エネルギー

流体：

```math
\frac{\partial \hat{T}_f}{\partial \hat{t}}+\frac{\partial \hat{T}_f}{\partial \hat{x}} =\frac{1}{\mathrm{Pe}_{Tf}}\frac{\partial^2 \hat{T}_f}{\partial \hat{x}^2} +H_f(\hat{T}_s-\hat{T}_f) +\gamma\,\mathrm{Da}\,\chi_f\,\mathrm{e}^{-\Theta/(1+\delta\hat{T}_s)}\,\hat{c}_A\hat{c}_B\,f_T(\hat{T}_s)
```

固体：

```math
\frac{\partial \hat{T}_s}{\partial \hat{t}} =\frac{1}{\mathrm{Pe}_{Ts}}\frac{\partial^2 \hat{T}_s}{\partial \hat{x}^2} +H_s(\hat{T}_f-\hat{T}_s) +(1-\gamma)\,\mathrm{Da}\,\chi_s\,\mathrm{e}^{-\Theta/(1+\delta\hat{T}_s)}\,\hat{c}_A\hat{c}_B\,f_T(\hat{T}_s)
```

---

## 使用方法

### コンパイル
```bash
g++ -std=c++17 -O2 -Wall -Wextra main.cpp -o sim
