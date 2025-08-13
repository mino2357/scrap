#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>

//-------------------- パラメータ --------------------
template<class T>
struct Params {
    // Grid
    int Nx = 101;  T L = T(1);

    // Flow / porous
    T u = T(0.1);    // m/s (>0)
    T eps = T(0.5);  // porosity

    // Species diffusion (effective)
    T DA = T(0), DB = T(0), DC = T(0);

    // Thermal (conductivity k, not alpha)
    T kf = T(1.0), ks_eff = T(1.0);     // W/(m*K)
    T rho_f = T(1.0), Cp_f = T(1.0);
    T rho_s = T(1.0), Cp_s = T(1.0);
    T h_sf = T(1.0); // W/(m^3 K)

    // Surface / reaction
    T a_s = T(1.0); // m^2/m^3
    T k0 = T(1.0), Ea = T(1.0), Rg = T(1.0); // Arrhenius
    T dH = T(-1.0); // J/mol (exothermic < 0)
    T gamma_heat_to_fluid = T(0.6); // 0..1

    // Inlet BC (Dirichlet)
    T cA_in = T(1), cB_in = T(0.4), cC_in = T(0), Tf_in = T(300);
};

//-------------------- 状態 --------------------
template<class T>
struct State {
    int Nx{}; T dx{};
    std::vector<T> x, cA, cB, cC, Tf, Ts;
};

//-------------------- 数学ユーティリティ --------------------
template<class T>
inline T arrhenius(T k0, T Ea, T Rg, T Ts) {
    const T Tmin = T(1);
    return k0 * std::exp(-Ea / (Rg * std::max(Ts, Tmin)));
}

// S字型温度因子: 0.5*(1 + (T-300)/sqrt(10000+(T-300)^2))
template<class T>
inline T temp_factor(T Ts) {
    const T T0 = T(300);
    const T w2 = T(10000);
    const T dT = Ts - T0;
    return T(0.5) * (T(1) + dT / std::sqrt(w2 + dT * dT));
}

// upwind (u>0)
template<class T>
inline T upwind_grad_posU(const std::vector<T>& q, T dx, int i) {
    return (i == 0) ? T(0) : (q[i] - q[i - 1]) / dx;
}

// Laplacian with one-sided Neumann(0) at both ends
template<class T>
inline T lap_central(const std::vector<T>& q, T dx, int i) {
    const int N = static_cast<int>(q.size());
    if (i == 0)      return (T(2)*(q[1] - q[0])) / (dx * dx);
    if (i == N - 1)  return (T(2)*(q[N - 2] - q[N - 1])) / (dx * dx);
    return (q[i + 1] - T(2) * q[i] + q[i - 1]) / (dx * dx);
}

//-------------------- 初期化 --------------------
// Params を参照で受け取り，入口条件を上書き可能にする
template<class T>
State<T> make_state(Params<T>& P,
                    // 初期条件
                    T cA0 = T(1), T cB0 = T(1), T cC0 = T(0),
                    T Tf0 = T(300), T Ts0 = T(300),
                    // 入口条件（Dirichlet）
                    T cA_in = T(1), T cB_in = T(0.4), T cC_in = T(0), T Tf_in = T(300)) {
    // inlet BC を Params に反映
    P.cA_in = cA_in; P.cB_in = cB_in; P.cC_in = cC_in; P.Tf_in = Tf_in;

    State<T> S;
    S.Nx = P.Nx;
    assert(S.Nx >= 2);
    S.dx = P.L / T(P.Nx - 1);
    S.x.resize(P.Nx);
    S.cA.assign(P.Nx, cA0);
    S.cB.assign(P.Nx, cB0);
    S.cC.assign(P.Nx, cC0);
    S.Tf.assign(P.Nx, Tf0);
    S.Ts.assign(P.Nx, Ts0);
    for (int i = 0; i < P.Nx; ++i) S.x[i] = T(i) * S.dx;
    return S;
}

//-------------------- 反応レート（表面→体積） --------------------
template<class T>
void compute_reaction(const Params<T>& P, const State<T>& S,
                      std::vector<T>& rs, std::vector<T>& Rvol) {
    const int N = S.Nx;
    rs.resize(N);
    Rvol.resize(N);
    for (int i = 0; i < N; ++i) {
        const T kA = arrhenius(P.k0, P.Ea, P.Rg, S.Ts[i]);
        const T cA = std::max(S.cA[i], T(0));
        const T cB = std::max(S.cB[i], T(0));
        const T fT = temp_factor(S.Ts[i]);
        rs[i]   = kA * cA * cB * fT;  // mol/(m^2 s)
        Rvol[i] = P.a_s * rs[i];      // mol/(m^3 s)
    }
}

//-------------------- 物質更新（完全陽） --------------------
template<class T>
void update_species_explicit(const Params<T>& P, const State<T>& S,
                             const std::vector<T>& Rvol, T dt,
                             std::vector<T>& cA_new,
                             std::vector<T>& cB_new,
                             std::vector<T>& cC_new) {
    const int N = S.Nx;
    const T dx = S.dx;
    cA_new = S.cA;
    cB_new = S.cB;
    cC_new = S.cC;
    const T inv_eps = T(1) / P.eps;

    for (int i = 0; i < N; ++i) {
        const T gcA = upwind_grad_posU(S.cA, dx, i);
        const T gcB = upwind_grad_posU(S.cB, dx, i);
        const T gcC = upwind_grad_posU(S.cC, dx, i);
        const T lcA = lap_central(S.cA, dx, i);
        const T lcB = lap_central(S.cB, dx, i);
        const T lcC = lap_central(S.cC, dx, i);

        cA_new[i] = S.cA[i] + dt * (-(P.u * inv_eps) * gcA + (P.DA * inv_eps) * lcA - inv_eps * Rvol[i]);
        cB_new[i] = S.cB[i] + dt * (-(P.u * inv_eps) * gcB + (P.DB * inv_eps) * lcB - inv_eps * Rvol[i]);
        cC_new[i] = S.cC[i] + dt * (-(P.u * inv_eps) * gcC + (P.DC * inv_eps) * lcC + inv_eps * Rvol[i]);
    }
}

//-------------------- エネルギー更新（完全陽、反応熱分配あり） --------------------
template<class T>
void update_energy_explicit(const Params<T>& P, const State<T>& S,
                            const std::vector<T>& rs, T dt,
                            std::vector<T>& Tf_new, std::vector<T>& Ts_new) {
    const int N = S.Nx;
    const T dx = S.dx;
    Tf_new = S.Tf;
    Ts_new = S.Ts;
    const T one = T(1);

    const T fluid_cap = P.eps * P.rho_f * P.Cp_f;             // ε ρ_f Cp_f
    const T solid_cap = (one - P.eps) * P.rho_s * P.Cp_s;      // (1-ε) ρ_s Cp_s

    for (int i = 0; i < N; ++i) {
        const T gTf = upwind_grad_posU(S.Tf, dx, i);
        const T lTf = lap_central(S.Tf, dx, i);
        const T lTs = lap_central(S.Ts, dx, i);

        // 固⇔流の熱交換（+: 固体→流体）
        const T q_fs = P.h_sf * (S.Ts[i] - S.Tf[i]); // W/m^3

        // 反応熱（総量）を分配
        const T q_rx_total = P.a_s * (-P.dH) * rs[i];      // W/m^3
        const T q_rx_f = P.gamma_heat_to_fluid * q_rx_total;
        const T q_rx_s = (one - P.gamma_heat_to_fluid) * q_rx_total;

        // 流体温度：dTf/dt = - (u/ε) * ∂Tf/∂x + (k_f/(ε ρ_f Cp_f)) ∂^2Tf/∂x^2 + ...
        Tf_new[i] = S.Tf[i] + dt * (
            -(P.u / P.eps) * gTf
            + (P.kf / fluid_cap) * lTf
            + q_fs / fluid_cap
            + q_rx_f / fluid_cap
        );

        // 固体温度：dTs/dt = (k_s^eff/((1-ε) ρ_s Cp_s)) ∂^2Ts/∂x^2 + ...
        Ts_new[i] = S.Ts[i] + dt * (
            (P.ks_eff / solid_cap) * lTs
            - q_fs / solid_cap
            + q_rx_s / solid_cap
        );
    }
}

//-------------------- 境界条件 --------------------
template<class T>
void apply_bc(const Params<T>& P,
              std::vector<T>& cA, std::vector<T>& cB, std::vector<T>& cC,
              std::vector<T>& Tf, std::vector<T>& Ts) {
    const int N = static_cast<int>(cA.size());
    assert(N >= 2);
    // Inlet: Dirichlet（流体）、固体は断熱（Neumann0）
    cA[0] = P.cA_in; cB[0] = P.cB_in; cC[0] = P.cC_in; Tf[0] = P.Tf_in;
    Ts[0] = Ts[1]; // insulated

    // Outlet: Neumann(0) ＝ コピーアウト
    cA[N - 1] = cA[N - 2];
    cB[N - 1] = cB[N - 2];
    cC[N - 1] = cC[N - 2];
    Tf[N - 1] = Tf[N - 2];
    Ts[N - 1] = Ts[N - 2];
}
