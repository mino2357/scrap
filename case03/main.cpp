//-----------------------------------------------
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstdio>   // std::rename
#include "porous_reactor.hpp"

// スナップショット（濃度・温度）を conc.dat / temp.dat に保存
template<class T>
static void dump_snapshot(const State<T>& S) {
    // concentration
    {
        std::ofstream ofs("conc.tmp", std::ios::binary);
        ofs << std::fixed << std::setprecision(15);
        for (int i = 0; i < S.Nx; ++i)
            ofs << S.x[i]  << ' ' << S.cA[i] << ' ' << S.cB[i] << ' ' << S.cC[i] << '\n';
        ofs.flush();
        ofs.close();
        std::rename("conc.tmp", "conc.dat");
    }
    // temperature
    {
        std::ofstream ofs("temp.tmp", std::ios::binary);
        ofs << std::fixed << std::setprecision(15);
        for (int i = 0; i < S.Nx; ++i)
            ofs << S.x[i]  << ' ' << S.Tf[i] << ' ' << S.Ts[i] << '\n';
        ofs.flush();
        ofs.close();
        std::rename("temp.tmp", "temp.dat");
    }
}

// 最終結果をファイルに保存（列: x cA cB cC Tf Ts）
template<class T>
static void save_result(const char* fname, const State<T>& S) {
    std::ofstream ofs(fname, std::ios::binary);
    ofs << std::fixed << std::setprecision(15);
    for (int i = 0; i < S.Nx; ++i)
        ofs << S.x[i]  << ' '
            << S.cA[i] << ' '
            << S.cB[i] << ' '
            << S.cC[i] << ' '
            << S.Tf[i] << ' '
            << S.Ts[i] << '\n';
}

// 1ステップ進める（反応→輸送→境界条件→コミット）
template<class T>
static void advance(const Params<T>& P, State<T>& S, T dt,
                    std::vector<T>& rs, std::vector<T>& Rvol,
                    std::vector<T>& cA_new, std::vector<T>& cB_new, std::vector<T>& cC_new,
                    std::vector<T>& Tf_new, std::vector<T>& Ts_new) {
    compute_reaction(P, S, rs, Rvol);
    update_species_explicit(P, S, Rvol, dt, cA_new, cB_new, cC_new);
    update_energy_explicit (P, S, rs,   dt, Tf_new, Ts_new);
    apply_bc(P, cA_new, cB_new, cC_new, Tf_new, Ts_new);
    S.cA.swap(cA_new); S.cB.swap(cB_new); S.cC.swap(cC_new);
    S.Tf.swap(Tf_new); S.Ts.swap(Ts_new);
}

int main() {
    using T = double;

    // （任意）標準出力のバッファリング制御
    std::ios::sync_with_stdio(false);
    std::cout.setf(std::ios::unitbuf);

    Params<T> P;

    const T   dt      = 1e-5;
    const int nsteps  = 100000;
    const int output_interval = 200; // 可視化スナップショットの出力間隔

    // 初期・入口条件
    State<T> S = make_state(
        P,
        T(0.0), T(0.0), T(0.0), T(300.0), T(300.0),   // 初期条件
        T(1.0), T(0.4), T(0.0), T(350.0)              // 入口（Dirichlet）
    );

    // 作業配列
    std::vector<T> rs, Rvol, cA_new, cB_new, cC_new, Tf_new, Ts_new;

    // 初期スナップショット
    dump_snapshot(S);

    for (int n = 0; n < nsteps; ++n) {
        advance(P, S, dt, rs, Rvol, cA_new, cB_new, cC_new, Tf_new, Ts_new);
        if (n % output_interval == 0)
            dump_snapshot(S);
    }

    save_result("result.dat", S);
    return 0;
}

