#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstdio>   // std::rename
#include "porous_reactor.hpp"

int main(){
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

    // gnuplot, conc.dat と temp.dat を原子的に更新
    auto dump_snapshot = [&](const State<T>& Sref){
        {
            std::ofstream ofs("conc.tmp", std::ios::binary);
            ofs << std::fixed << std::setprecision(15);
            for(int i=0;i<Sref.Nx;i++){
                ofs << Sref.x[i]  << " "
                    << Sref.cA[i] << " "
                    << Sref.cB[i] << " "
                    << Sref.cC[i] << "\n";
            }
            ofs.flush();
            ofs.close();
            std::rename("conc.tmp","conc.dat");
        }
        {
            std::ofstream ofs("temp.tmp", std::ios::binary);
            ofs << std::fixed << std::setprecision(15);
            for(int i=0;i<Sref.Nx;i++){
                ofs << Sref.x[i] << " "
                    << Sref.Tf[i] << " "
                    << Sref.Ts[i] << "\n";
            }
            ofs.flush();
            ofs.close();
            std::rename("temp.tmp","temp.dat");
        }
    };

    // 初期スナップショット
    dump_snapshot(S);

    for(int n=0; n<nsteps; ++n){
        // (1) 反応
        compute_reaction(P, S, rs, Rvol);

        // (2) 種・エネルギー更新（完全陽）
        update_species_explicit(P, S, Rvol, dt, cA_new, cB_new, cC_new);
        update_energy_explicit (P, S, rs,   dt, Tf_new, Ts_new);

        // (3) 境界条件
        apply_bc(P, cA_new, cB_new, cC_new, Tf_new, Ts_new);

        // (4) commit
        S.cA.swap(cA_new); S.cB.swap(cB_new); S.cC.swap(cC_new);
        S.Tf.swap(Tf_new); S.Ts.swap(Ts_new);

        // (5) 可視化スナップショット
        if(n % output_interval == 0){
            dump_snapshot(S);
        }
    }

    // 最終結果（列: x cA cB cC Tf Ts）を保存（15桁）
    {
        std::ofstream ofs("result.dat", std::ios::binary);
        ofs << std::fixed << std::setprecision(15);
        for(int i=0;i<S.Nx;i++){
            ofs << S.x[i]  << " "
                << S.cA[i] << " "
                << S.cB[i] << " "
                << S.cC[i] << " "
                << S.Tf[i] << " "
                << S.Ts[i] << "\n";
        }
    }

    return 0;
}
