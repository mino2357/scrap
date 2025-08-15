#include "porous_reactor.hpp"
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <limits>
#include <fstream>
#include <string>
#include <cstdio>   // std::rename

int main(){
    // バッファ抑止（パイプでの遅延回避）
    std::ios::sync_with_stdio(false);
    std::cout.setf(std::ios::unitbuf);
    using T=double;
    Params<T> P;
    // グリッドと時間はあなたが管理
    const T dt = 1e-5;
    const int nsteps = 1000000000;
    const int output_interval = 100000;

    /*
    State<T> make_state(
        Params<T>& P, // ← const を外して、ここで inlet BC を書き換え可能にする
        // 初期条件
        T cA0 = T(1), T cB0 = T(1), T cC0 = T(0),
        T Tf0 = T(300), T Ts0 = T(300),
        // 流入条件（入口 Dirichlet）
        T cA_in = T(1), T cB_in = T(0.4), T cC_in = T(0), T Tf_in = T(300)
    )
     */
    State<double> S = make_state(P, 0.0, 0.0, 0.0, 300.0, 300.0, 1.0, 0.8, 0.0, 301.0);

    std::vector<T> rs, Rvol, cA_new, cB_new, cC_new, Tf_new, Ts_new;

    
    auto dump_snapshot = [&](const State<double>& S){
        // conc
        {
            std::ofstream ofs("conc.tmp", std::ios::binary);
            ofs << std::fixed << std::setprecision(15);
            ofs.setf(std::ios::fmtflags(0), std::ios::floatfield);
            ofs.setf(std::ios::fmtflags(0), std::ios::adjustfield);
            for (int i=0;i<S.Nx;i++)
                ofs << S.x[i]  << " " << S.cA[i] << " " << S.cB[i] << " " << S.cC[i] << "\n";
            ofs.flush();
            ofs.close();
            std::rename("conc.tmp","conc.dat");
        }
        // temp
        {
            std::ofstream ofs("temp.tmp", std::ios::binary);
            ofs << std::fixed << std::setprecision(15);
            for (int i=0;i<S.Nx;i++)
                ofs << S.x[i]  << " " << S.Tf[i] << " " << S.Ts[i] << "\n";
            ofs.flush();
            ofs.close();
            std::rename("temp.tmp","temp.dat");
        }
    };

    for(int n=0;n<nsteps;++n){
        compute_reaction(P, S, rs, Rvol);
        update_species_explicit(P, S, Rvol, dt, cA_new, cB_new, cC_new);
        update_energy_explicit (P, S, rs,   dt, Tf_new, Ts_new);
        apply_bc(P, cA_new, cB_new, cC_new, Tf_new, Ts_new);

        // commit
        S.cA.swap(cA_new); S.cB.swap(cB_new); S.cC.swap(cC_new);
        S.Tf.swap(Tf_new); S.Ts.swap(Ts_new);

        if (n % output_interval == 0) {
            dump_snapshot(S);
        }
    }

    return 0;
}