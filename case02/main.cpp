#include "porous_reactor.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstdio>   // std::rename

// スナップショットを conc.dat, temp.dat に保存
static void dump_snapshot(const State<double>& S) {
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

    State<double> S = make_state(P, 0.0, 0.0, 0.0, 300.0, 300.0, 1.0, 0.8, 0.0, 301.0);

    std::vector<T> rs, Rvol, cA_new, cB_new, cC_new, Tf_new, Ts_new;

    // 初期状態を出力
    dump_snapshot(S);

    for(int n=0;n<nsteps;++n){
        compute_reaction(P, S, rs, Rvol);
        update_species_explicit(P, S, Rvol, dt, cA_new, cB_new, cC_new);
        update_energy_explicit (P, S, rs,   dt, Tf_new, Ts_new);
        apply_bc(P, cA_new, cB_new, cC_new, Tf_new, Ts_new);

        // commit
        S.cA.swap(cA_new); S.cB.swap(cB_new); S.cC.swap(cC_new);
        S.Tf.swap(Tf_new); S.Ts.swap(Ts_new);

        if (n % output_interval == 0)
            dump_snapshot(S);
    }

    return 0;
}