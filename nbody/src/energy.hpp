#pragma once
#include <cstddef>
#include <cmath>
#include <omp.h>

/*
  エネルギーと角運動量の計算ユーティリティ
  ----------------------------------------
  - 目的: 運動エネルギー T、ポテンシャルエネルギー V、角運動量 L を
    スカラ実装で計算します（OpenMP による単純並列化）。
  - データは SoA 配置（x[],y[],z[], vx[],vy[],vz[], m[]）。
  - いずれの関数も入力配列の長さは N であることが前提です。
  - 並列化: 各 for ループを OpenMP の reduction で合算。
*/

// 運動エネルギー T = Σ_i 1/2 m_i |v_i|^2
inline double kinetic_energy(const double* vx, const double* vy, const double* vz,
                             const double* m, std::size_t N){
  double T = 0.0;
  #pragma omp parallel for reduction(+:T) schedule(static)
  for (ptrdiff_t i=0;i<(ptrdiff_t)N;++i){
    double v2 = vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    T += 0.5 * m[i] * v2;
  }
  return T;
}

// ポテンシャルエネルギー V = Σ_{i<j} - m_i m_j / sqrt(|r_j - r_i|^2 + eps2)
// - 二重和の両側を数えないため、i<j のみを走査。
// - O(N^2) の計算で、チェック用/統計用を想定。
inline double potential_energy(const double* x, const double* y, const double* z,
                               const double* m, std::size_t N, double eps2){
  double V = 0.0;
  #pragma omp parallel for reduction(+:V) schedule(static)
  for (ptrdiff_t i=0;i<(ptrdiff_t)N;++i){
    double Vi = 0.0;
    for (std::size_t j=i+1;j<N;++j){
      double dx = x[j]-x[i];
      double dy = y[j]-y[i];
      double dz = z[j]-z[i];
      double r2 = dx*dx + dy*dy + dz*dz + eps2;
      double invr = 1.0/std::sqrt(r2);
      Vi += - m[i]*m[j]*invr;
    }
    V += Vi;
  }
  return V;
}

// 角運動量 L = Σ_i r_i × (m_i v_i)
// - L は 3 要素配列 [Lx, Ly, Lz] に出力。
inline void angular_momentum(const double* x,const double* y,const double* z,
                             const double* vx,const double* vy,const double* vz,
                             const double* m, std::size_t N, double L[3]){
  double Lx=0, Ly=0, Lz=0;
  #pragma omp parallel for reduction(+:Lx,Ly,Lz) schedule(static)
  for (ptrdiff_t i=0;i<(ptrdiff_t)N;++i){
    double px = m[i]*vx[i], py = m[i]*vy[i], pz = m[i]*vz[i];
    Lx += y[i]*pz - z[i]*py;
    Ly += z[i]*px - x[i]*pz;
    Lz += x[i]*py - y[i]*px;
  }
  L[0]=Lx; L[1]=Ly; L[2]=Lz;
}
