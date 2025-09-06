#pragma once
#include <cstddef>
#include <cmath>
#include <omp.h>

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
