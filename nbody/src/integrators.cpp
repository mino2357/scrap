#include "integrators.hpp"
#include "force_kernel_avx2.hpp"
#include <algorithm>
#include <vector>
#include <immintrin.h>
#include <omp.h>

static bool g_kernel_fused = true;
extern "C" bool* __get_kernel_flag(){ return &g_kernel_fused; } // optional

static inline void drift(double* x,double* y,double* z,
                         const double* vx,const double* vy,const double* vz,
                         std::size_t N, double tau){
  #pragma omp parallel for schedule(static)
  for (ptrdiff_t i=0;i<(ptrdiff_t)N;++i){
    x[i] += tau * vx[i];
    y[i] += tau * vy[i];
    z[i] += tau * vz[i];
  }
}

static inline void kick(double* vx,double* vy,double* vz,
                        double* ax,double* ay,double* az,
                        const double* x,const double* y,const double* z,const double* m,
                        std::size_t N, double tau, double eps2, std::size_t Bj){
  if (g_kernel_fused){
    (void)ax; (void)ay; (void)az;
    kick_accumulate_dp_avx2(x,y,z,m, vx,vy,vz, N, eps2, tau, Bj);
  } else {
    compute_accel_dp_avx2(x,y,z,m, ax,ay,az, N, eps2, Bj);
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t i=0;i<(ptrdiff_t)N;++i){
      vx[i] += tau * ax[i];
      vy[i] += tau * ay[i];
      vz[i] += tau * az[i];
    }
  }
}

// Yoshida4 coefficients
static constexpr double c1 =  0.6756035959798289;
static constexpr double c2 = -0.17560359597982877;
static constexpr double c3 =  c2;
static constexpr double c4 =  c1;
static constexpr double d1 =  1.3512071919596578;
static constexpr double d2 = -1.7024143839193153;
static constexpr double d3 =  d1;

void step_yoshida4(double dt, Arrays& S, double eps2, std::size_t Bj){
  kick(S.vx,S.vy,S.vz, S.ax,S.ay,S.az, S.x,S.y,S.z,S.m, S.N, d1*dt, eps2, Bj);
  drift(S.x,S.y,S.z, S.vx,S.vy,S.vz, S.N, c1*dt);
  kick(S.vx,S.vy,S.vz, S.ax,S.ay,S.az, S.x,S.y,S.z,S.m, S.N, d2*dt, eps2, Bj);
  drift(S.x,S.y,S.z, S.vx,S.vy,S.vz, S.N, c2*dt);
  kick(S.vx,S.vy,S.vz, S.ax,S.ay,S.az, S.x,S.y,S.z,S.m, S.N, d3*dt, eps2, Bj);
  drift(S.x,S.y,S.z, S.vx,S.vy,S.vz, S.N, (c3+c4)*dt);
}

void step_rk4(double dt, Arrays& S, double eps2, std::size_t Bj){
  const std::size_t N = S.N;
  static thread_local std::vector<double> k1x, k1y, k1z, k1vx, k1vy, k1vz;
  static thread_local std::vector<double> k2x, k2y, k2z, k2vx, k2vy, k2vz;
  static thread_local std::vector<double> k3x, k3y, k3z, k3vx, k3vy, k3vz;
  static thread_local std::vector<double> k4x, k4y, k4z, k4vx, k4vy, k4vz;

  auto ensure = [&](auto& v){ if (v.size()!=N) v.assign(N,0.0); };
  ensure(k1x); ensure(k1y); ensure(k1z); ensure(k1vx); ensure(k1vy); ensure(k1vz);
  ensure(k2x); ensure(k2y); ensure(k2z); ensure(k2vx); ensure(k2vy); ensure(k2vz);
  ensure(k3x); ensure(k3y); ensure(k3z); ensure(k3vx); ensure(k3vy); ensure(k3vz);
  ensure(k4x); ensure(k4y); ensure(k4z); ensure(k4vx); ensure(k4vy); ensure(k4vz);

  compute_accel_dp_avx2(S.x,S.y,S.z,S.m, S.ax,S.ay,S.az, N, eps2, Bj);
  #pragma omp parallel for schedule(static)
  for (ptrdiff_t i=0;i<(ptrdiff_t)N;++i){
    k1x[i]=S.vx[i]; k1y[i]=S.vy[i]; k1z[i]=S.vz[i];
    k1vx[i]=S.ax[i]; k1vy[i]=S.ay[i]; k1vz[i]=S.az[i];
  }

  std::vector<double> x2(N),y2(N),z2(N),vx2(N),vy2(N),vz2(N);
  #pragma omp parallel for schedule(static)
  for (ptrdiff_t i=0;i<(ptrdiff_t)N;++i){
    x2[i]=S.x[i]+0.5*dt*k1x[i];
    y2[i]=S.y[i]+0.5*dt*k1y[i];
    z2[i]=S.z[i]+0.5*dt*k1z[i];
    vx2[i]=S.vx[i]+0.5*dt*k1vx[i];
    vy2[i]=S.vy[i]+0.5*dt*k1vy[i];
    vz2[i]=S.vz[i]+0.5*dt*k1vz[i];
  }
  compute_accel_dp_avx2(x2.data(),y2.data(),z2.data(),S.m, S.ax,S.ay,S.az, N, eps2, Bj);
  #pragma omp parallel for schedule(static)
  for (ptrdiff_t i=0;i<(ptrdiff_t)N;++i){
    k2x[i]=vx2[i]; k2y[i]=vy2[i]; k2z[i]=vz2[i];
    k2vx[i]=S.ax[i]; k2vy[i]=S.ay[i]; k2vz[i]=S.az[i];
  }

  std::vector<double> x3(N),y3(N),z3(N),vx3(N),vy3(N),vz3(N);
  #pragma omp parallel for schedule(static)
  for (ptrdiff_t i=0;i<(ptrdiff_t)N;++i){
    x3[i]=S.x[i]+0.5*dt*k2x[i];
    y3[i]=S.y[i]+0.5*dt*k2y[i];
    z3[i]=S.z[i]+0.5*dt*k2z[i];
    vx3[i]=S.vx[i]+0.5*dt*k2vx[i];
    vy3[i]=S.vy[i]+0.5*dt*k2vy[i];
    vz3[i]=S.vz[i]+0.5*dt*k2vz[i];
  }
  compute_accel_dp_avx2(x3.data(),y3.data(),z3.data(),S.m, S.ax,S.ay,S.az, N, eps2, Bj);
  #pragma omp parallel for schedule(static)
  for (ptrdiff_t i=0;i<(ptrdiff_t)N;++i){
    k3x[i]=vx3[i]; k3y[i]=vy3[i]; k3z[i]=vz3[i];
    k3vx[i]=S.ax[i]; k3vy[i]=S.ay[i]; k3vz[i]=S.az[i];
  }

  std::vector<double> x4(N),y4(N),z4(N),vx4(N),vy4(N),vz4(N);
  #pragma omp parallel for schedule(static)
  for (ptrdiff_t i=0;i<(ptrdiff_t)N;++i){
    x4[i]=S.x[i]+dt*k3x[i];
    y4[i]=S.y[i]+dt*k3y[i];
    z4[i]=S.z[i]+dt*k3z[i];
    vx4[i]=S.vx[i]+dt*k3vx[i];
    vy4[i]=S.vy[i]+dt*k3vy[i];
    vz4[i]=S.vz[i]+dt*k3vz[i];
  }
  compute_accel_dp_avx2(x4.data(),y4.data(),z4.data(),S.m, S.ax,S.ay,S.az, N, eps2, Bj);
  #pragma omp parallel for schedule(static)
  for (ptrdiff_t i=0;i<(ptrdiff_t)N;++i){
    k4x[i]=vx4[i]; k4y[i]=vy4[i]; k4z[i]=vz4[i];
    k4vx[i]=S.ax[i]; k4vy[i]=S.ay[i]; k4vz[i]=S.az[i];
  }

  #pragma omp parallel for schedule(static)
  for (ptrdiff_t i=0;i<(ptrdiff_t)N;++i){
    S.x[i]  += dt*(k1x[i] + 2*k2x[i] + 2*k3x[i] + k4x[i])*(1.0/6.0);
    S.y[i]  += dt*(k1y[i] + 2*k2y[i] + 2*k3y[i] + k4y[i])*(1.0/6.0);
    S.z[i]  += dt*(k1z[i] + 2*k2z[i] + 2*k3z[i] + k4z[i])*(1.0/6.0);
    S.vx[i] += dt*(k1vx[i]+ 2*k2vx[i]+ 2*k3vx[i]+ k4vx[i])*(1.0/6.0);
    S.vy[i] += dt*(k1vy[i]+ 2*k2vy[i]+ 2*k3vy[i]+ k4vy[i])*(1.0/6.0);
    S.vz[i] += dt*(k1vz[i]+ 2*k2vz[i]+ 2*k3vz[i]+ k4vz[i])*(1.0/6.0);
  }
}

void step_verlet(double dt, Arrays& S, double eps2, std::size_t Bj){
  kick(S.vx,S.vy,S.vz, S.ax,S.ay,S.az, S.x,S.y,S.z,S.m, S.N, 0.5*dt, eps2, Bj);
  drift(S.x,S.y,S.z, S.vx,S.vy,S.vz, S.N, dt);
  kick(S.vx,S.vy,S.vz, S.ax,S.ay,S.az, S.x,S.y,S.z,S.m, S.N, 0.5*dt, eps2, Bj);
}
