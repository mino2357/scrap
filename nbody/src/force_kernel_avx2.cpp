#include "force_kernel_avx2.hpp"
#include <immintrin.h>
#include <algorithm>
#include <cmath>
#include <omp.h>

static inline __m256d rsqrt_dp_avx2(__m256d x) {
    __m128  xf   = _mm256_cvtpd_ps(x);
    __m128  y0f  = _mm_rsqrt_ps(xf);
    __m256d y    = _mm256_cvtps_pd(y0f);
    const __m256d half  = _mm256_set1_pd(0.5);
    const __m256d three = _mm256_set1_pd(1.5);
    __m256d yy = _mm256_mul_pd(y, y);
    __m256d t  = _mm256_fnmadd_pd(x, yy, three);
    y = _mm256_mul_pd(y, _mm256_mul_pd(half, t));
    yy = _mm256_mul_pd(y, y);
    t  = _mm256_fnmadd_pd(x, yy, three);
    y  = _mm256_mul_pd(y, _mm256_mul_pd(half, t));
    return y;
}

void compute_accel_dp_avx2(const double* x, const double* y, const double* z, const double* m,
                           double* ax, double* ay, double* az,
                           std::size_t N, double eps2, std::size_t Bj)
{
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t i=0;i<(ptrdiff_t)N;++i){ ax[i]=ay[i]=az[i]=0.0; }

    const __m256d veps2 = _mm256_set1_pd(eps2);
    ptrdiff_t Ivec = (ptrdiff_t)(N & ~std::size_t(3));

    #pragma omp parallel for schedule(static)
    for (ptrdiff_t ib=0; ib<Ivec; ib+=4){
        __m256d xi = _mm256_load_pd(x+ib);
        __m256d yi = _mm256_load_pd(y+ib);
        __m256d zi = _mm256_load_pd(z+ib);
        __m256d axi = _mm256_setzero_pd();
        __m256d ayi = _mm256_setzero_pd();
        __m256d azi = _mm256_setzero_pd();

        for (std::size_t jb=0; jb<N; jb+=Bj){
            std::size_t jn = std::min(jb+Bj, N);
            std::size_t j  = jb;
            for (; j+3<jn; j+=4){
                // j+0
                {
                    __m256d xj = _mm256_broadcast_sd(x + j + 0);
                    __m256d yj = _mm256_broadcast_sd(y + j + 0);
                    __m256d zj = _mm256_broadcast_sd(z + j + 0);
                    __m256d mj = _mm256_broadcast_sd(m + j + 0);
                    __m256d rx = _mm256_sub_pd(xj, xi);
                    __m256d ry = _mm256_sub_pd(yj, yi);
                    __m256d rz = _mm256_sub_pd(zj, zi);
                    __m256d r2 = _mm256_fmadd_pd(rx, rx, veps2);
                    r2         = _mm256_fmadd_pd(ry, ry, r2);
                    r2         = _mm256_fmadd_pd(rz, rz, r2);
                    __m256d inv = rsqrt_dp_avx2(r2);
                    __m256d inv3= _mm256_mul_pd(_mm256_mul_pd(inv,inv), inv);
                    __m256d a   = _mm256_mul_pd(mj, inv3);
                    axi = _mm256_fmadd_pd(a, rx, axi);
                    ayi = _mm256_fmadd_pd(a, ry, ayi);
                    azi = _mm256_fmadd_pd(a, rz, azi);
                }
                // j+1
                {
                    __m256d xj = _mm256_broadcast_sd(x + j + 1);
                    __m256d yj = _mm256_broadcast_sd(y + j + 1);
                    __m256d zj = _mm256_broadcast_sd(z + j + 1);
                    __m256d mj = _mm256_broadcast_sd(m + j + 1);
                    __m256d rx = _mm256_sub_pd(xj, xi);
                    __m256d ry = _mm256_sub_pd(yj, yi);
                    __m256d rz = _mm256_sub_pd(zj, zi);
                    __m256d r2 = _mm256_fmadd_pd(rx, rx, veps2);
                    r2         = _mm256_fmadd_pd(ry, ry, r2);
                    r2         = _mm256_fmadd_pd(rz, rz, r2);
                    __m256d inv = rsqrt_dp_avx2(r2);
                    __m256d inv3= _mm256_mul_pd(_mm256_mul_pd(inv,inv), inv);
                    __m256d a   = _mm256_mul_pd(mj, inv3);
                    axi = _mm256_fmadd_pd(a, rx, axi);
                    ayi = _mm256_fmadd_pd(a, ry, ayi);
                    azi = _mm256_fmadd_pd(a, rz, azi);
                }
                // j+2
                {
                    __m256d xj = _mm256_broadcast_sd(x + j + 2);
                    __m256d yj = _mm256_broadcast_sd(y + j + 2);
                    __m256d zj = _mm256_broadcast_sd(z + j + 2);
                    __m256d mj = _mm256_broadcast_sd(m + j + 2);
                    __m256d rx = _mm256_sub_pd(xj, xi);
                    __m256d ry = _mm256_sub_pd(yj, yi);
                    __m256d rz = _mm256_sub_pd(zj, zi);
                    __m256d r2 = _mm256_fmadd_pd(rx, rx, veps2);
                    r2         = _mm256_fmadd_pd(ry, ry, r2);
                    r2         = _mm256_fmadd_pd(rz, rz, r2);
                    __m256d inv = rsqrt_dp_avx2(r2);
                    __m256d inv3= _mm256_mul_pd(_mm256_mul_pd(inv,inv), inv);
                    __m256d a   = _mm256_mul_pd(mj, inv3);
                    axi = _mm256_fmadd_pd(a, rx, axi);
                    ayi = _mm256_fmadd_pd(a, ry, ayi);
                    azi = _mm256_fmadd_pd(a, rz, azi);
                }
                // j+3
                {
                    __m256d xj = _mm256_broadcast_sd(x + j + 3);
                    __m256d yj = _mm256_broadcast_sd(y + j + 3);
                    __m256d zj = _mm256_broadcast_sd(z + j + 3);
                    __m256d mj = _mm256_broadcast_sd(m + j + 3);
                    __m256d rx = _mm256_sub_pd(xj, xi);
                    __m256d ry = _mm256_sub_pd(yj, yi);
                    __m256d rz = _mm256_sub_pd(zj, zi);
                    __m256d r2 = _mm256_fmadd_pd(rx, rx, veps2);
                    r2         = _mm256_fmadd_pd(ry, ry, r2);
                    r2         = _mm256_fmadd_pd(rz, rz, r2);
                    __m256d inv = rsqrt_dp_avx2(r2);
                    __m256d inv3= _mm256_mul_pd(_mm256_mul_pd(inv,inv), inv);
                    __m256d a   = _mm256_mul_pd(mj, inv3);
                    axi = _mm256_fmadd_pd(a, rx, axi);
                    ayi = _mm256_fmadd_pd(a, ry, ayi);
                    azi = _mm256_fmadd_pd(a, rz, azi);
                }
            }
            for (; j<jn; ++j){
                __m256d xj = _mm256_broadcast_sd(x + j);
                __m256d yj = _mm256_broadcast_sd(y + j);
                __m256d zj = _mm256_broadcast_sd(z + j);
                __m256d mj = _mm256_broadcast_sd(m + j);
                __m256d rx = _mm256_sub_pd(xj, xi);
                __m256d ry = _mm256_sub_pd(yj, yi);
                __m256d rz = _mm256_sub_pd(zj, zi);
                __m256d r2 = _mm256_fmadd_pd(rx, rx, veps2);
                r2         = _mm256_fmadd_pd(ry, ry, r2);
                r2         = _mm256_fmadd_pd(rz, rz, r2);
                __m256d inv = rsqrt_dp_avx2(r2);
                __m256d inv3= _mm256_mul_pd(_mm256_mul_pd(inv,inv), inv);
                __m256d a   = _mm256_mul_pd(mj, inv3);
                axi = _mm256_fmadd_pd(a, rx, axi);
                ayi = _mm256_fmadd_pd(a, ry, ayi);
                azi = _mm256_fmadd_pd(a, rz, azi);
            }
        }
        _mm256_store_pd(ax+ib, axi);
        _mm256_store_pd(ay+ib, ayi);
        _mm256_store_pd(az+ib, azi);
    }

    for (std::size_t i=Ivec; i<N; ++i){
        double xi=x[i], yi=y[i], zi=z[i];
        double axi=0,ayi=0,azi=0;
        for (std::size_t jb=0; jb<N; jb+=Bj){
            std::size_t jn = std::min(jb+Bj, N);
            for (std::size_t j=jb; j<jn; ++j){
                double rx=x[j]-xi, ry=y[j]-yi, rz=z[j]-zi;
                double r2=rx*rx+ry*ry+rz*rz+eps2;
                double inv=1.0/std::sqrt(r2);
                double inv3=inv*inv*inv;
                double a=m[j]*inv3;
                axi += a*rx; ayi += a*ry; azi += a*rz;
            }
        }
        ax[i]=axi; ay[i]=ayi; az[i]=azi;
    }
}

void kick_accumulate_dp_avx2(const double* x, const double* y, const double* z, const double* m,
                             double* vx, double* vy, double* vz,
                             std::size_t N, double eps2, double tau, std::size_t Bj)
{
    const __m256d veps2 = _mm256_set1_pd(eps2);
    const __m256d vtau  = _mm256_set1_pd(tau);
    ptrdiff_t Ivec = (ptrdiff_t)(N & ~std::size_t(3));

    #pragma omp parallel for schedule(static)
    for (ptrdiff_t ib=0; ib<Ivec; ib+=4){
        __m256d xi = _mm256_load_pd(x+ib);
        __m256d yi = _mm256_load_pd(y+ib);
        __m256d zi = _mm256_load_pd(z+ib);
        __m256d vxi = _mm256_load_pd(vx+ib);
        __m256d vyi = _mm256_load_pd(vy+ib);
        __m256d vzi = _mm256_load_pd(vz+ib);

        for (std::size_t jb=0; jb<N; jb+=Bj){
            std::size_t jn = std::min(jb+Bj, N);
            std::size_t j  = jb;
            for (; j+3<jn; j+=4){
                // j+0..3
                for (int k=0;k<4;++k){
                    __m256d xj = _mm256_broadcast_sd(x + j + k);
                    __m256d yj = _mm256_broadcast_sd(y + j + k);
                    __m256d zj = _mm256_broadcast_sd(z + j + k);
                    __m256d mj = _mm256_broadcast_sd(m + j + k);
                    __m256d rx = _mm256_sub_pd(xj, xi);
                    __m256d ry = _mm256_sub_pd(yj, yi);
                    __m256d rz = _mm256_sub_pd(zj, zi);
                    __m256d r2 = _mm256_fmadd_pd(rx, rx, veps2);
                    r2         = _mm256_fmadd_pd(ry, ry, r2);
                    r2         = _mm256_fmadd_pd(rz, rz, r2);
                    __m256d inv = rsqrt_dp_avx2(r2);
                    __m256d inv3= _mm256_mul_pd(_mm256_mul_pd(inv,inv), inv);
                    __m256d a   = _mm256_mul_pd(mj, inv3);
                    __m256d tau_a = _mm256_mul_pd(vtau, a);
                    vxi = _mm256_fmadd_pd(tau_a, rx, vxi);
                    vyi = _mm256_fmadd_pd(tau_a, ry, vyi);
                    vzi = _mm256_fmadd_pd(tau_a, rz, vzi);
                }
            }
            for (; j<jn; ++j){
                __m256d xj = _mm256_broadcast_sd(x + j);
                __m256d yj = _mm256_broadcast_sd(y + j);
                __m256d zj = _mm256_broadcast_sd(z + j);
                __m256d mj = _mm256_broadcast_sd(m + j);
                __m256d rx = _mm256_sub_pd(xj, xi);
                __m256d ry = _mm256_sub_pd(yj, yi);
                __m256d rz = _mm256_sub_pd(zj, zi);
                __m256d r2 = _mm256_fmadd_pd(rx, rx, veps2);
                r2         = _mm256_fmadd_pd(ry, ry, r2);
                r2         = _mm256_fmadd_pd(rz, rz, r2);
                __m256d inv = rsqrt_dp_avx2(r2);
                __m256d inv3= _mm256_mul_pd(_mm256_mul_pd(inv,inv), inv);
                __m256d a   = _mm256_mul_pd(mj, inv3);
                __m256d tau_a = _mm256_mul_pd(vtau, a);
                vxi = _mm256_fmadd_pd(tau_a, rx, vxi);
                vyi = _mm256_fmadd_pd(tau_a, ry, vyi);
                vzi = _mm256_fmadd_pd(tau_a, rz, vzi);
            }
        }
        _mm256_store_pd(vx+ib, vxi);
        _mm256_store_pd(vy+ib, vyi);
        _mm256_store_pd(vz+ib, vzi);
    }

    for (std::size_t i=Ivec; i<N; ++i){
        double xi_=x[i], yi_=y[i], zi_=z[i];
        double vxi=vx[i], vyi=vy[i], vzi=vz[i];
        for (std::size_t jb=0; jb<N; jb+=Bj){
            std::size_t jn = std::min(jb+Bj, N);
            for (std::size_t j=jb; j<jn; ++j){
                double rx=x[j]-xi_, ry=y[j]-yi_, rz=z[j]-zi_;
                double r2=rx*rx+ry*ry+rz*rz+eps2;
                double inv=1.0/std::sqrt(r2);
                double inv3=inv*inv*inv;
                double a=m[j]*inv3, tau_a=tau*a;
                vxi += tau_a*rx; vyi += tau_a*ry; vzi += tau_a*rz;
            }
        }
        vx[i]=vxi; vy[i]=vyi; vz[i]=vzi;
    }
}
