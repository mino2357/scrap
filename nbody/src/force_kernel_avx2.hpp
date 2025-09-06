#pragma once
#include <cstddef>

void compute_accel_dp_avx2(const double* x, const double* y, const double* z, const double* m,
                           double* ax, double* ay, double* az,
                           std::size_t N, double eps2, std::size_t Bj);

// Fused kick: v += tau * a(x)
void kick_accumulate_dp_avx2(const double* x, const double* y, const double* z, const double* m,
                             double* vx, double* vy, double* vz,
                             std::size_t N, double eps2, double tau, std::size_t Bj);
