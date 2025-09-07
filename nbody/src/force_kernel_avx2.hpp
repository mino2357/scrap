#pragma once
#include <cstddef>

/*
  AVX2 ベクトル化された力計算カーネル（インタフェース）
  ----------------------------------------
  前提
  - 配列は SoA 形式（x[],y[],z[],m[] など）。
  - 可能なら 32B/64B 整列（_mm256_load_pd/_store_pd を使用）。
  - N は任意。内部で 4 要素ずつ処理し、端数はスカラで処理します。
  - eps2 はソフトニング項 ε^2。
  - Bj は j タイル幅（キャッシュ局所性向上）。

  関数
  - compute_accel_dp_avx2:
      a_i = Σ_j m_j inv^3 (r_j - r_i) を計算して ax, ay, az に格納。
  - kick_accumulate_dp_avx2:
      v_i += tau * a_i(x) を融合実行（a を明示配列に持たない）。
*/

void compute_accel_dp_avx2(const double* x, const double* y, const double* z, const double* m,
                           double* ax, double* ay, double* az,
                           std::size_t N, double eps2, std::size_t Bj);

// Fused kick: v += tau * a(x)
void kick_accumulate_dp_avx2(const double* x, const double* y, const double* z, const double* m,
                             double* vx, double* vy, double* vz,
                             std::size_t N, double eps2, double tau, std::size_t Bj);
