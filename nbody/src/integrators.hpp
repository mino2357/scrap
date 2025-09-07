#pragma once
#include <cstddef>

/*
  時間積分インタフェース
  ----------------------------------------
  - Arrays: SoA 形式で N 粒子の状態を保持。
    - 位置: x[], y[], z[]
    - 速度: vx[], vy[], vz[]
    - 質量: m[]
    - 加速度ワーク: ax[], ay[], az[]（two-pass で使用）
  - メモリ整列: 可能なら 32B/64B 整列（AVX2 の整列ロード/ストアに有利）。
  - eps2: ソフトニング ε^2。
  - Bj: 力計算カーネルの j タイル幅（ベンチにより調整）。
*/
struct Arrays {
  std::size_t N;
  double *x,*y,*z;
  double *vx,*vy,*vz;
  double *m;
  double *ax,*ay,*az;
};

// 4 段 Yoshida（シンプレクティック, 4 次）
void step_yoshida4(double dt, Arrays& S, double eps2, std::size_t Bj);
// 古典 4 次 Runge-Kutta（非シンプレクティック）
void step_rk4(double dt, Arrays& S, double eps2, std::size_t Bj);
// Velocity-Verlet（シンプレクティック, 2 次）
void step_verlet(double dt, Arrays& S, double eps2, std::size_t Bj);
