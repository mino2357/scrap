#pragma once
#include <cstddef>

struct Arrays {
  std::size_t N;
  double *x,*y,*z;
  double *vx,*vy,*vz;
  double *m;
  double *ax,*ay,*az;
};

void step_yoshida4(double dt, Arrays& S, double eps2, std::size_t Bj);
void step_rk4(double dt, Arrays& S, double eps2, std::size_t Bj);
void step_verlet(double dt, Arrays& S, double eps2, std::size_t Bj);
