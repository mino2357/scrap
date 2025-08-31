#include <cmath>
#include <vector>
#include <cassert>

#define compute_rhs compute_rhs_chem
#include "chem.hpp"
#undef compute_rhs
#define RK_CLAMP_NEGATIVE 0
#include "integrators.hpp"

// RHS for simple harmonic oscillator x'' = -x
// y[0] = x, y[1] = v

template <typename T>
void compute_rhs(const std::vector<Reaction<T>>&, const std::vector<ThermoData<T>>&, T /*P*/, const std::vector<T>& y, std::vector<T>& dy) {
    dy[0] = y[1];
    dy[1] = -y[0];
}

template <typename Integrator>
void run_tests(Integrator integrator) {
    using T = double;
    std::vector<Reaction<T>> reactions; // unused
    std::vector<ThermoData<T>> thermo;  // unused
    const T P = 0;
    const T rtol = 1e-9;
    const T atol = 1e-12;
    const T dt = 1e-3;

    auto check = [&](T t_end, T x_expect, T v_expect) {
        std::vector<T> y = {T(0), T(1)};
        integrator(y, T(0), t_end, dt, P, rtol, atol, reactions, thermo);
        assert(std::abs(y[0] - x_expect) < 1e-5);
        assert(std::abs(y[1] - v_expect) < 1e-5);
    };

    check(M_PI/2, 1.0, 0.0);
    check(M_PI,   0.0, -1.0);
    check(3*M_PI/2, -1.0, 0.0);
    check(2*M_PI, 0.0, 1.0);
}

int main(){
    run_tests(rk4<double>);
    run_tests(rk45<double>);
    run_tests(rk78<double>);
    return 0;
}
