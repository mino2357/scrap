//
// integrators.hpp
// -----------------
// Lightweight explicit Runge–Kutta integrators used by the case04 driver.
// The functions here advance the state vector `y` from time `t0` to `t1`
// under constant pressure `P`.  Adaptive schemes accept relative and
// absolute tolerances and internally use an L2-norm based error control.
//

#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include "chem.hpp"

// Generic function pointer type for a time integrator.  The method updates
// `y` in-place while keeping track of the provided tolerances.
template<typename T>
using Integrator = void(*)(std::vector<T>&, T, T, T,
                          T, T,
                          const std::vector<Reaction<T>>&, const std::vector<ThermoData<T>>&);

// Classical fourth-order Runge–Kutta with a fixed number of steps.
// The tolerance arguments are ignored; this routine is mainly used for
// reference and debugging.
template<typename T>
inline void rk4(std::vector<T>& y, T t0, T t1, T P,
               T /*rtol*/, T /*atol*/,
               const std::vector<Reaction<T>>& reactions,
               const std::vector<ThermoData<T>>& thermo){
    // Use a simple constant step size based on 1000 sub-intervals.
    T h = (t1 - t0) / T(1000);
    T t = t0;
    size_t m = y.size();
    std::vector<T> k1(m),k2(m),k3(m),k4(m),yt(m);
    while(t < t1){
        if(t + h > t1) h = t1 - t;        // Last step hits t1 exactly
        compute_rhs(reactions, thermo, P, y, k1);
        for(size_t i=0;i<m;++i) yt[i] = y[i] + T(0.5)*h*k1[i];
        compute_rhs(reactions, thermo, P, yt, k2);
        for(size_t i=0;i<m;++i) yt[i] = y[i] + T(0.5)*h*k2[i];
        compute_rhs(reactions, thermo, P, yt, k3);
        for(size_t i=0;i<m;++i) yt[i] = y[i] + h*k3[i];
        compute_rhs(reactions, thermo, P, yt, k4);
        for(size_t i=0;i<m;++i)
            y[i] += (h/T(6))*(k1[i] + T(2)*k2[i] + T(2)*k3[i] + k4[i]);
        for(size_t i=0;i<m-1;++i)
            if(y[i] < T(0)) y[i] = T(0);
        t += h;
    }
}

// Fifth-order Dormand–Prince method with adaptive step-size control.
// The fourth- and fifth-order solutions are compared to estimate the
// local truncation error.  A scaled L2 norm is used to decide if the
// step is accepted and how to adapt the step size.
template<typename T>
inline void rk45(std::vector<T>& y, T t0, T t1, T P,
                T rtol, T atol,
                const std::vector<Reaction<T>>& reactions,
                const std::vector<ThermoData<T>>& thermo){
    const T safety = T(0.9);              // Conservative step-size factor
    T h = (t1 - t0) / T(1000);            // Initial guess for step size
    T t = t0;
    const size_t m = y.size();
    std::vector<T> k1(m),k2(m),k3(m),k4(m),k5(m),k6(m),yt(m),y4(m),y5(m),errv(m);
    while(t < t1){
        if(t + h > t1) h = t1 - t;
        T err;
        do{
            // Evaluate intermediate stages
            compute_rhs(reactions, thermo, P, y, k1);
            for(size_t i=0;i<m;++i) yt[i] = y[i] + h*(T(1.0)/T(4.0))*k1[i];
            compute_rhs(reactions, thermo, P, yt, k2);
            for(size_t i=0;i<m;++i) yt[i] = y[i] + h*(T(3.0)/T(32.0)*k1[i] + T(9.0)/T(32.0)*k2[i]);
            compute_rhs(reactions, thermo, P, yt, k3);
            for(size_t i=0;i<m;++i)
                yt[i] = y[i] + h*(T(1932.0)/T(2197.0)*k1[i] - T(7200.0)/T(2197.0)*k2[i] + T(7296.0)/T(2197.0)*k3[i]);
            compute_rhs(reactions, thermo, P, yt, k4);
            for(size_t i=0;i<m;++i)
                yt[i] = y[i] + h*(T(439.0)/T(216.0)*k1[i] - T(8.0)*k2[i] + T(3680.0)/T(513.0)*k3[i] - T(845.0)/T(4104.0)*k4[i]);
            compute_rhs(reactions, thermo, P, yt, k5);
            for(size_t i=0;i<m;++i)
                yt[i] = y[i] + h*( -T(8.0)/T(27.0)*k1[i] + T(2.0)*k2[i] - T(3544.0)/T(2565.0)*k3[i]
                                   + T(1859.0)/T(4104.0)*k4[i] - T(11.0)/T(40.0)*k5[i]);
            compute_rhs(reactions, thermo, P, yt, k6);

            // Form both orders of the solution and their difference
            for(size_t i=0;i<m;++i){
                y4[i] = y[i] + h*( T(25.0)/T(216.0)*k1[i] + T(1408.0)/T(2565.0)*k3[i]
                                   + T(2197.0)/T(4104.0)*k4[i] - T(1.0)/T(5.0)*k5[i] );
                y5[i] = y[i] + h*( T(16.0)/T(135.0)*k1[i] + T(6656.0)/T(12825.0)*k3[i]
                                   + T(28561.0)/T(56430.0)*k4[i] - T(9.0)/T(50.0)*k5[i]
                                   + T(2.0)/T(55.0)*k6[i] );
                errv[i] = y5[i] - y4[i];
            }

            // Compute the scaled L2 error norm
            err = T(0);
            for(size_t i=0;i<m;++i){
                T sc = atol + rtol * std::max(std::abs(y[i]), std::abs(y5[i]));
                T ratio = (sc>0)? errv[i] / sc : errv[i];
                err += ratio * ratio;
            }
            err = std::sqrt(err / T(m));
            if(err> T(1)){
                // Reject step and reduce the step size
                T fac = (err>0)? safety*std::pow(T(1)/err, T(0.2)) : T(0.5);
                fac = std::min(T(5.0), std::max(T(0.1), fac));
                h *= fac;
                if(t + h > t1) h = t1 - t;
            }
        }while(err> T(1));

        // Accept the step and propose a new step size
        y = y5;
        for(size_t i=0;i<m-1;++i)
            if(y[i] < T(0)) y[i] = T(0);
        t += h;

        T fac = (err>0)? safety*std::pow(T(1)/err, T(0.2)) : T(5.0);
        fac = std::min(T(5.0), std::max(T(0.1), fac));
        h *= fac;
    }
}

// High-order Dormand–Prince 7/8 method taken from Boost.Odeint.
// Similar to rk45 but using a higher-order pair and 13 stages.
template<typename T>
inline void rk78(std::vector<T>& y, T t0, T t1, T P,
                T rtol, T atol,
                const std::vector<Reaction<T>>& reactions,
                const std::vector<ThermoData<T>>& thermo){
    const T safety = T(0.9);
    T h = (t1 - t0) / T(1000);
    T t = t0;
    const size_t m = y.size();
    std::vector<std::vector<T>> k(13, std::vector<T>(m));
    std::vector<T> yt(m), y8(m), errv(m);
    // Coefficients from Boost.Odeint
    const T a[13][13] = {
        {0},
        { T(2.0)/T(27.0) },
        { T(1.0)/T(36.0), T(1.0)/T(12.0) },
        { T(1.0)/T(24.0), T(0), T(1.0)/T(8.0) },
        { T(5.0)/T(12.0), T(0), T(-25.0)/T(16.0), T(25.0)/T(16.0) },
        { T(1.0)/T(20.0), T(0), T(0), T(1.0)/T(4.0), T(1.0)/T(5.0) },
        { T(-25.0)/T(108.0), T(0), T(0), T(125.0)/T(108.0), T(-65.0)/T(27.0), T(125.0)/T(54.0) },
        { T(31.0)/T(300.0), T(0), T(0), T(0), T(61.0)/T(225.0), T(-2.0)/T(9.0), T(13.0)/T(900.0) },
        { T(2.0), T(0), T(0), T(-53.0)/T(6.0), T(704.0)/T(45.0), T(-107.0)/T(9.0), T(67.0)/T(90.0), T(3.0) },
        { T(-91.0)/T(108.0), T(0), T(0), T(23.0)/T(108.0), T(-976.0)/T(135.0), T(311.0)/T(54.0), T(-19.0)/T(60.0), T(17.0)/T(6.0), T(-1.0)/T(12.0) },
        { T(2383.0)/T(4100.0), T(0), T(0), T(-341.0)/T(164.0), T(4496.0)/T(1025.0), T(-301.0)/T(82.0), T(2133.0)/T(4100.0), T(45.0)/T(82.0), T(45.0)/T(164.0), T(18.0)/T(41.0) },
        { T(3.0)/T(205.0), T(0), T(0), T(0), T(0), T(-6.0)/T(41.0), T(-3.0)/T(205.0), T(-3.0)/T(41.0), T(3.0)/T(41.0), T(6.0)/T(41.0), T(0) },
        { T(-1777.0)/T(4100.0), T(0), T(0), T(-341.0)/T(164.0), T(4496.0)/T(1025.0), T(-289.0)/T(82.0), T(2193.0)/T(4100.0), T(51.0)/T(82.0), T(33.0)/T(164.0), T(12.0)/T(41.0), T(0), T(1.0) }
    };
    const T b[13] = { T(0), T(0), T(0), T(0), T(0), T(34.0)/T(105.0), T(9.0)/T(35.0), T(9.0)/T(35.0),
                      T(9.0)/T(280.0), T(9.0)/T(280.0), T(0), T(41.0)/T(840.0), T(41.0)/T(840.0) };
    const T db[13] = { -T(41.0)/T(840.0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), -T(41.0)/T(840.0), T(41.0)/T(840.0), T(41.0)/T(840.0) };

    while(t < t1){
        if(t + h > t1) h = t1 - t;
        T err;
        do{
            // Build all intermediate stages
            compute_rhs(reactions, thermo, P, y, k[0]);
            for(int s=1; s<13; ++s){
                for(size_t i=0;i<m;++i){
                    yt[i] = y[i];
                    for(int j=0;j<s; ++j) yt[i] += h * a[s][j] * k[j][i];
                }
                compute_rhs(reactions, thermo, P, yt, k[s]);
            }
            // Combine stages for 8th order solution and embedded error estimate
            for(size_t i=0;i<m;++i){
                T sum_b = T(0), sum_e = T(0);
                for(int s=0; s<13; ++s){
                    sum_b += b[s]*k[s][i];
                    sum_e += db[s]*k[s][i];
                }
                y8[i] = y[i] + h*sum_b;
                errv[i] = h*sum_e;
            }
            // Scaled L2 error norm
            err = T(0);
            for(size_t i=0;i<m;++i){
                T sc = atol + rtol * std::max(std::abs(y[i]), std::abs(y8[i]));
                T ratio = (sc>0)? errv[i] / sc : errv[i];
                err += ratio * ratio;
            }
            err = std::sqrt(err / T(m));
            if(err> T(1)){
                T fac = (err>0)? safety*std::pow(T(1)/err, T(1.0)/T(8.0)) : T(0.5);
                fac = std::min(T(4.0), std::max(T(0.1), fac));
                h *= fac;
                if(t + h > t1) h = t1 - t;
            }
        }while(err> T(1));

        // Accept the step and adjust step size for the next iteration
        y = y8;
        for(size_t i=0;i<m-1;++i)
            if(y[i] < T(0)) y[i] = T(0);
        t += h;

        T fac = (err>0)? safety*std::pow(T(1)/err, T(1.0)/T(8.0)) : T(4.0);
        fac = std::min(T(4.0), std::max(T(0.1), fac));
        h *= fac;
    }
}

