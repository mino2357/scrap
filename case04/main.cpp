//
// main.cpp
// ---------
// Simple driver for case04.  Loads a small reaction mechanism and integrates
// it using one of the Runge–Kutta methods defined in integrators.hpp.  Results
// are periodically written to `case04.dat` for plotting.

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <sstream>
#include <limits>

#include "chem.hpp"
#include "integrators.hpp"

int main(int argc, char** argv){
    using T = double;

    // Allow overriding mechanism files from the command line ----------------
    std::string mech_file  = "chem.inp";
    std::string thermo_file = "therm.dat";
    if(argc > 2) mech_file  = argv[2];
    if(argc > 3) thermo_file = argv[3];

    // Parse mechanism files -------------------------------------------------
    std::vector<std::string> species;
    std::map<std::string,int> idx;
    std::vector<Reaction<T>> reactions;
    if(!load_chemkin<T>(mech_file, species, idx, reactions)){
        std::cerr << mech_file << " not found\n";
        return 1;
    }
    std::vector<ThermoData<T>> thermo;
    if(!load_thermo<T>(thermo_file, idx, thermo)){
        std::cerr << thermo_file << " not found\n";
        return 1;
    }

    // Default initial composition: H2 1, O2 1, N2 3.76 ----------------------
    std::vector<T> X(species.size(), T(1e-8));
    auto set_init = [&](const std::string& name, T val){
        auto it = idx.find(name);
        if(it!=idx.end()) X[it->second] = val;
    };
    set_init("H2", T(1.0));
    set_init("O2", T(1.0));
    set_init("N2", T(3.76));

    // Simulation parameters with defaults matching the README ----------------
    T P = T(202650.0); // Pa
    T T0 = T(1000.0);
    T t_end = T(1e-3);
    T output_interval = T(1e-5);
    T dt = output_interval / T(1000);   // Initial step size
    T rtol = T(1e-6);                   // Relative tolerance
    T atol = T(1e-12);                  // Absolute tolerance

    // Optional overrides from input.start -----------------------------------
    std::ifstream cfg("input.start");
    if(cfg){
        std::string key;
        while(cfg >> key){
            if(key[0]=='#'){ cfg.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); continue; }
            if(key=="END") break;
            if(key=="PRES")      cfg >> P;
            else if(key=="TEMP") cfg >> T0;
            else if(key=="TIME") cfg >> t_end;
            else if(key=="DELT") { cfg >> output_interval; dt = output_interval / T(1000); }
            else if(key=="DT")   cfg >> dt;
            else if(key=="RTOL") cfg >> rtol;
            else if(key=="ATOL") cfg >> atol;
            else if(key=="ENERG" || key=="ENERGY"){
                cfg.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            } else{
                T val; cfg >> val; set_init(key, val);
            }
        }
    }

    // Normalize mole fractions
    T X_sum = T(0);
    for(T v : X) X_sum += v;
    if(X_sum > T(0)) for(T& v : X) v /= X_sum;

    std::vector<T> y = X;
    y.push_back(T0); // Temperature appended as last component

    // Select integrator based on command line argument ----------------------
    Integrator<T> integrator = rk4<T>;
    if(argc > 1){
        std::string method = argv[1];
        if(method == "rk45") integrator = rk45<T>;
        else if(method == "rk78") integrator = rk78<T>;
    }

    // March solution in time ------------------------------------------------
    T t = T(0);
    std::vector<T> times{t};
    std::vector<std::vector<T>> history{y};
    while(t < t_end){
        T next = std::min(t + output_interval, t_end);
        integrator(y, t, next, dt, P, rtol, atol, reactions, thermo);
        t = next;
        times.push_back(t);
        history.push_back(y);
    }

    // Write history to file for gnuplot ------------------------------------
    const size_t n = species.size();
    const size_t m = n + 1;
    std::vector<T> dy(m);
    std::ofstream ofs("case04.dat");
    std::ofstream ofs_conc("case04_conc.dat");
    ofs.setf(std::ios::scientific);
    ofs.precision(15);
    ofs_conc.setf(std::ios::scientific);
    ofs_conc.precision(15);
    ofs << "time";
    ofs_conc << "time";
    for(const auto& s : species){
        ofs << ' ' << s;
        ofs_conc << ' ' << s;
    }
    ofs << " T dTdt\n";
    ofs_conc << '\n';
    const T R = T(8.3144621);
    for(size_t k=0;k<times.size();++k){
        T sum=T(0); for(size_t i=0;i<n;++i) sum+=history[k][i];
        compute_rhs(reactions, thermo, P, history[k], dy);
        ofs << times[k];
        ofs_conc << times[k];
        T Ctot = P/(R*history[k][n]);
        for(size_t i=0;i<n;++i){
            T Xi = (sum>0? history[k][i]/sum : T(0));
            ofs << ' ' << Xi;
            ofs_conc << ' ' << Xi*Ctot;
        }
        ofs << ' ' << history[k][n] << ' ' << dy[n] << '\n';
        ofs_conc << '\n';
    }

    std::cout.setf(std::ios::scientific);
    std::cout.precision(15);
    for(size_t i=0;i<n;++i)
        std::cout << species[i] << " " << y[i] << "\n";
    return 0;
}

