//
// chem.hpp
// ---------
// Minimal utilities for parsing CHEMKIN style mechanism files and evaluating
// reaction source terms.  The routines are intentionally lightweight to keep
// the focus on the ODE integrators.

#pragma once

#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cmath>

// Basic chemical reaction and thermodynamic data structures -----------------

// Elementary reaction: reactants and products are stored as (index, stoich)
// pairs.  The Arrhenius parameters follow the common A * T^b * exp(-E/RT)
// convention.
template<typename T>
struct Reaction {
    std::vector<std::pair<int,T>> react; // species index and coefficient
    std::vector<std::pair<int,T>> prod;  // species index and coefficient
    T A{}, b{}, E{};                     // Arrhenius parameters
};

// NASA polynomial coefficients used for evaluating specific heats etc.
template<typename T>
struct ThermoData{
    T t_low{}, t_mid{}, t_high{};        // Temperature bounds
    T high[7]{};                         // High-temperature coefficients
    T low[7]{};                          // Low-temperature coefficients
};

// Helper to trim whitespace from both ends of a string
inline std::string trim(const std::string& s){
    size_t b = s.find_first_not_of(" \t\r\n");
    if(b==std::string::npos) return "";
    size_t e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e-b+1);
}

// Parse one side of a reaction expression such as "2H2 + O2".  Each species
// token is converted to an index using `idx` and stored with its coefficient.
template<typename T>
inline std::vector<std::pair<int,T>>
parse_side(const std::string& side, const std::map<std::string,int>& idx){
    std::vector<std::pair<int,T>> res;
    std::string token;
    for(size_t i=0;i<=side.size();++i){
        if(i==side.size() || side[i]=='+'){
            if(!token.empty()){
                size_t p = token.find('(');
                if(p!=std::string::npos) token = token.substr(0,p);
                if(token=="M" || token.empty()){
                    token.clear();
                    continue;           // Ignore third bodies and blanks
                }
                T coeff = static_cast<T>(1.0);
                size_t j=0;
                while(j<token.size() && (std::isdigit(token[j])||token[j]=='.')) j++;
                if(j>0) coeff = static_cast<T>(std::stod(token.substr(0,j)));
                std::string name = token.substr(j);
                auto it = idx.find(name);
                if(it!=idx.end()) res.push_back({it->second, coeff});
                token.clear();
            }
        } else if(!std::isspace(static_cast<unsigned char>(side[i])))
            token.push_back(side[i]);
    }
    return res;
}

// Read a very small subset of the CHEMKIN format: species names and simple
// reactions with Arrhenius parameters.  Lines starting with '!' are treated as
// comments and ignored.  Reaction arrows "=>", "<=>" and the common reversible
// form "=" are recognized.
template<typename T>
inline bool load_chemkin(const std::string& fname,
                         std::vector<std::string>& species,
                         std::map<std::string,int>& idx,
                         std::vector<Reaction<T>>& reactions){
    std::ifstream ifs(fname);
    if(!ifs) return false;
    species.clear();
    idx.clear();
    reactions.clear();
    std::string line;
    enum Section{NONE,SPECIES,REACTIONS};
    Section sec=NONE;
    while(std::getline(ifs,line)){
        line = trim(line);
        if(line.empty() || line[0]=='!') continue; // Skip comments

        if(line.rfind("SPECIES",0)==0){
            sec = SPECIES;
            std::string rest = trim(line.substr(7));
            if(!rest.empty()){
                std::stringstream ss(rest);
                std::string sp;
                while(ss>>sp){
                    if(idx.find(sp)==idx.end()){
                        idx[sp]=species.size();
                        species.push_back(sp);
                    }
                }
            }
            continue;
        }
        if(line.rfind("REACTIONS",0)==0){ sec=REACTIONS; continue; }
        if(line=="END"){ sec=NONE; continue; }

        if(sec==SPECIES){
            std::stringstream ss(line);
            std::string sp;
            while(ss>>sp){
                if(idx.find(sp)==idx.end()){
                    idx[sp]=species.size();
                    species.push_back(sp);
                }
            }
        } else if(sec==REACTIONS){
            if(line.find('=')==std::string::npos) continue;
            std::string uline = line;
            std::transform(uline.begin(), uline.end(), uline.begin(),
                           [](unsigned char c){ return std::toupper(c); });
            if(uline.rfind("LOW",0)==0 || uline.rfind("TROE",0)==0 ||
               uline.rfind("HIGH",0)==0) continue; // ignore falloff
            if(uline.rfind("DUP",0)==0) continue;   // ignore duplicates
            size_t excl = line.find('!');
            std::string nocmt = (excl==std::string::npos)? line : line.substr(0,excl);
            std::stringstream ss(nocmt);
            std::string expr; double A,b,E;
            if(!(ss>>expr>>A>>b>>E)) continue;
            std::string lhs, rhs;
            size_t pos;
            if((pos = expr.find("<=>")) != std::string::npos){
                lhs = expr.substr(0,pos);
                rhs = expr.substr(pos+3);
            } else if((pos = expr.find("=>")) != std::string::npos){
                lhs = expr.substr(0,pos);
                rhs = expr.substr(pos+2);
            } else if((pos = expr.find('=')) != std::string::npos){
                lhs = expr.substr(0,pos);
                rhs = expr.substr(pos+1);
            } else continue;
            Reaction<T> r;
            r.react = parse_side<T>(lhs, idx);
            r.prod  = parse_side<T>(rhs, idx);
            r.A = static_cast<T>(A);
            r.b = static_cast<T>(b);
            r.E = static_cast<T>(E);
            reactions.push_back(r);
        }
    }
    return true;
}

// Load NASA polynomial thermodynamic coefficients.
template<typename T>
inline bool load_thermo(const std::string& fname,
                        const std::map<std::string,int>& idx,
                        std::vector<ThermoData<T>>& thermo){
    std::ifstream ifs(fname);
    if(!ifs) return false;
    std::string line;
    T t_low = 0, t_mid = 0, t_high = 0;
    while(std::getline(ifs,line)){
        line = trim(line);
        if(line.empty()) continue;
        if(line.rfind("THERMO",0)==0){
            if(std::getline(ifs,line)){
                std::stringstream ss(line);
                ss >> t_low >> t_mid >> t_high;
            }
            break;
        }
    }
    if(t_high <= t_low) return false;
    thermo.assign(idx.size(), ThermoData<T>{t_low, t_mid, t_high, {}, {}});
    while(std::getline(ifs,line)){
        if(line.rfind("END",0)==0) break;
        if(line.size()<1) continue;
        std::string name = trim(line.substr(0,16));
        std::string l2,l3,l4;
        if(!std::getline(ifs,l2)) break;
        if(!std::getline(ifs,l3)) break;
        if(!std::getline(ifs,l4)) break;
        if(l2.size()>=2) l2.erase(l2.size()-2);
        if(l3.size()>=2) l3.erase(l3.size()-2);
        if(l4.size()>=2) l4.erase(l4.size()-2);
        std::stringstream ss;
        ss << l2 << " " << l3 << " " << l4;
        std::vector<T> coef; T val;
        while(ss >> val) coef.push_back(val);
        if(coef.size()!=14) continue;
        auto it = idx.find(name);
        if(it!=idx.end()){
            auto& td = thermo[it->second];
            std::copy(coef.begin(), coef.begin()+7, td.high);
            std::copy(coef.begin()+7, coef.end(), td.low);
        }
    }
    return true;
}

// Compute d/dt of species mole fractions and temperature given current state.
// A more complete energy equation is used where temperature evolution is
// driven by species enthalpies obtained from NASA polynomial coefficients
// in `thermo`.
template<typename T>
inline void compute_rhs(const std::vector<Reaction<T>>& reactions,
                        const std::vector<ThermoData<T>>& thermo,
                        T P, const std::vector<T>& y, std::vector<T>& dy){
    const T R = static_cast<T>(8.3144621);
    size_t n = y.size()-1;                 // last entry is temperature
    T Tval = y[n];
    if(!std::isfinite(Tval) || Tval < T(1e-6)) Tval = T(1e-6);
    if(Tval > T(1e4)) Tval = T(1e4);

    // initialise derivatives
    for(size_t i=0;i<=n;++i) dy[i]=T(0);

    // species concentrations (ideal gas, constant pressure)
    std::vector<T> conc(n);
    for(size_t i=0;i<n;++i){
        T yi = y[i];
        if(!std::isfinite(yi) || yi < T(0)) yi = T(0);
        if(yi > T(1)) yi = T(1);
        conc[i] = yi*P/(R*Tval);
    }

    // reaction source terms for species
    for(const auto& r : reactions){
        T k = r.A*std::exp(r.b*std::log(Tval)-r.E/(R*Tval));
        T rate = k;
        for(auto [idx,nu] : r.react) rate *= std::pow(conc[idx],nu);
        if(!std::isfinite(rate)) continue;
        const T max_rate = T(1e6);
        if(rate > max_rate) rate = max_rate;
        if(rate < -max_rate) rate = -max_rate;
        for(auto [idx,nu] : r.react) dy[idx] -= nu*rate;
        for(auto [idx,nu] : r.prod)  dy[idx] += nu*rate;
    }

    // prevent negative time derivatives when species is already depleted
    for(size_t i=0;i<n;++i)
        if(y[i] <= T(0) && dy[i] < T(0)) dy[i] = T(0);

    // Temperature equation using species enthalpies from thermo data
    T sumY = T(0);
    for(size_t i=0;i<n;++i) sumY += y[i];
    T cp_mix = T(0);
    T heat = T(0);
    for(size_t i=0;i<n;++i){
        const auto& td = thermo[i];
        const T* c = (Tval >= td.t_mid) ? td.high : td.low;
        T t = Tval;
        T t2 = t*t;
        T t3 = t2*t;
        T t4 = t3*t;
        // cp and enthalpy (J/mol/K and J/mol)
        T cp_i = R*(c[0] + c[1]*t + c[2]*t2 + c[3]*t3 + c[4]*t4);
        T h_i = R*t*(c[0] + c[1]*t/T(2) + c[2]*t2/T(3) +
                     c[3]*t3/T(4) + c[4]*t4/T(5) + c[5]/t);
        T Xi = (sumY>0)? y[i]/sumY : T(0);
        cp_mix += Xi * cp_i;
        heat += h_i * dy[i];
    }
    T Ctot = P/(R*Tval); // total molar concentration
    if(cp_mix <= T(0)) cp_mix = T(1); // prevent divide-by-zero
    dy[n] = -heat / (Ctot * cp_mix);
}

