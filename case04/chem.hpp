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
// irreversible reactions with Arrhenius parameters.  Lines starting with '!' are
// treated as comments and ignored.
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
// A very coarse constant-
// Cp energy equation is used for temperature evolution.
template<typename T>
inline void compute_rhs(const std::vector<Reaction<T>>& reactions,
                        const std::vector<ThermoData<T>>& /*thermo*/,
                        T P, const std::vector<T>& y, std::vector<T>& dy){
    const T R = static_cast<T>(8.3144621);
    size_t n = y.size()-1;
    T Tval = y[n];
    for(size_t i=0;i<=n;++i) dy[i]=T(0);
    std::vector<T> conc(n);
    for(size_t i=0;i<n;++i) conc[i] = y[i]*P/(R*Tval);
    for(const auto& r : reactions){
        T k = r.A*std::exp(r.b*std::log(Tval)-r.E/(R*Tval));
        T rate = k;
        for(auto [idx,nu] : r.react) rate *= std::pow(conc[idx],nu);
        for(auto [idx,nu] : r.react) dy[idx] -= nu*rate;
        for(auto [idx,nu] : r.prod)  dy[idx] += nu*rate;
    }
    // Temperature equation: assume constant Cp
    const T Cp = static_cast<T>(1000.0);
    T omega=0;
    for(size_t i=0;i<n;++i) omega += dy[i];
    dy[n] = -omega*R*Tval/Cp;
}

