#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <functional>

template<typename T>
struct Reaction {
    std::vector<std::pair<int,T>> react;
    std::vector<std::pair<int,T>> prod;
    T A{}, b{}, E{};
};

template<typename T>
struct ThermoData{
    T t_low{}, t_mid{}, t_high{};
    T high[7]{};
    T low[7]{};
};

static std::string trim(const std::string& s){
    size_t b = s.find_first_not_of(" \t\r\n");
    if(b==std::string::npos) return "";
    size_t e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e-b+1);
}

template<typename T>
static std::vector<std::pair<int,T>>
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
                    continue;
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

template<typename T>
static bool load_chemkin(const std::string& fname,
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
        if(line.empty() || line[0]=='!') continue;

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
               uline.rfind("HIGH",0)==0) continue;
            if(uline.rfind("DUP",0)==0) continue;
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

template<typename T>
static bool load_thermo(const std::string& fname,
                        const std::map<std::string,int>& idx,
                        std::vector<ThermoData<T>>& thermo){
    std::ifstream ifs(fname);
    if(!ifs) return false;
    std::string line;
    T t_low=0, t_mid=0, t_high=0;
    while(std::getline(ifs,line)){
        line = trim(line);
        if(line.empty()) continue;
        if(line.rfind("THERMO",0)==0){
            if(std::getline(ifs,line)){
                std::stringstream ss(line);
                ss>>t_low>>t_mid>>t_high;
            }
            break;
        }
    }
    if(t_high<=t_low) return false;
    thermo.assign(idx.size(), ThermoData<T>{t_low,t_mid,t_high,{},{}});
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
        ss<<l2<<" "<<l3<<" "<<l4;
        std::vector<T> coef; T val;
        while(ss>>val) coef.push_back(val);
        if(coef.size()!=14) continue;
        auto it = idx.find(name);
        if(it==idx.end()) continue;
        ThermoData<T> td; td.t_low=t_low; td.t_mid=t_mid; td.t_high=t_high;
        for(int i=0;i<7;++i) td.high[i]=coef[i];
        for(int i=0;i<7;++i) td.low[i]=coef[7+i];
        thermo[it->second]=td;
    }
    return true;
}

template<typename T>
static T cp_from_thermo(const ThermoData<T>& td, T Tgas){
    const T R = T(8.314462618);
    const T* a = (Tgas>td.t_mid)? td.high : td.low;
    T t=Tgas;
    return R*(a[0] + a[1]*t + a[2]*t*t + a[3]*t*t*t + a[4]*t*t*t*t);
}

template<typename T>
static T h_from_thermo(const ThermoData<T>& td, T Tgas){
    const T R = T(8.314462618);
    const T* a = (Tgas>td.t_mid)? td.high : td.low;
    T t=Tgas;
    return R*t*(a[0] + a[1]*t/T(2) + a[2]*t*t/T(3) + a[3]*t*t*t/T(4) + a[4]*t*t*t*t/T(5) + a[5]/t);
}

template<typename T>
static T s_from_thermo(const ThermoData<T>& td, T Tgas){
    const T R = T(8.314462618);
    const T* a = (Tgas>td.t_mid)? td.high : td.low;
    T t=Tgas;
    return R*(a[0]*std::log(t) + a[1]*t + a[2]*t*t/T(2) + a[3]*t*t*t/T(3) + a[4]*t*t*t*t/T(4) + a[6]);
}

template<typename T>
static void compute_rates_conc(const std::vector<Reaction<T>>& reactions,
                               const std::vector<ThermoData<T>>& thermo,
                               T Tgas,
                               const std::vector<T>& c,
                               std::vector<T>& dc){
    const T Rcal = T(1.987); // cal/mol/K
    const T R = T(8.314462618);
    std::fill(dc.begin(), dc.end(), T(0));
    std::vector<T> g(thermo.size());
    for(size_t i=0;i<thermo.size();++i){
        T h = h_from_thermo(thermo[i], Tgas);
        T s = s_from_thermo(thermo[i], Tgas);
        g[i] = h - Tgas*s;
    }
    for(const auto& r: reactions){
        T kf = r.A * std::pow(Tgas, r.b) * std::exp(-r.E / (Rcal*Tgas));
        T conc_f = T(1), conc_r = T(1);
        T delta_nu = T(0);
        T dg = T(0);
        for(auto [i,nu]: r.react){
            conc_f *= std::pow(std::max(c[i], T(0)), nu);
            delta_nu -= nu;
            dg -= nu * g[i];
        }
        for(auto [i,nu]: r.prod){
            conc_r *= std::pow(std::max(c[i], T(0)), nu);
            delta_nu += nu;
            dg += nu * g[i];
        }
        T Kc = std::exp(-dg/(R*Tgas)) * std::pow(R*Tgas*T(1e-6), delta_nu);
        T kr = (Kc>0)? kf / Kc : T(0);
        T rate = kf*conc_f - kr*conc_r;
        for(auto [i,nu]: r.react) dc[i] -= nu*rate;
        for(auto [i,nu]: r.prod)  dc[i] += nu*rate;
    }
}

template<typename T>
static void compute_rates(const std::vector<Reaction<T>>& reactions,
                          const std::vector<ThermoData<T>>& thermo,
                          T Tgas, T P,
                          const std::vector<T>& X, std::vector<T>& omega){
    const T R = T(8.314462618);
    T Ctot = P/(R*Tgas);
    T Ctot_c = Ctot/T(1e6); // mol/cm^3
    std::vector<T> c(X.size());
    for(size_t i=0;i<X.size();++i) c[i] = X[i]*Ctot_c;
    compute_rates_conc(reactions, thermo, Tgas, c, omega);
    for(size_t i=0;i<omega.size();++i) omega[i] *= T(1e6); // back to mol/m^3/s
}
// Evaluate species rates and temperature derivative.
template<typename T>
static void compute_rhs(const std::vector<Reaction<T>>& reactions,
                        const std::vector<ThermoData<T>>& thermo,
                        T P,
                        const std::vector<T>& y,
                        std::vector<T>& dy){
    const size_t n = thermo.size();
    std::vector<T> X(y.begin(), y.begin()+n);
    std::vector<T> omega(n);
    T Tgas = std::max(y[n], T(1));
    compute_rates(reactions, thermo, Tgas, P, X, omega);
    const T R = T(8.314462618);
    T Ctot = P/(R*Tgas);
    T cp_mix = T(0);
    T enth_rate = T(0);
    for(size_t i=0;i<n;++i){
        T cp = cp_from_thermo(thermo[i], Tgas);
        T h  = h_from_thermo(thermo[i], Tgas);
        cp_mix += X[i]*cp;
        enth_rate += h*omega[i];
    }
    dy[n] = (cp_mix>0)? -enth_rate/(Ctot*cp_mix) : T(0);
    for(size_t i=0;i<n;++i)
        dy[i] = omega[i]/Ctot + X[i]/Tgas * dy[n];
}


template<typename T>
using Callback = std::function<void(T,const std::vector<T>&)>;
template<typename T>
using Integrator = void(*)(std::vector<T>&, T, T, T,
                          const std::vector<Reaction<T>>&, const std::vector<ThermoData<T>>&, Callback<T>);

template<typename T>
static void rk4(std::vector<T>& y, T t0, T t1, T P,
                const std::vector<Reaction<T>>& reactions,
                const std::vector<ThermoData<T>>& thermo,
                Callback<T> cb){
    T h = (t1 - t0) / T(1000);
    T t = t0;
    size_t m = y.size();
    std::vector<T> k1(m),k2(m),k3(m),k4(m),yt(m);
    while(t < t1){
        if(t + h > t1) h = t1 - t;
        compute_rhs(reactions, thermo, P, y, k1);
        for(size_t i=0;i<m;++i) yt[i] = y[i] + T(0.5)*h*k1[i];
        compute_rhs(reactions, thermo, P, yt, k2);
        for(size_t i=0;i<m;++i) yt[i] = y[i] + T(0.5)*h*k2[i];
        compute_rhs(reactions, thermo, P, yt, k3);
        for(size_t i=0;i<m;++i) yt[i] = y[i] + h*k3[i];
        compute_rhs(reactions, thermo, P, yt, k4);
        for(size_t i=0;i<m;++i)
            y[i] += (h/T(6))*(k1[i] + T(2)*k2[i] + T(2)*k3[i] + k4[i]);
        T sum = T(0);
        for(size_t i=0;i<m-1;++i){
            if(y[i] < T(0)) y[i] = T(0);
            sum += y[i];
        }
        if(sum>0) for(size_t i=0;i<m-1;++i) y[i] /= sum;
        t += h;
        if(cb) cb(t, y);
    }
}

// Runge-Kutta-Fehlberg 4(5) with adaptive step size.
template<typename T>
static void rk45(std::vector<T>& y, T t0, T t1, T P,
                 const std::vector<Reaction<T>>& reactions,
                 const std::vector<ThermoData<T>>& thermo,
                 Callback<T> cb){
    const T tol = T(1e-6);
    const T safety = T(0.9);
    T h = (t1 - t0) / T(1000);
    T t = t0;
    const size_t m = y.size();
    std::vector<T> k1(m),k2(m),k3(m),k4(m),k5(m),k6(m),yt(m),y4(m),y5(m),errv(m);
    while(t < t1){
        if(t + h > t1) h = t1 - t;
        T err;
        do{
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

            for(size_t i=0;i<m;++i){
                y4[i] = y[i] + h*( T(25.0)/T(216.0)*k1[i] + T(1408.0)/T(2565.0)*k3[i]
                                   + T(2197.0)/T(4104.0)*k4[i] - T(1.0)/T(5.0)*k5[i] );
                y5[i] = y[i] + h*( T(16.0)/T(135.0)*k1[i] + T(6656.0)/T(12825.0)*k3[i]
                                   + T(28561.0)/T(56430.0)*k4[i] - T(9.0)/T(50.0)*k5[i]
                                   + T(2.0)/T(55.0)*k6[i] );
                errv[i] = y5[i] - y4[i];
            }
            err = T(0);
            for(size_t i=0;i<m;++i) err = std::max(err, std::abs(errv[i]));
            if(err>tol){
                T fac = (err>0)? safety*std::pow(tol/err, T(0.2)) : T(0.5);
                fac = std::min(T(5.0), std::max(T(0.1), fac));
                h *= fac;
                if(t + h > t1) h = t1 - t;
            }
        }while(err>tol);

        y = y5;
        T sum = T(0);
        for(size_t i=0;i<m-1;++i){
            if(y[i] < T(0)) y[i] = T(0);
            sum += y[i];
        }
        if(sum>0) for(size_t i=0;i<m-1;++i) y[i] /= sum;
        t += h;
        if(cb) cb(t, y);

        T fac = (err>0)? safety*std::pow(tol/err, T(0.2)) : T(5.0);
        fac = std::min(T(5.0), std::max(T(0.1), fac));
        h *= fac;
    }
}

// Runge-Kutta-Fehlberg 7(8) with adaptive step size.
template<typename T>
static void rk78(std::vector<T>& y, T t0, T t1, T P,
                 const std::vector<Reaction<T>>& reactions,
                 const std::vector<ThermoData<T>>& thermo,
                 Callback<T> cb){
    const T tol = T(1e-7);
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
            compute_rhs(reactions, thermo, P, y, k[0]);
            for(int s=1; s<13; ++s){
                for(size_t i=0;i<m;++i){
                    yt[i] = y[i];
                    for(int j=0;j<s; ++j) yt[i] += h * a[s][j] * k[j][i];
                }
                compute_rhs(reactions, thermo, P, yt, k[s]);
            }
            for(size_t i=0;i<m;++i){
                T sum_b = T(0), sum_e = T(0);
                for(int s=0; s<13; ++s){
                    sum_b += b[s]*k[s][i];
                    sum_e += db[s]*k[s][i];
                }
                y8[i] = y[i] + h*sum_b;
                errv[i] = h*sum_e;
            }
            err = T(0);
            for(size_t i=0;i<m;++i) err = std::max(err, std::abs(errv[i]));
            if(err>tol){
                T fac = (err>0)? safety*std::pow(tol/err, T(1.0)/T(8.0)) : T(0.5);
                fac = std::min(T(4.0), std::max(T(0.1), fac));
                h *= fac;
                if(t + h > t1) h = t1 - t;
            }
        }while(err>tol);

        y = y8;
        T sum = T(0);
        for(size_t i=0;i<m-1;++i){
            if(y[i] < T(0)) y[i] = T(0);
            sum += y[i];
        }
        if(sum>0) for(size_t i=0;i<m-1;++i) y[i] /= sum;
        t += h;
        if(cb) cb(t, y);

        T fac = (err>0)? safety*std::pow(tol/err, T(1.0)/T(8.0)) : T(4.0);
        fac = std::min(T(4.0), std::max(T(0.1), fac));
        h *= fac;
    }
}

int main(int argc, char** argv){
    using T = double;
    std::vector<std::string> species;
    std::map<std::string,int> idx;
    std::vector<Reaction<T>> reactions;
    if(!load_chemkin<T>("chem.inp", species, idx, reactions)){
        std::cerr << "chem.inp not found\n";
        return 1;
    }
    std::vector<ThermoData<T>> thermo;
    if(!load_thermo<T>("therm.dat", idx, thermo)){
        std::cerr << "therm.dat not found\n";
        return 1;
    }

    // Initial composition: H2 1, O2 1, N2 3.76
    std::vector<T> X(species.size(), T(1e-8));
    auto set_init = [&](const std::string& name, T val){
        auto it = idx.find(name);
        if(it!=idx.end()) X[it->second] = val;
    };
    set_init("H2", T(1.0));
    set_init("O2", T(1.0));
    set_init("N2", T(3.76));

    T X_sum = T(0);
    for(T v : X) X_sum += v;
    if(X_sum > T(0)) for(T& v : X) v /= X_sum;

    const T P = T(202650.0); // Pa
    T T0 = T(1000.0);
    std::vector<T> y = X;
    y.push_back(T0);

    Integrator<T> integrator = rk4<T>;
    if(argc > 1){
        std::string method = argv[1];
        if(method == "rk45") integrator = rk45<T>;
        else if(method == "rk78") integrator = rk78<T>;
    }

    std::vector<T> times{T(0)};
    std::vector<std::vector<T>> history{y};
    const T output_interval = T(1e-5);
    T next_output = output_interval;
    auto cb = [&](T t,const std::vector<T>& ystate){
        while(t >= next_output){
            times.push_back(next_output);
            history.push_back(ystate);
            next_output += output_interval;
        }
    };

    integrator(y, T(0), T(1e-3), P, reactions, thermo, cb);

    const size_t n = species.size();
    std::ofstream ofs("case04.dat");
    ofs << "time";
    for(const auto& s : species) ofs << ' ' << s;
    ofs << " T\n";
    for(size_t k=0;k<times.size();++k){
        T sum=T(0); for(size_t i=0;i<n;++i) sum+=history[k][i];
        ofs << times[k];
        for(size_t i=0;i<n;++i) ofs << ' ' << (sum>0? history[k][i]/sum : T(0));
        ofs << ' ' << history[k][n] << '\n';
    }

    for(size_t i=0;i<n;++i)
        std::cout<<species[i]<<" "<<y[i]<<"\n";
    return 0;
}
