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
static void compute_rates(const std::vector<Reaction<T>>& reactions, T Tgas,
                          const std::vector<T>& c, std::vector<T>& dc){
    std::fill(dc.begin(), dc.end(), T(0));
    for(const auto& r: reactions){
        T k = r.A * std::pow(Tgas, r.b) * std::exp(-r.E / Tgas);
        T rate = k;
        for(auto [i,nu]: r.react)
            rate *= std::pow(std::max(c[i], T(0)), nu);
        for(auto [i,nu]: r.react)
            dc[i] -= nu*rate;
        for(auto [i,nu]: r.prod)
            dc[i] += nu*rate;
    }
}


template<typename T>
using Callback = std::function<void(T,const std::vector<T>&)>;
template<typename T>
using Integrator = void(*)(std::vector<T>&, T, T, T,
                          const std::vector<Reaction<T>>&, Callback<T>);

template<typename T>
static void euler(std::vector<T>& y, T t0, T t1, T Tgas,
                  const std::vector<Reaction<T>>& reactions, Callback<T> cb){
    T h = (t1 - t0) / T(1000);
    T t = t0;
    std::vector<T> dy(y.size());
    while(t < t1){
        if(t + h > t1) h = t1 - t;
        compute_rates(reactions, Tgas, y, dy);
        for(size_t i=0;i<y.size();++i) y[i] += h*dy[i];
        for(auto& v : y) if(v < T(0)) v = T(0);
        t += h;
        if(cb) cb(t, y);
    }
}

template<typename T>
static void rk45(std::vector<T>& y, T t0, T t1, T Tgas,
                 const std::vector<Reaction<T>>& reactions, Callback<T> cb){
    const T tol = T(1e-8);
    T h = (t1 - t0) / T(1000);
    T t = t0;
    size_t n = y.size();
    std::vector<T> k1(n),k2(n),k3(n),k4(n),k5(n),k6(n),yt(n),y4(n),y5(n);
    while(t < t1){
        if(t + h > t1) h = t1 - t;
        compute_rates(reactions, Tgas, y, k1);
        for(size_t i=0;i<n;++i) yt[i] = y[i] + h*(T(1.0/5.0)*k1[i]);
        compute_rates(reactions, Tgas, yt, k2);
        for(size_t i=0;i<n;++i) yt[i] = y[i] + h*(T(3.0/40.0)*k1[i] + T(9.0/40.0)*k2[i]);
        compute_rates(reactions, Tgas, yt, k3);
        for(size_t i=0;i<n;++i) yt[i] = y[i] + h*(T(3.0/10.0)*k1[i] - T(9.0/10.0)*k2[i] + T(6.0/5.0)*k3[i]);
        compute_rates(reactions, Tgas, yt, k4);
        for(size_t i=0;i<n;++i) yt[i] = y[i] + h*(T(-11.0/54.0)*k1[i] + T(2.5)*k2[i] - T(70.0/27.0)*k3[i] + T(35.0/27.0)*k4[i]);
        compute_rates(reactions, Tgas, yt, k5);
        for(size_t i=0;i<n;++i) yt[i] = y[i] + h*(T(1631.0/55296.0)*k1[i] + T(175.0/512.0)*k2[i] + T(575.0/13824.0)*k3[i] + T(44275.0/110592.0)*k4[i] + T(253.0/4096.0)*k5[i]);
        compute_rates(reactions, Tgas, yt, k6);
        T err = T(0);
        for(size_t i=0;i<n;++i){
            y5[i] = y[i] + h*(T(37.0/378.0)*k1[i] + T(250.0/621.0)*k3[i] + T(125.0/594.0)*k4[i] + T(512.0/1771.0)*k6[i]);
            y4[i] = y[i] + h*(T(2825.0/27648.0)*k1[i] + T(18575.0/48384.0)*k3[i] + T(13525.0/55296.0)*k4[i] + T(277.0/14336.0)*k5[i] + T(0.25)*k6[i]);
            err = std::max(err, std::abs(y5[i]-y4[i]));
        }
        if(err <= tol){
            y = y5;
            for(auto& v : y) if(v < T(0)) v = T(0);
            t += h;
            if(cb) cb(t, y);
        }
        T scale = (err==T(0) ? T(2) : T(0.9)*std::pow(tol/err, T(0.2)));
        scale = std::min(T(5), std::max(T(0.2), scale));
        h *= scale;
    }
}

template<typename T>
static void rk78(std::vector<T>& y, T t0, T t1, T Tgas,
                 const std::vector<Reaction<T>>& reactions, Callback<T> cb){
    const T tol = T(1e-8);
    T h = (t1 - t0) / T(1000);
    T t = t0;
    size_t n = y.size();
    std::vector<T> k1(n),k2(n),k3(n),k4(n),k5(n),k6(n),k7(n),k8(n),
                   k9(n),k10(n),k11(n),k12(n),k13(n),yt(n),y8(n);
    while(t < t1){
        if(t + h > t1) h = t1 - t;

        compute_rates(reactions, Tgas, y, k1);
        for(size_t i=0;i<n;++i) yt[i] = y[i] + h*(T(2.0/27.0)*k1[i]);
        compute_rates(reactions, Tgas, yt, k2);

        for(size_t i=0;i<n;++i) yt[i] = y[i] + h*(T(1.0/36.0)*k1[i] + T(1.0/12.0)*k2[i]);
        compute_rates(reactions, Tgas, yt, k3);

        for(size_t i=0;i<n;++i) yt[i] = y[i] + h*(T(1.0/24.0)*k1[i] + T(1.0/8.0)*k3[i]);
        compute_rates(reactions, Tgas, yt, k4);

        for(size_t i=0;i<n;++i) yt[i] = y[i] + h*(T(5.0/12.0)*k1[i] - T(25.0/16.0)*k3[i] + T(25.0/16.0)*k4[i]);
        compute_rates(reactions, Tgas, yt, k5);

        for(size_t i=0;i<n;++i) yt[i] = y[i] + h*(T(1.0/20.0)*k1[i] + T(1.0/4.0)*k4[i] + T(1.0/5.0)*k5[i]);
        compute_rates(reactions, Tgas, yt, k6);

        for(size_t i=0;i<n;++i) yt[i] = y[i] + h*(T(-25.0/108.0)*k1[i] + T(125.0/108.0)*k4[i]
                                               - T(65.0/27.0)*k5[i] + T(125.0/54.0)*k6[i]);
        compute_rates(reactions, Tgas, yt, k7);

        for(size_t i=0;i<n;++i) yt[i] = y[i] + h*(T(31.0/300.0)*k1[i] + T(61.0/225.0)*k5[i]
                                               - T(2.0/9.0)*k6[i] + T(13.0/900.0)*k7[i]);
        compute_rates(reactions, Tgas, yt, k8);

        for(size_t i=0;i<n;++i) yt[i] = y[i] + h*(T(2)*k1[i] - T(53.0/6.0)*k4[i]
                                               + T(704.0/45.0)*k5[i] - T(107.0/9.0)*k6[i]
                                               + T(67.0/90.0)*k7[i] + T(3)*k8[i]);
        compute_rates(reactions, Tgas, yt, k9);

        for(size_t i=0;i<n;++i) yt[i] = y[i] + h*(T(-91.0/108.0)*k1[i] + T(23.0/108.0)*k4[i]
                                               - T(976.0/135.0)*k5[i] + T(311.0/54.0)*k6[i]
                                               - T(19.0/60.0)*k7[i] + T(17.0/6.0)*k8[i]
                                               - T(1.0/12.0)*k9[i]);
        compute_rates(reactions, Tgas, yt, k10);

        for(size_t i=0;i<n;++i) yt[i] = y[i] + h*(T(2383.0/4100.0)*k1[i]
                                               - T(341.0/164.0)*k4[i] + T(4496.0/1025.0)*k5[i]
                                               - T(301.0/82.0)*k6[i] + T(2133.0/4100.0)*k7[i]
                                               + T(45.0/82.0)*k8[i] + T(45.0/164.0)*k9[i]
                                               + T(18.0/41.0)*k10[i]);
        compute_rates(reactions, Tgas, yt, k11);

        for(size_t i=0;i<n;++i) yt[i] = y[i] + h*(T(3.0/205.0)*k1[i] - T(6.0/41.0)*k6[i]
                                               - T(3.0/205.0)*k7[i] - T(3.0/41.0)*k8[i]
                                               + T(3.0/41.0)*k9[i] + T(6.0/41.0)*k10[i]);
        compute_rates(reactions, Tgas, yt, k12);

        for(size_t i=0;i<n;++i) yt[i] = y[i] + h*(T(-1777.0/4100.0)*k1[i]
                                               - T(341.0/164.0)*k4[i] + T(4496.0/1025.0)*k5[i]
                                               - T(289.0/82.0)*k6[i] + T(2193.0/4100.0)*k7[i]
                                               + T(51.0/82.0)*k8[i] + T(33.0/164.0)*k9[i]
                                               + T(12.0/41.0)*k10[i] + k12[i]);
        compute_rates(reactions, Tgas, yt, k13);

        T err = T(0);
        for(size_t i=0;i<n;++i){
            y8[i] = y[i] + h*(T(34.0/105.0)*k6[i] + T(9.0/35.0)*(k7[i]+k8[i])
                               + T(9.0/280.0)*(k9[i]+k10[i])
                               + T(41.0/840.0)*(k12[i]+k13[i]));
            T e = h*(T(-41.0/840.0)*k1[i] + T(-41.0/840.0)*k11[i]
                     + T(41.0/840.0)*k12[i] + T(41.0/840.0)*k13[i]);
            err = std::max(err, std::abs(e));
        }
        if(err <= tol){
            y = y8;
            for(auto& v : y) if(v < T(0)) v = T(0);
            t += h;
            if(cb) cb(t, y);
        }
        T scale = (err==T(0) ? T(2) : T(0.9)*std::pow(tol/err, T(1.0/8.0)));
        scale = std::min(T(5), std::max(T(0.2), scale));
        h *= scale;
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

    const T Tgas = T(1000.0);
    std::vector<T> c(species.size(), T(1e-8));
    if(!species.empty()) c[0] = T(1.0);
    if(species.size() > 1) c[1] = T(0.5);

    // Normalize initial composition to obtain mole fractions
    T c_sum = T(0);
    for(T v : c) c_sum += v;
    if(c_sum > T(0)){
        for(T& v : c) v /= c_sum;
    }

    Integrator<T> integrator = rk45<T>;
    if(argc>1){
        std::string method = argv[1];
        if(method=="euler")
            integrator = euler<T>;
        else if(method=="rk78")
            integrator = rk78<T>;
    }

    std::vector<T> times{T(0)};
    std::vector<std::vector<T>> history{c};
    const T output_interval = T(1e-5);
    T next_output = output_interval;
    auto cb = [&](T t,const std::vector<T>& y){
        while(t >= next_output){
            times.push_back(next_output);
            history.push_back(y);
            next_output += output_interval;
        }
    };

    integrator(c, T(0), T(1e-3), Tgas, reactions, cb);

    std::ofstream ofs("case04.dat");
    ofs << "time";
    for(const auto& s : species) ofs << ' ' << s;
    ofs << '\n';
    for(size_t k=0;k<times.size();++k){
        T sum=T(0); for(T v: history[k]) sum+=v;
        ofs << times[k];
        for(T v: history[k]) ofs << ' ' << (sum>0? v/sum : T(0));
        ofs << '\n';
    }

    for(size_t i=0;i<species.size();++i)
        std::cout<<species[i]<<" "<<c[i]<<"\n";
    return 0;
}
