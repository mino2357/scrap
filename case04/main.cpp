#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <cctype>
#include <cmath>

struct Reaction {
    std::vector<std::pair<int,double>> react;
    std::vector<std::pair<int,double>> prod;
    double A{}, b{}, E{};
};

static std::string trim(const std::string& s){
    size_t b = s.find_first_not_of(" \t\r\n");
    if(b==std::string::npos) return "";
    size_t e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e-b+1);
}

static std::vector<std::pair<int,double>>
parse_side(const std::string& side, const std::map<std::string,int>& idx){
    std::vector<std::pair<int,double>> res;
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
                double coeff = 1.0;
                size_t j=0;
                while(j<token.size() && (std::isdigit(token[j])||token[j]=='.')) j++;
                if(j>0) coeff = std::stod(token.substr(0,j));
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

static bool load_chemkin(const std::string& fname,
                         std::vector<std::string>& species,
                         std::map<std::string,int>& idx,
                         std::vector<Reaction>& reactions){
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
        if(line=="SPECIES"){ sec=SPECIES; continue; }
        if(line=="REACTIONS"){ sec=REACTIONS; continue; }
        if(line=="END"){ sec=NONE; continue; }
        if(sec==SPECIES){
            std::stringstream ss(line);
            std::string sp;
            while(ss>>sp){
                idx[sp]=species.size();
                species.push_back(sp);
            }
        } else if(sec==REACTIONS){
            if(line.find('=')==std::string::npos) continue;
            if(line.find("LOW /")!=std::string::npos || line.find("TROE")!=std::string::npos) continue;
            if(line=="DUP") continue;
            size_t excl = line.find('!');
            std::string nocmt = (excl==std::string::npos)? line : line.substr(0,excl);
            std::stringstream ss(nocmt);
            std::string expr; double A,b,E;
            if(!(ss>>expr>>A>>b>>E)) continue;
            size_t eq = expr.find('=');
            std::string lhs = expr.substr(0,eq);
            std::string rhs = expr.substr(eq+1);
            Reaction r;
            r.react = parse_side(lhs, idx);
            r.prod  = parse_side(rhs, idx);
            r.A=A; r.b=b; r.E=E;
            reactions.push_back(r);
        }
    }
    return true;
}

static void compute_rates(const std::vector<Reaction>& reactions, double T,
                          const std::vector<double>& c, std::vector<double>& dc){
    std::fill(dc.begin(), dc.end(), 0.0);
    for(const auto& r: reactions){
        double k = r.A * std::pow(T, r.b) * std::exp(-r.E / T);
        double rate = k;
        for(auto [i,nu]: r.react)
            rate *= std::pow(std::max(c[i],0.0), nu);
        for(auto [i,nu]: r.react)
            dc[i] -= nu*rate;
        for(auto [i,nu]: r.prod)
            dc[i] += nu*rate;
    }
}

static void load_initial(const std::string& fname,
                         const std::map<std::string,int>& idx,
                         std::vector<double>& c){
    std::ifstream ifs(fname);
    if(!ifs) return;
    std::string name; double val;
    while(ifs >> name >> val){
        auto it = idx.find(name);
        if(it != idx.end()) c[it->second] = val;
    }
}

using Integrator = void(*)(std::vector<double>&, double, double, double, const std::vector<Reaction>&);

static void euler(std::vector<double>& y, double t0, double t1, double T,
                  const std::vector<Reaction>& reactions){
    double h = (t1 - t0) / 1000.0;
    double t = t0;
    std::vector<double> dy(y.size());
    while(t < t1){
        if(t + h > t1) h = t1 - t;
        compute_rates(reactions, T, y, dy);
        for(size_t i=0;i<y.size();++i) y[i] += h*dy[i];
        for(auto& v : y) if(v < 0) v = 0;
        t += h;
    }
}

static void rk45(std::vector<double>& y, double t0, double t1, double T,
                 const std::vector<Reaction>& reactions){
    const double tol = 1e-8;
    double h = (t1 - t0) / 1000.0;
    double t = t0;
    size_t n = y.size();
    std::vector<double> k1(n),k2(n),k3(n),k4(n),k5(n),k6(n),yt(n),y4(n),y5(n);
    while(t < t1){
        if(t + h > t1) h = t1 - t;
        compute_rates(reactions, T, y, k1);
        for(size_t i=0;i<n;++i) yt[i] = y[i] + h*(1.0/5.0*k1[i]);
        compute_rates(reactions, T, yt, k2);
        for(size_t i=0;i<n;++i) yt[i] = y[i] + h*(3.0/40.0*k1[i] + 9.0/40.0*k2[i]);
        compute_rates(reactions, T, yt, k3);
        for(size_t i=0;i<n;++i) yt[i] = y[i] + h*(3.0/10.0*k1[i] - 9.0/10.0*k2[i] + 6.0/5.0*k3[i]);
        compute_rates(reactions, T, yt, k4);
        for(size_t i=0;i<n;++i) yt[i] = y[i] + h*(-11.0/54.0*k1[i] + 2.5*k2[i] - 70.0/27.0*k3[i] + 35.0/27.0*k4[i]);
        compute_rates(reactions, T, yt, k5);
        for(size_t i=0;i<n;++i) yt[i] = y[i] + h*(1631.0/55296.0*k1[i] + 175.0/512.0*k2[i] + 575.0/13824.0*k3[i] + 44275.0/110592.0*k4[i] + 253.0/4096.0*k5[i]);
        compute_rates(reactions, T, yt, k6);
        double err = 0.0;
        for(size_t i=0;i<n;++i){
            y5[i] = y[i] + h*(37.0/378.0*k1[i] + 250.0/621.0*k3[i] + 125.0/594.0*k4[i] + 512.0/1771.0*k6[i]);
            y4[i] = y[i] + h*(2825.0/27648.0*k1[i] + 18575.0/48384.0*k3[i] + 13525.0/55296.0*k4[i] + 277.0/14336.0*k5[i] + 0.25*k6[i]);
            err = std::max(err, std::abs(y5[i]-y4[i]));
        }
        if(err <= tol){
            y = y5;
            for(auto& v : y) if(v < 0) v = 0;
            t += h;
        }
        double scale = (err==0.0 ? 2.0 : 0.9*std::pow(tol/err, 0.2));
        scale = std::min(5.0, std::max(0.2, scale));
        h *= scale;
    }
}

int main(int argc, char** argv){
    std::vector<std::string> species;
    std::map<std::string,int> idx;
    std::vector<Reaction> reactions;
    if(!load_chemkin("chem.inp", species, idx, reactions)){
        std::cerr << "chem.inp not found\n";
        return 1;
    }

    const double T = 1000.0;
    std::vector<double> c(species.size(),0.0);
    load_initial("init.inp", idx, c);

    Integrator integrator = rk45;
    if(argc>1 && std::string(argv[1])=="euler")
        integrator = euler;

    integrator(c, 0.0, 1e-3, T, reactions);

    for(size_t i=0;i<species.size();++i)
        std::cout<<species[i]<<" "<<c[i]<<"\n";
    return 0;
}
