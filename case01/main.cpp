// 1D Euler (ideal gas), MUSCL + HLLC(+HLLE fallback), TVD-RK2
// Build: g++ -O3 -std=c++17 -Wall -Wextra -pedantic euler1d.cpp -o euler1d

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

// ---- Types ----
struct State { double rho, mom, E; };
struct Prim  { double rho, u,   p; };

struct Params {
    int    N = 800; double xL = 0.0, xR = 1.0; double CFL = 0.5;
    double t_end = 0.2; double gamma = 1.4; double p_floor = 1e-12, rho_floor = 1e-12;
    bool   use_muscl = true; bool output_every_step = false;
};

// ---- Utils ----
static inline double sqr(double x) noexcept { return x*x; }
static inline double sound_speed(const Prim& W, double g) noexcept { return std::sqrt(std::max(0.0, g*W.p/std::max(W.rho,1e-300))); }
static inline double minmod(double a,double b) noexcept { if(a*b<=0.0) return 0.0; return (std::abs(a)<std::abs(b))?a:b; }

static inline State prim2cons(const Prim& W, double g) noexcept {
    State U{std::max(W.rho,0.0), 0.0, 0.0};
    U.mom = U.rho*W.u; U.E = W.p/(g-1.0) + 0.5*U.rho*sqr(W.u); return U;
}
static inline Prim cons2prim(const State& U, double g, double pf, double rf) noexcept {
    Prim W{}; W.rho=std::max(U.rho, rf); W.u=U.mom/W.rho; const double kin=0.5*W.rho*sqr(W.u);
    W.p=std::max((g-1.0)*(U.E-kin), pf); return W;
}
static inline State flux_from_prim(const Prim& W, double g) noexcept {
    State F{}; const double E=W.p/(g-1.0)+0.5*W.rho*sqr(W.u); F.rho=W.rho*W.u; F.mom=W.rho*sqr(W.u)+W.p; F.E=(E+W.p)*W.u; return F;
}

// ---- Riemann solvers ----
static inline State hlle_flux(const Prim& WL, const Prim& WR, double g) noexcept {
    const State UL=prim2cons(WL,g), UR=prim2cons(WR,g);
    const State FL=flux_from_prim(WL,g), FR=flux_from_prim(WR,g);
    const double aL=sound_speed(WL,g), aR=sound_speed(WR,g);
    const double SL=std::min(WL.u - aL, WR.u - aR);
    const double SR=std::max(WL.u + aL, WR.u + aR);
    if (SL >= 0.0) return FL;
    if (SR <= 0.0) return FR;
    const double inv = 1.0/(SR - SL);
    return State{ (SR*FL.rho - SL*FR.rho + SL*SR*(UR.rho-UL.rho))*inv,
                  (SR*FL.mom - SL*FR.mom + SL*SR*(UR.mom-UL.mom))*inv,
                  (SR*FL.E   - SL*FR.E   + SL*SR*(UR.E  -UL.E  ))*inv };
}

static inline State hllc_flux(const Prim& WL, const Prim& WR, double g) noexcept {
    const double aL = sound_speed(WL,g), aR = sound_speed(WR,g);
    const double SL = std::min(WL.u - aL, WR.u - aR);
    const double SR = std::max(WL.u + aL, WR.u + aR);
    const State  UL = prim2cons(WL,g),  UR = prim2cons(WR,g);
    const State  FL = flux_from_prim(WL,g), FR = flux_from_prim(WR,g);
    if (SL >= 0.0) return FL;
    if (SR <= 0.0) return FR;

    const double num = (WR.p - WL.p) + WL.rho*WL.u*(SL - WL.u) - WR.rho*WR.u*(SR - WR.u);
    const double den = WL.rho*(SL - WL.u) - WR.rho*(SR - WR.u);
    const double Sstar = (std::abs(den) < 1e-14) ? 0.0 : (num/den);

    auto star_flux = [&](const Prim& W, const State& U, const State& F, double S) noexcept -> State {
        const double denom = S - Sstar;
        if (std::abs(denom) < 1e-12) return hlle_flux(WL, WR, g); // safeguard
        const double rho   = W.rho, u = W.u, p = W.p;
        const double pstar = p + rho*(S - u)*(Sstar - u);
        const double rho_star = rho*(S - u)/denom;
        if (!(pstar>0.0) || !(rho_star>0.0)) return hlle_flux(WL, WR, g); // positivity fallback
        const State Ustar{rho_star, rho_star*Sstar,
                          ((S - u)*U.E - p*u + pstar*Sstar)/denom};
        return State{ F.rho + S*(Ustar.rho - U.rho),
                      F.mom + S*(Ustar.mom - U.mom),
                      F.E   + S*(Ustar.E   - U.E) };
    };

    const State FLs = star_flux(WL, UL, FL, SL);
    const State FRs = star_flux(WR, UR, FR, SR);
    return (Sstar >= 0.0) ? FLs : FRs;
}

// ---- Exact Sod solution for comparison ----
struct StarState { double p, u; };

static inline void pressure_function(double p, double rho, double p0, double c,
                                     double g, double& f, double& df) noexcept {
    if (p > p0) {
        const double A = 2.0 / ((g + 1.0) * rho);
        const double B = (g - 1.0) / (g + 1.0) * p0;
        f  = (p - p0) * std::sqrt(A / (p + B));
        df = std::sqrt(A / (p + B)) * (1.0 - 0.5 * (p - p0) / (p + B));
    } else {
        const double pr = p / p0;
        f  = 2.0 * c / (g - 1.0) * (std::pow(pr, (g - 1.0) / (2.0 * g)) - 1.0);
        df = (1.0 / (rho * c)) * std::pow(pr, -(g + 1.0) / (2.0 * g));
    }
}

static inline StarState star_region(double rhoL, double uL, double pL,
                                    double rhoR, double uR, double pR,
                                    double g) noexcept {
    const double cL = std::sqrt(g * pL / rhoL);
    const double cR = std::sqrt(g * pR / rhoR);
    double p = std::max(1e-6, 0.5 * (pL + pR));
    double fL, dfL, fR, dfR;
    for (int k = 0; k < 100; ++k) {
        pressure_function(p, rhoL, pL, cL, g, fL, dfL);
        pressure_function(p, rhoR, pR, cR, g, fR, dfR);
        const double num = fL + fR + (uR - uL);
        const double den = dfL + dfR;
        const double p_new = p - num / den;
        if (std::abs(2.0 * (p_new - p) / (p_new + p)) < 1e-8) { p = p_new; break; }
        p = p_new;
    }
    const double u = 0.5 * (uL + uR + fR - fL);
    return StarState{p, u};
}

static inline double density_star(double p, double rho, double p0, double g) noexcept {
    if (p > p0) {
        return rho * ((p / p0 + (g - 1.0) / (g + 1.0)) /
                      ((g - 1.0) / (g + 1.0) * p / p0 + 1.0));
    } else {
        return rho * std::pow(p / p0, 1.0 / g);
    }
}

static inline Prim sample(double xi, const StarState& s,
                          double rhoL, double uL, double pL,
                          double rhoR, double uR, double pR,
                          double g) noexcept {
    const double cL = std::sqrt(g * pL / rhoL);
    const double cR = std::sqrt(g * pR / rhoR);
    const double p_star = s.p;
    const double u_star = s.u;
    if (xi < u_star) {
        if (p_star > pL) {
            const double SL = uL - cL * std::sqrt((g + 1.0)/(2.0*g) * p_star/pL + (g - 1.0)/(2.0*g));
            if (xi < SL) return Prim{rhoL, uL, pL};
            const double rho = density_star(p_star, rhoL, pL, g);
            return Prim{rho, u_star, p_star};
        } else {
            const double SHL = uL - cL;
            const double cL_star = cL * std::pow(p_star/pL, (g - 1.0)/(2.0*g));
            const double STL = u_star - cL_star;
            if (xi < SHL) {
                return Prim{rhoL, uL, pL};
            } else if (xi < STL) {
                const double u = (2.0/(g + 1.0))*(cL + 0.5*(g - 1.0)*uL + xi);
                const double c = (2.0/(g + 1.0))*(cL + 0.5*(g - 1.0)*(uL - xi));
                const double rho = rhoL * std::pow(c/cL, 2.0/(g - 1.0));
                const double p = pL * std::pow(c/cL, 2.0*g/(g - 1.0));
                return Prim{rho, u, p};
            } else {
                const double rho = density_star(p_star, rhoL, pL, g);
                return Prim{rho, u_star, p_star};
            }
        }
    } else {
        if (p_star > pR) {
            const double SR = uR + cR * std::sqrt((g + 1.0)/(2.0*g) * p_star/pR + (g - 1.0)/(2.0*g));
            if (xi > SR) return Prim{rhoR, uR, pR};
            const double rho = density_star(p_star, rhoR, pR, g);
            return Prim{rho, u_star, p_star};
        } else {
            const double SHR = uR + cR;
            const double cR_star = cR * std::pow(p_star/pR, (g - 1.0)/(2.0*g));
            const double STR = u_star + cR_star;
            if (xi > SHR) {
                return Prim{rhoR, uR, pR};
            } else if (xi > STR) {
                const double u = (2.0/(g + 1.0))*(-cR + 0.5*(g - 1.0)*uR + xi);
                const double c = (2.0/(g + 1.0))*(cR - 0.5*(g - 1.0)*(xi - uR));
                const double rho = rhoR * std::pow(c/cR, 2.0/(g - 1.0));
                const double p = pR * std::pow(c/cR, 2.0*g/(g - 1.0));
                return Prim{rho, u, p};
            } else {
                const double rho = density_star(p_star, rhoR, pR, g);
                return Prim{rho, u_star, p_star};
            }
        }
    }
}

int main(){
    std::ios::sync_with_stdio(false); std::cin.tie(nullptr);
    Params P; const int ng=2; const int Ntot=P.N+2*ng; const double L=P.xR-P.xL; const double dx=L/P.N;

    std::vector<double> x(Ntot); for(int i=0;i<Ntot;++i) x[i]=P.xL+(i-ng+0.5)*dx;
    std::vector<State> U(Ntot), U0(Ntot), RHS(Ntot); std::vector<Prim> W(Ntot);
    const int Nf=P.N+1; std::vector<Prim> WLf(Nf), WRf(Nf); std::vector<State> F(Nf);

    auto to_prim = [&](){ for(int i=0;i<Ntot;++i) W[i]=cons2prim(U[i],P.gamma,P.p_floor,P.rho_floor); };

    auto init = [&](){
        for(int i=0;i<Ntot;++i){ const bool left = (x[i] < 0.5*(P.xL+P.xR)); Prim w = left? Prim{1.0,0.0,1.0} : Prim{0.125,0.0,0.1}; U[i]=prim2cons(w,P.gamma);} };

    auto apply_bc = [&](){ for(int k=0;k<ng;++k){ U[k]=U[ng]; U[Ntot-1-k]=U[Ntot-1-ng]; } };

    // MUSCL: face i+1/2 gets left= i^R, right= (i+1)^L
    auto reconstruct_interfaces = [&](){
        // slopes per cell
        static std::vector<Prim> S; S.resize(Ntot);
        for(int i=1;i<Ntot-1;++i){ // safe due to ghosts
            S[i].rho = minmod(W[i].rho - W[i-1].rho, W[i+1].rho - W[i].rho);
            S[i].u   = minmod(W[i].u   - W[i-1].u,   W[i+1].u   - W[i].u  );
            S[i].p   = minmod(W[i].p   - W[i-1].p,   W[i+1].p   - W[i].p  );
        }
        for(int i=ng;i<ng+P.N;++i){ int f=i-ng; // face between i and i+1
            const Prim WLc{ std::max(W[i].rho + 0.5*S[i].rho, P.rho_floor), W[i].u + 0.5*S[i].u, std::max(W[i].p + 0.5*S[i].p, P.p_floor) };
            const Prim WRc{ std::max(W[i+1].rho - 0.5*S[i+1].rho, P.rho_floor), W[i+1].u - 0.5*S[i+1].u, std::max(W[i+1].p - 0.5*S[i+1].p, P.p_floor) };
            WLf[f]=WLc; WRf[f]=WRc;
        }
        // also construct the last face at ng+P.N (right boundary)
        {
            int i = ng+P.N-1; int f = P.N; // last interior cell i and face f=P.N between i and i+1
            const Prim WLc{ std::max(W[i].rho + 0.5*S[i].rho, P.rho_floor), W[i].u + 0.5*S[i].u, std::max(W[i].p + 0.5*S[i].p, P.p_floor) };
            const Prim WRc{ std::max(W[i+1].rho - 0.5*S[i+1].rho, P.rho_floor), W[i+1].u - 0.5*S[i+1].u, std::max(W[i+1].p - 0.5*S[i+1].p, P.p_floor) };
            WLf[f]=WLc; WRf[f]=WRc;
        }
    };

    auto compute_rhs = [&](){
        apply_bc(); to_prim();
        if(P.use_muscl) reconstruct_interfaces(); else { for(int f=0; f<Nf; ++f){ WLf[f]=W[ng+f-1]; WRf[f]=W[ng+f]; } }
        for(int f=0; f<Nf; ++f) F[f]=hllc_flux(WLf[f],WRf[f],P.gamma);
        std::fill(RHS.begin(), RHS.end(), State{0,0,0});
        for(int i=ng;i<ng+P.N;++i){ int f=i-ng; RHS[i].rho=-(F[f+1].rho-F[f].rho)/dx; RHS[i].mom=-(F[f+1].mom-F[f].mom)/dx; RHS[i].E=-(F[f+1].E-F[f].E)/dx; }
    };

    auto dump = [&](const std::string& fname){ std::ofstream ofs(fname); ofs.setf(std::ios::scientific); ofs<<std::setprecision(15)<<"x,rho,u,p,E\n"; to_prim(); for(int i=ng;i<ng+P.N;++i) ofs<<x[i]<<","<<W[i].rho<<","<<W[i].u<<","<<W[i].p<<","<<U[i].E<<"\n"; };

    init(); apply_bc();
    double t=0.0; int step=0;
    auto max_wave = [&](){ to_prim(); double amax=0.0; for(int i=ng;i<ng+P.N;++i) amax=std::max(amax, std::abs(W[i].u)+sound_speed(W[i],P.gamma)); return amax; };

    while(t < P.t_end){ double amax=max_wave(); if(amax<=0.0) break; double dt=P.CFL*dx/amax; if(t+dt>P.t_end) dt=P.t_end-t;
        U0=U; compute_rhs(); for(int i=ng;i<ng+P.N;++i){ U[i].rho=U0[i].rho+dt*RHS[i].rho; U[i].mom=U0[i].mom+dt*RHS[i].mom; U[i].E=U0[i].E+dt*RHS[i].E; }
        compute_rhs(); for(int i=ng;i<ng+P.N;++i){ U[i].rho=0.5*(U0[i].rho+U[i].rho+dt*RHS[i].rho); U[i].mom=0.5*(U0[i].mom+U[i].mom+dt*RHS[i].mom); U[i].E=0.5*(U0[i].E+U[i].E+dt*RHS[i].E); }
        // simple positivity fix (limited): floor rho,p and re-sync E if needed
        to_prim();
        for(int i=ng;i<ng+P.N;++i){ if(W[i].rho<P.rho_floor || W[i].p<P.p_floor){ W[i].rho=std::max(W[i].rho,P.rho_floor); W[i].p=std::max(W[i].p,P.p_floor); U[i]=prim2cons(W[i],P.gamma);} }
        t+=dt; ++step; if(P.output_every_step){ char buf[64]; std::snprintf(buf,sizeof(buf),"solution_%06d.csv",step); dump(buf);} }
    dump("solution.csv");

    // compute exact Sod solution and L1 error in density
    std::vector<Prim> exact;
    auto write_exact = [&](const std::string& fname){
        const double rhoL=1.0, uL=0.0, pL=1.0;
        const double rhoR=0.125, uR=0.0, pR=0.1;
        StarState st = star_region(rhoL,uL,pL,rhoR,uR,pR,P.gamma);
        std::ofstream ofs(fname);
        ofs.setf(std::ios::scientific);
        ofs<<std::setprecision(15)<<"x,rho,u,p\n";
        exact.resize(P.N);
        for(int i=0;i<P.N;++i){
            double xi=(x[i+ng]-0.5*(P.xL+P.xR))/P.t_end;
            Prim s=sample(xi,st,rhoL,uL,pL,rhoR,uR,pR,P.gamma);
            exact[i]=s;
            ofs<<x[i+ng]<<","<<s.rho<<","<<s.u<<","<<s.p<<"\n";
        }
    };
    write_exact("exact.csv");

    to_prim();
    double l1=0.0;
    for(int i=ng;i<ng+P.N;++i) l1+=std::abs(W[i].rho-exact[i-ng].rho);
    l1/=P.N;

    std::cerr.setf(std::ios::scientific);
    std::cerr<<std::setprecision(6)<<"Finished t="<<t<<", steps="<<step<<", L1_rho="<<l1<<"\n";

    return 0;
}

