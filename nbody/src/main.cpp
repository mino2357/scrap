#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <omp.h>

#include "util.hpp"
#include "integrators.hpp"
#include "energy.hpp"

// --- Simple live plotting via gnuplot (optional) ---
#include <cstdio>
struct Plotter {
  FILE* gp = nullptr;
  bool enabled = false;
  int every = 50;
  std::size_t limit = 4096;
  int axes = 0; // 0:xy,1:xz,2:yz
  void open(){
    if (enabled && !gp){
      gp = popen("gnuplot -persist", "w");
      if (gp){
        fprintf(gp, "unset key\nset size ratio -1\nset xlabel 'x'\nset ylabel 'y'\n");
        fflush(gp);
      }
    }
  }
  void close(){ if (gp){ fflush(gp); pclose(gp); gp=nullptr; } }
  void frame(const Arrays& S, int step, double t){
    if (!gp) return;
    const double* a = (axes==0)? S.x : (axes==1? S.x : S.y);
    const double* b = (axes==0)? S.y : (axes==1? S.z : S.z);
    std::size_t M = std::min<std::size_t>(S.N, limit);
    double amin=1e300, amax=-1e300, bmin=1e300, bmax=-1e300;
    for (std::size_t i=0;i<M;++i){
      double ai = a[i], bi = b[i];
      if (ai<amin) amin=ai; if (ai>amax) amax=ai;
      if (bi<bmin) bmin=bi; if (bi>bmax) bmax=bi;
    }
    double da = (amax-amin); double db = (bmax-bmin);
    if (da<=0) da = 1.0; if (db<=0) db = 1.0;
    double pad=0.05;
    double ax0=amin-pad*da, ax1=amax+pad*da;
    double by0=bmin-pad*db, by1=bmax+pad*db;
    fprintf(gp, "set title 'N-body (step=%d, t=%.6g)'\n", step, t);
    fprintf(gp, "set xrange [%.15g:%.15g]\n", ax0, ax1);
    fprintf(gp, "set yrange [%.15g:%.15g]\n", by0, by1);
    fprintf(gp, "plot '-' with points pt 7 ps 0.3 notitle\n");
    for (std::size_t i=0;i<M;++i) fprintf(gp, "%.15g %.15g\n", a[i], b[i]);
    fprintf(gp, "e\n"); fflush(gp);
  }
};

extern "C" bool* __get_kernel_flag(); // from integrators.cpp

static void usage(){
  std::puts("Usage: ./nbody --N <int> --steps <int> --dt <float> --eps <float> "
            "--method {yoshida4|rk4|verlet} --Bj <int> --threads <int> "
            "[--seed <uint64>] [--check 0|1] [--selftest] "
            "[--plot live|none] [--plot_every k] [--plot_limit m] [--plot_axes xy|xz|yz] "
            "[--kernel fused|two_pass] [--autosweep CSV]");
}

static void init_random(Arrays& S, uint64_t seed){
  RNG rng(seed);
  for (std::size_t i=0;i<S.N;++i){
    S.x[i]=rng.uniform(-1.0,1.0);
    S.y[i]=rng.uniform(-1.0,1.0);
    S.z[i]=rng.uniform(-1.0,1.0);
    S.vx[i]=rng.uniform(-0.1,0.1);
    S.vy[i]=rng.uniform(-0.1,0.1);
    S.vz[i]=rng.uniform(-0.1,0.1);
    S.m[i]=1.0/S.N;
  }
}

static void init_two_body(Arrays& S){
  if (S.N < 2) return;
  for (std::size_t i=0;i<S.N;++i){
    S.x[i]=S.y[i]=S.z[i]=S.vx[i]=S.vy[i]=S.vz[i]=0.0; S.m[i]=0.0;
  }
  double m1=1.0, m2=3e-6;
  S.m[0]=m1; S.m[1]=m2;
  S.x[0]=0.0; S.y[0]=0.0; S.z[0]=0.0;
  S.x[1]=1.0; S.y[1]=0.0; S.z[1]=0.0;
  double M=m1+m2;
  double v = std::sqrt(M/1.0);
  S.vx[0]=0.0; S.vy[0]= v*m2/M; S.vz[0]=0.0;
  S.vx[1]=0.0; S.vy[1]=-v*m1/M; S.vz[1]=0.0;
}

int main(int argc, char** argv){
  enable_FTZ_DAZ();

  std::size_t N=4096;
  int steps=100;
  double dt=1e-3;
  double eps=1e-3;
  std::size_t Bj=512;
  std::string method="yoshida4";
  int threads = std::max(1, omp_get_max_threads());
  int check=1;
  uint64_t seed=42;
  bool selftest=false;
  std::string plot_mode="none";
  int plot_every=50;
  std::size_t plot_limit=4096;
  std::string plot_axes="xy";
  std::string kernel="fused";
  std::string autosweep="";

  for (int i=1;i<argc;i++){
    auto eq = [&](const char* k){ return std::strcmp(argv[i],k)==0; };
    if      (eq("--N") && i+1<argc) N = std::strtoull(argv[++i],nullptr,10);
    else if (eq("--steps") && i+1<argc) steps = std::atoi(argv[++i]);
    else if (eq("--dt") && i+1<argc) dt = std::atof(argv[++i]);
    else if (eq("--eps") && i+1<argc) eps = std::atof(argv[++i]);
    else if (eq("--Bj") && i+1<argc) Bj = std::strtoull(argv[++i],nullptr,10);
    else if (eq("--method") && i+1<argc) method = argv[++i];
    else if (eq("--threads") && i+1<argc) threads = std::atoi(argv[++i]);
    else if (eq("--check") && i+1<argc) check = std::atoi(argv[++i]);
    else if (eq("--seed") && i+1<argc) seed = std::strtoull(argv[++i],nullptr,10);
    else if (eq("--selftest")) selftest=true;
    else if (eq("--plot") && i+1<argc) plot_mode = argv[++i];
    else if (eq("--plot_every") && i+1<argc) plot_every = std::atoi(argv[++i]);
    else if (eq("--plot_limit") && i+1<argc) plot_limit = std::strtoull(argv[++i],nullptr,10);
    else if (eq("--plot_axes") && i+1<argc) plot_axes = argv[++i];
    else if (eq("--kernel") && i+1<argc) kernel = argv[++i];
    else if (eq("--autosweep") && i+1<argc) autosweep = argv[++i];
    else { usage(); return 0; }
  }

  omp_set_num_threads(threads);

  Arrays S;
  S.N=N;
  S.x  = static_cast<double*>(aligned_malloc(sizeof(double)*N));
  S.y  = static_cast<double*>(aligned_malloc(sizeof(double)*N));
  S.z  = static_cast<double*>(aligned_malloc(sizeof(double)*N));
  S.vx = static_cast<double*>(aligned_malloc(sizeof(double)*N));
  S.vy = static_cast<double*>(aligned_malloc(sizeof(double)*N));
  S.vz = static_cast<double*>(aligned_malloc(sizeof(double)*N));
  S.m  = static_cast<double*>(aligned_malloc(sizeof(double)*N));
  S.ax = static_cast<double*>(aligned_malloc(sizeof(double)*N));
  S.ay = static_cast<double*>(aligned_malloc(sizeof(double)*N));
  S.az = static_cast<double*>(aligned_malloc(sizeof(double)*N));

  if (!S.x||!S.y||!S.z||!S.vx||!S.vy||!S.vz||!S.m||!S.ax||!S.ay||!S.az){
    std::fprintf(stderr,"Allocation failed\n");
    return 1;
  }

  // Plotter
  Plotter plt;
  if (plot_mode=="live"){
    plt.enabled = true;
    plt.every = std::max(1, plot_every);
    plt.limit = plot_limit;
    plt.axes = (plot_axes=="xy"?0:(plot_axes=="xz"?1:2));
    plt.open();
  }

  // Kernel mode flag
  *(__get_kernel_flag()) = (kernel!="two_pass");
  std::printf("[INFO] kernel=%s  Bj=%zu  threads=%d\n", kernel.c_str(), Bj, threads);

  // Helper for quick bench (used by autosweep)
  auto run_once = [&](int steps_bench)->double{
    Timer t; t.tic();
    for (int s=0;s<steps_bench;++s){
      if (method=="yoshida4") step_yoshida4(dt, S, eps*eps, Bj);
      else if (method=="verlet") step_verlet(dt, S, eps*eps, Bj);
      else step_yoshida4(dt, S, eps*eps, Bj);
    }
    double sec = t.toc();
    double pairs = double(S.N) * double(S.N) * steps_bench;
    return pairs/sec;
  };

  // Selftest or autosweep or normal run
  if (!autosweep.empty()){
    // parse CSV list
    std::vector<int> list; int cur=0; bool in=false; bool neg=false;
    for (char c: autosweep){
      if (c=='-' && !in){ in=true; neg=true; cur=0; }
      else if (c>='0' && c<='9'){ in=true; cur = cur*10 + (c-'0'); }
      else if (c==',' && in){ list.push_back(neg?-cur:cur); in=false; neg=false; cur=0; }
    }
    if (in) list.push_back(neg?-cur:cur);
    if (list.empty()){ std::fprintf(stderr,"Bad --autosweep list\n"); return 2; }

    init_random(S, seed);
    std::printf("[SWEEP] N=%zu steps=50 dt=%.3g eps=%.3g kernel=%s\n",
                S.N, dt, eps, kernel.c_str());
    double best_pps=0.0; int best_bj=(int)Bj;
    for (int b: list){
      if (b<=0) continue;
      Bj = (std::size_t)b;
      for (std::size_t i=0;i<S.N;++i){ S.vx[i]=S.vy[i]=S.vz[i]=0.0; }
      double pps = run_once(50);
      std::printf("[SWEEP] Bj=%d  Pairs/s=%.3e\n", b, pps);
      if (pps>best_pps){ best_pps=pps; best_bj=b; }
    }
    std::printf("[SWEEP] BEST Bj=%d  Pairs/s=%.3e\n", best_bj, best_pps);
    plt.close();
    aligned_free(S.x); aligned_free(S.y); aligned_free(S.z);
    aligned_free(S.vx); aligned_free(S.vy); aligned_free(S.vz);
    aligned_free(S.m); aligned_free(S.ax); aligned_free(S.ay); aligned_free(S.az);
    return 0;
  }

  if (selftest){
    S.N = N = 2;
    init_two_body(S);
    double eps2 = eps*eps;
    double H0 = kinetic_energy(S.vx,S.vy,S.vz,S.m,S.N) + potential_energy(S.x,S.y,S.z,S.m,S.N,eps2);
    double L0[3]; angular_momentum(S.x,S.y,S.z,S.vx,S.vy,S.vz,S.m,S.N,L0);

    Timer tim; tim.tic();
    for (int s=0;s<steps;++s){
      if (plt.enabled && (s%plt.every==0)) plt.frame(S, s, s*dt);
      if (method=="yoshida4") step_yoshida4(dt, S, eps2, Bj);
      else if (method=="rk4") step_rk4(dt, S, eps2, Bj);
      else if (method=="verlet") step_verlet(dt, S, eps2, Bj);
      else { std::fprintf(stderr,"Unknown method\n"); return 3; }
    }
    double sec = tim.toc();
    double H1 = kinetic_energy(S.vx,S.vy,S.vz,S.m,S.N) + potential_energy(S.x,S.y,S.z,S.m,S.N,eps2);
    double L1[3]; angular_momentum(S.x,S.y,S.z,S.vx,S.vy,S.vz,S.m,S.N,L1);
    double dH = std::abs(H1-H0)/std::max(1.0,std::abs(H0));
    double dL = std::sqrt((L1[0]-L0[0])*(L1[0]-L0[0])+(L1[1]-L0[1])*(L1[1]-L0[1])+(L1[2]-L0[2])*(L1[2]-L0[2]));
    double L0n= std::sqrt(L0[0]*L0[0]+L0[1]*L0[1]+L0[2]*L0[2]);
    double dLrel = dL/std::max(1.0,L0n);

    std::printf("[SELFTEST] method=%s steps=%d dt=%.6g eps=%.1e  time=%.3fs\n",
                method.c_str(), steps, dt, eps, sec);
    std::printf("[SELFTEST] rel|ΔH|=%.3e  rel|ΔL|=%.3e\n", dH, dLrel);
    plt.close();
  } else {
    init_random(S, seed);
    double eps2 = eps*eps;
    const double pairs_per_step = double(S.N) * double(S.N);
    Timer tim; tim.tic();
    for (int s=0;s<steps;++s){
      if (plt.enabled && (s%plt.every==0)) plt.frame(S, s, s*dt);
      if (method=="yoshida4") step_yoshida4(dt, S, eps2, Bj);
      else if (method=="rk4") step_rk4(dt, S, eps2, Bj);
      else if (method=="verlet") step_verlet(dt, S, eps2, Bj);
      else { std::fprintf(stderr,"Unknown method\n"); return 3; }
      if (check && (s%50==0 || s==steps-1)){
        double H = kinetic_energy(S.vx,S.vy,S.vz,S.m,S.N) + potential_energy(S.x,S.y,S.z,S.m,S.N,eps2);
        std::printf("step %d  H=%.12e\n", s+1, H);
      }
    }
    double sec = tim.toc();
    double pairs = pairs_per_step * steps;
    double pairs_per_sec = pairs / sec;
    double gflops = (pairs_per_sec * 32.0) / 1e9;
    std::printf("Done: N=%zu steps=%d dt=%.3g eps=%.3g method=%s\n", S.N, steps, dt, eps, method.c_str());
    std::printf("Time = %.3fs  Pairs/s = %.3e  ~GFLOP/s ≈ %.2f (rough)\n", sec, pairs_per_sec, gflops);
    plt.close();
  }

  aligned_free(S.x); aligned_free(S.y); aligned_free(S.z);
  aligned_free(S.vx); aligned_free(S.vy); aligned_free(S.vz);
  aligned_free(S.m); aligned_free(S.ax); aligned_free(S.ay); aligned_free(S.az);
  return 0;
}
