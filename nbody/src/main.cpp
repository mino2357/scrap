//
// N-body ミニアプリ（AVX2 + OpenMP）
// ----------------------------------------
// - SoA 配列に粒子状態を格納し、力計算は AVX2/FMA ベクトル化カーネルを使用。
// - 時間積分は Yoshida4 / RK4 / Verlet を選択可能。
// - オプションでライブプロット（gnuplot）や Bj のパラメータスイープが可能。
// - コマンドライン引数で N, steps, dt, eps, Bj, threads 等を指定。
//   例: ./nbody --N 8192 --steps 1000 --dt 1e-3 --eps 1e-3 --method yoshida4 --Bj 512
//
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

// --- Two-body helpers: elliptical ICs and orbital diagnostics ---
static inline void rel_vec2(const double* x,const double* y,const double* z,
                            const double* vx,const double* vy,const double* vz,
                            double r[3], double v[3]){
  // r = r2 - r1, v = v2 - v1 (assumes N>=2)
  r[0] = x[1]-x[0]; r[1] = y[1]-y[0]; r[2] = z[1]-z[0];
  v[0] = vx[1]-vx[0]; v[1] = vy[1]-vy[0]; v[2] = vz[1]-vz[0];
}
static inline void cross3(const double a[3], const double b[3], double c[3]){
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}
static inline double dot3(const double a[3], const double b[3]){
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
static inline double norm3(const double a[3]){
  return std::sqrt(dot3(a,a));
}

static void init_two_body_ellipse(Arrays& S, double m1, double m2,
                                  double a, double e, const std::string& plane){
  // Place two bodies on an ellipse at pericenter with COM at origin
  // r_p = a(1-e), v_p = sqrt( mu*(1+e)/(a*(1-e)) ), with mu = m1+m2
  if (S.N < 2) return;
  const double mu = m1 + m2;
  const double rp = a * (1.0 - e);
  const double vp = std::sqrt(mu * (1.0 + e) / (a * (1.0 - e)));

  // Relative state at pericenter: r along +x, v along +y (counter-clockwise)
  double rrel[3] = {rp, 0.0, 0.0};
  double vrel[3] = {0.0, vp, 0.0};
  // Rotate into requested plane if needed
  // plane: xy (default), xz, yz. For xz, map y->z; for yz, map x->y, y->z
  if (plane == "xz"){
    rrel[2] = rrel[1]; rrel[1] = 0.0;
    vrel[2] = vrel[1]; vrel[1] = 0.0;
  } else if (plane == "yz"){
    rrel[0] = 0.0; rrel[1] = rp; rrel[2] = 0.0;
    vrel[0] = -vp; vrel[1] = 0.0; vrel[2] = 0.0;
  }

  // COM split: r1 = -m2/M * rrel, r2 = +m1/M * rrel; same for velocities
  double M = mu;
  double c1 = -m2 / M;
  double c2 =  m1 / M;

  for (std::size_t i=0;i<S.N;++i){
    S.x[i]=S.y[i]=S.z[i]=S.vx[i]=S.vy[i]=S.vz[i]=0.0; S.m[i]=0.0;
  }
  S.m[0] = m1; S.m[1] = m2;
  S.x[0] = c1 * rrel[0]; S.y[0] = c1 * rrel[1]; S.z[0] = c1 * rrel[2];
  S.x[1] = c2 * rrel[0]; S.y[1] = c2 * rrel[1]; S.z[1] = c2 * rrel[2];
  S.vx[0]= c1 * vrel[0]; S.vy[0]= c1 * vrel[1]; S.vz[0]= c1 * vrel[2];
  S.vx[1]= c2 * vrel[0]; S.vy[1]= c2 * vrel[1]; S.vz[1]= c2 * vrel[2];
}

static void orbital_elements_from_state(const Arrays& S, double m1, double m2,
                                        double& a, double& e){
  double r[3], v[3];
  rel_vec2(S.x,S.y,S.z, S.vx,S.vy,S.vz, r, v);
  const double rnorm = norm3(r);
  const double mu = m1 + m2; // G=1
  // h = r x v
  double h[3]; cross3(r,v,h);
  // evec = (v x h)/mu - r/|r|
  double vxh[3]; cross3(v,h,vxh);
  double evec[3] = { vxh[0]/mu - r[0]/rnorm,
                     vxh[1]/mu - r[1]/rnorm,
                     vxh[2]/mu - r[2]/rnorm };
  e = norm3(evec);
  // specific orbital energy
  const double v2 = dot3(v,v);
  const double eps = 0.5*v2 - mu/rnorm;
  a = - mu / (2.0*eps);
}

// --- Simple live plotting via gnuplot (optional) ---
// 軽量なプロッタ。粒子数を上限で間引き、ステップごとに座標を散布図で描画。
// 注意: popen/gnuplot 依存のため、環境によっては使用不可の場合があります。
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
    // TODO: 自動スケーリングを使うなら ax0/ax1, by0/by1 を xrange/yrange に反映
    fprintf(gp, "set xrange [-2:2]\n");
    fprintf(gp, "set yrange [-2:2]\n");
    fprintf(gp, "plot '-' with points pt 7 ps 0.1 notitle\n");
    for (std::size_t i=0;i<M;++i) fprintf(gp, "%.15g %.15g\n", a[i], b[i]);
    fprintf(gp, "e\n"); fflush(gp);
  }
};

extern "C" bool* __get_kernel_flag(); // from integrators.cpp
extern "C" int* __get_vector_mode();  // from force_kernel_avx2.cpp
extern "C" int* __get_scalar_tiled(); // from force_kernel_avx2.cpp

static void usage(){
  std::puts("Usage: ./nbody --N <int> --steps <int> --dt <float> --eps <float> "
            "--method {yoshida4|rk4|verlet} --Bj <int> --threads <int> "
            "[--seed <uint64>] [--check 0|1] [--selftest] [--two_body] "
            "[--m1 <float>] [--m2 <float>] [--a <float>] [--e <float>] [--plane xy|xz|yz] "
            "[--out_csv <path>] [--check_every <int>] "
            "[--plot live|none] [--plot_every k] [--plot_limit m] [--plot_axes xy|xz|yz] "
            "[--kernel fused|two_pass] [--mode vector|scalar|scalar_tiled] [--autosweep CSV]");
}

static void init_random(Arrays& S, uint64_t seed){
  RNG rng(seed);
  for (std::size_t i=0;i<S.N;++i){
    S.x[i]=rng.uniform(-1.0,1.0);
    S.y[i]=rng.uniform(-1.0,1.0);
    S.z[i]=rng.uniform(-1.0,1.0);
    // 初期速度は 0。ランダムにしたい場合は下のコメントアウトを有効化
    S.vx[i]=0.0; //rng.uniform(-0.1,0.1);
    S.vy[i]=0.0; //rng.uniform(-0.1,0.1);
    S.vz[i]=0.0; //rng.uniform(-0.1,0.1);
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
  bool two_body=false;
  double tb_m1=1.0, tb_m2=3e-6, tb_a=1.0, tb_e=0.0;
  std::string tb_plane="xy";
  std::string out_csv="";
  int check_every=10;
  double rand_vel=0.0; // if >0, initialize random velocities in [-rand_vel,rand_vel]
  std::string plot_mode="none";
  int plot_every=50;
  std::size_t plot_limit=4096;
  std::string plot_axes="xy";
  std::string kernel="fused";
  std::string mode="vector"; // vector | scalar | scalar_tiled
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
    else if (eq("--two_body")) two_body=true;
    else if (eq("--m1") && i+1<argc) tb_m1 = std::atof(argv[++i]);
    else if (eq("--m2") && i+1<argc) tb_m2 = std::atof(argv[++i]);
    else if (eq("--a") && i+1<argc) tb_a = std::atof(argv[++i]);
    else if (eq("--e") && i+1<argc) tb_e = std::atof(argv[++i]);
    else if (eq("--plane") && i+1<argc) tb_plane = argv[++i];
    else if (eq("--out_csv") && i+1<argc) out_csv = argv[++i];
    else if (eq("--check_every") && i+1<argc) check_every = std::max(1, std::atoi(argv[++i]));
    else if (eq("--rand_vel") && i+1<argc) rand_vel = std::atof(argv[++i]);
    else if (eq("--plot") && i+1<argc) plot_mode = argv[++i];
    else if (eq("--plot_every") && i+1<argc) plot_every = std::atoi(argv[++i]);
    else if (eq("--plot_limit") && i+1<argc) plot_limit = std::strtoull(argv[++i],nullptr,10);
    else if (eq("--plot_axes") && i+1<argc) plot_axes = argv[++i];
    else if (eq("--kernel") && i+1<argc) kernel = argv[++i];
    else if (eq("--mode") && i+1<argc) mode = argv[++i];
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
  // Vectorization/scalar mode flags
  if (mode=="vector"){
    *(__get_vector_mode()) = 1;
    *(__get_scalar_tiled()) = 0;
  } else if (mode=="scalar_tiled"){
    *(__get_vector_mode()) = 0;
    *(__get_scalar_tiled()) = 1;
  } else { // scalar
    *(__get_vector_mode()) = 0;
    *(__get_scalar_tiled()) = 0;
  }
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
    // parse CSV list （例: "128,256,512,1024"）→ Bj 候補として順にベンチ
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
    // 2 体問題の保存量（エネルギー・角運動量）を確認
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
  } else if (two_body) {
    // Two-body elliptical test with diagnostics
    S.N = N = 2;
    init_two_body_ellipse(S, tb_m1, tb_m2, tb_a, tb_e, tb_plane);
    const double eps2 = eps*eps;

    // Initial diagnostics
    double a0,e0; orbital_elements_from_state(S, tb_m1, tb_m2, a0, e0);
    double H0 = kinetic_energy(S.vx,S.vy,S.vz,S.m,S.N) + potential_energy(S.x,S.y,S.z,S.m,S.N,eps2);
    double L0[3]; angular_momentum(S.x,S.y,S.z,S.vx,S.vy,S.vz,S.m,S.N,L0);
    const double L0n = std::sqrt(L0[0]*L0[0]+L0[1]*L0[1]+L0[2]*L0[2]);

    // Optional CSV
    FILE* fp = nullptr;
    if (!out_csv.empty()){
      fp = std::fopen(out_csv.c_str(), "w");
      if (fp) std::fprintf(fp, "step,t,r,e,a,rel_dE,rel_dL\n");
    }

    std::printf("[2BODY] m1=%.6g m2=%.6g  a=%.6g e=%.6g  plane=%s\n",
                tb_m1,tb_m2,a0,e0,tb_plane.c_str());
    Timer tim; tim.tic();
    for (int s=0;s<steps;++s){
      if (plt.enabled && (s%plt.every==0)) plt.frame(S, s, s*dt);
      if (method=="yoshida4") step_yoshida4(dt, S, eps2, Bj);
      else if (method=="rk4") step_rk4(dt, S, eps2, Bj);
      else if (method=="verlet") step_verlet(dt, S, eps2, Bj);
      else { std::fprintf(stderr,"Unknown method\n"); if (fp) std::fclose(fp); return 3; }

      if ((s%check_every)==0 || s==steps-1){
        double a1,e1; orbital_elements_from_state(S, tb_m1, tb_m2, a1, e1);
        double rvec[3], vvec[3]; rel_vec2(S.x,S.y,S.z, S.vx,S.vy,S.vz, rvec, vvec);
        double rnorm = std::sqrt(rvec[0]*rvec[0]+rvec[1]*rvec[1]+rvec[2]*rvec[2]);
        double H1 = kinetic_energy(S.vx,S.vy,S.vz,S.m,S.N) + potential_energy(S.x,S.y,S.z,S.m,S.N,eps2);
        double L1[3]; angular_momentum(S.x,S.y,S.z,S.vx,S.vy,S.vz,S.m,S.N,L1);
        double dH = std::abs(H1-H0)/std::max(1.0,std::abs(H0));
        double dL = std::sqrt((L1[0]-L0[0])*(L1[0]-L0[0])+(L1[1]-L0[1])*(L1[1]-L0[1])+(L1[2]-L0[2])*(L1[2]-L0[2]));
        double dLrel = dL/std::max(1.0,L0n);
        double t = (s+1)*dt;
        std::printf("step %7d t=%.6g  e=%.6g  a=%.6g  rel|ΔH|=%.3e  rel|ΔL|=%.3e\n",
                    s+1, t, e1, a1, dH, dLrel);
        if (fp){ std::fprintf(fp, "%d,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g\n", s+1, t, rnorm, e1, a1, dH, dLrel); }
      }
    }
    double sec = tim.toc();
    std::printf("[2BODY] Done. steps=%d dt=%.3g eps=%.3g method=%s  time=%.3fs\n",
                steps, dt, eps, method.c_str(), sec);
    if (fp) { std::fclose(fp); std::printf("[2BODY] CSV saved: %s\n", out_csv.c_str()); }
    plt.close();
  } else {
    init_random(S, seed);
    if (rand_vel>0){
      RNG rng(seed+1);
      // random velocities in [-rand_vel, rand_vel]
      for (std::size_t i=0;i<S.N;++i){
        S.vx[i] = rng.uniform(-rand_vel, rand_vel);
        S.vy[i] = rng.uniform(-rand_vel, rand_vel);
        S.vz[i] = rng.uniform(-rand_vel, rand_vel);
      }
      // remove net linear momentum to avoid COM drift
      double Px=0, Py=0, Pz=0, Mtot=0;
      for (std::size_t i=0;i<S.N;++i){
        Px += S.m[i]*S.vx[i];
        Py += S.m[i]*S.vy[i];
        Pz += S.m[i]*S.vz[i];
        Mtot += S.m[i];
      }
      double vxcm = Px / std::max(1e-300, Mtot);
      double vycm = Py / std::max(1e-300, Mtot);
      double vzcm = Pz / std::max(1e-300, Mtot);
      for (std::size_t i=0;i<S.N;++i){
        S.vx[i] -= vxcm;
        S.vy[i] -= vycm;
        S.vz[i] -= vzcm;
      }
    }
    double eps2 = eps*eps;
    const double pairs_per_step = double(S.N) * double(S.N);
    Timer tim; tim.tic();
    // Baseline invariants for relative error
    double H0 = kinetic_energy(S.vx,S.vy,S.vz,S.m,S.N) + potential_energy(S.x,S.y,S.z,S.m,S.N,eps2);
    double L0[3]; angular_momentum(S.x,S.y,S.z,S.vx,S.vy,S.vz,S.m,S.N,L0);
    const double L0n = std::sqrt(L0[0]*L0[0]+L0[1]*L0[1]+L0[2]*L0[2]);
    FILE* fp = nullptr;
    if (!out_csv.empty()){
      fp = std::fopen(out_csv.c_str(), "w");
      if (fp) std::fprintf(fp, "step,t,H,Lnorm,rel_dE,rel_dL\n");
    }
    for (int s=0;s<steps;++s){
      if (plt.enabled && (s%plt.every==0)) plt.frame(S, s, s*dt);
      if (method=="yoshida4") step_yoshida4(dt, S, eps2, Bj);
      else if (method=="rk4") step_rk4(dt, S, eps2, Bj);
      else if (method=="verlet") step_verlet(dt, S, eps2, Bj);
      else { std::fprintf(stderr,"Unknown method\n"); return 3; }
      if (check && ((s%check_every)==0 || s==steps-1)){
        double H = kinetic_energy(S.vx,S.vy,S.vz,S.m,S.N) + potential_energy(S.x,S.y,S.z,S.m,S.N,eps2);
        double L[3]; angular_momentum(S.x,S.y,S.z,S.vx,S.vy,S.vz,S.m,S.N,L);
        double Lnorm = std::sqrt(L[0]*L[0]+L[1]*L[1]+L[2]*L[2]);
        double dH = std::abs(H-H0)/std::max(1.0,std::abs(H0));
        double dL = std::sqrt((L[0]-L0[0])*(L[0]-L0[0])+(L[1]-L0[1])*(L[1]-L0[1])+(L[2]-L0[2])*(L[2]-L0[2]));
        double dLrel = dL/std::max(1.0,L0n);
        double t = (s+1)*dt;
        std::printf("step %7d t=%.6g  H=%.12e  |L|=%.6e  rel|ΔH|=%.3e  rel|ΔL|=%.3e\n",
                    s+1, t, H, Lnorm, dH, dLrel);
        if (fp){ std::fprintf(fp, "%d,%.15g,%.15g,%.15g,%.15g,%.15g\n", s+1, t, H, Lnorm, dH, dLrel); }
      }
    }
    double sec = tim.toc();
    double pairs = pairs_per_step * steps;
    double pairs_per_sec = pairs / sec;
    // 粗い換算: 1 ペアあたり ~32 FLOPs と仮定して GFLOP/s を見積り
    double gflops = (pairs_per_sec * 32.0) / 1e9;
    std::printf("Done: N=%zu steps=%d dt=%.3g eps=%.3g method=%s\n", S.N, steps, dt, eps, method.c_str());
    std::printf("Time = %.3fs  Pairs/s = %.3e  ~GFLOP/s ≈ %.2f (rough)\n", sec, pairs_per_sec, gflops);
    if (fp) { std::fclose(fp); std::printf("[REPORT] CSV saved: %s\n", out_csv.c_str()); }
    plt.close();
  }

  aligned_free(S.x); aligned_free(S.y); aligned_free(S.z);
  aligned_free(S.vx); aligned_free(S.vy); aligned_free(S.vz);
  aligned_free(S.m); aligned_free(S.ax); aligned_free(S.ay); aligned_free(S.az);
  return 0;
}
