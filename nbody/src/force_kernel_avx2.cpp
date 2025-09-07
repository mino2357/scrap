#include "force_kernel_avx2.hpp"
#include <immintrin.h>
#include <algorithm>
#include <cmath>
#include <omp.h>

/*
  このファイルについて
  ----------------------------------------
  - 目的: N体相互作用で必要な加速度（および速度更新）を AVX2 を使って
    倍精度(double)で4要素同時に計算するカーネルを実装します。
  - OpenMP による並列 for と、AVX2/FMA によるベクトル化を併用します。

  SIMD ベクトル型のメモ
  ----------------------------------------
  - __m256d: 256bit 幅のベクトル = 倍精度(double) 4要素。
    例) _mm256_set1_pd(a) は a を4レーン全てに複製します。
  - __m256:  256bit 幅のベクトル = 単精度(float) 8要素。
  - __m128:  128bit 幅のベクトル = 単精度(float) 4要素。
  - レーン(要素)は [0..3] の4本で、ロード/ストアは配列の連続4要素に対応。

  主な命令のメモ（_mm256_* 系）
  ----------------------------------------
  - load/store:   _mm256_load_pd, _mm256_store_pd（32バイト境界に整列している前提）
                  整列が不明なら _mm256_loadu_pd, _mm256_storeu_pd を使用します。
  - set:          _mm256_set1_pd(a) は全レーンを a で埋める。
  - 算術:        _mm256_add_pd, _mm256_sub_pd, _mm256_mul_pd など。
  - FMA:          _mm256_fmadd_pd(a,b,c)   = a*b + c（丸め一回）
                  _mm256_fnmadd_pd(a,b,c)  = -(a*b) + c
  - ブロードキャスト: _mm256_broadcast_sd(ptr) = *ptr を4レーンに複製。
  - 変換:        _mm256_cvtpd_ps(double→float, 4要素) / _mm256_cvtps_pd(float→double, 4要素)
                  精度を落として単精度にした後、再び倍精度に戻す用途に使います。

  物理式のメモ
  ----------------------------------------
  - r^2 = |r_j - r_i|^2 + eps2（ソフトニング付き）
  - inv = 1/sqrt(r^2), inv3 = inv^3
  - a_i += m_j * inv3 * (r_j - r_i)
    自分自身(j==i)の項は rx=ry=rz=0 なので寄与ゼロです。

  rsqrt 近似のメモ
  ----------------------------------------
  - 倍精度の逆平方根 1/sqrt(x) は、単精度の rsqrt で初期値を作り、
    Newton-Raphson を2回かけて精度を上げています。
  - 更新式: y ← 0.5*y*(3 - x*y*y)。FMA で (3 - x*y*y) を評価します。
*/

// Mode flags
// g_vector_mode: 1=ベクトル化(AVX2)を使う, 0=スカラ版にフォールバック
// g_scalar_tiled: スカラ版のとき j 方向にタイル分割(Bj)するかの指定（1: する, 0: しない）
static int g_vector_mode = 1;
static int g_scalar_tiled = 0;
extern "C" int* __get_vector_mode(){ return &g_vector_mode; }
extern "C" int* __get_scalar_tiled(){ return &g_scalar_tiled; }

// 倍精度 4要素に対する 1/sqrt(x) の近似を返す。
// 手順:
//  1) x(double×4) を単精度(float×4)に縮退し、_mm_rsqrt_ps で粗い近似を得る。
//  2) 倍精度に戻し、Newton-Raphson を2回適用して精度を改善。
//  メモ: FMA の _fnmadd は -(a*b)+c を計算するため、(3 - x*y*y) 等を効率よく評価できます。
static inline __m256d rsqrt_dp_avx2(__m256d x) {
    __m128  xf   = _mm256_cvtpd_ps(x);
    __m128  y0f  = _mm_rsqrt_ps(xf);
    __m256d y    = _mm256_cvtps_pd(y0f);
    const __m256d half  = _mm256_set1_pd(0.5);
    const __m256d three = _mm256_set1_pd(3.0);
    __m256d yy = _mm256_mul_pd(y, y);
    __m256d t  = _mm256_fnmadd_pd(x, yy, three);
    y = _mm256_mul_pd(y, _mm256_mul_pd(half, t));
    yy = _mm256_mul_pd(y, y);
    t  = _mm256_fnmadd_pd(x, yy, three);
    y  = _mm256_mul_pd(y, _mm256_mul_pd(half, t));
    return y;
}

// 加速度 a(x) を計算して ax, ay, az に格納する。
// - 入力: 位置 (x,y,z), 質量 m, 粒子数 N, ソフトニング eps2, j方向タイル幅 Bj
// - ベクトル版: i を4つまとめて、j をスカラで走査（j の1要素を4レーンへブロードキャスト）
// - スカラ版: g_vector_mode==0 のときに実行。g_scalar_tiled==1 なら j をタイル分割。
void compute_accel_dp_avx2(const double* x, const double* y, const double* z, const double* m,
                           double* ax, double* ay, double* az,
                           std::size_t N, double eps2, std::size_t Bj)
{
    // Scalar fallback
    if (!g_vector_mode){
        if (g_scalar_tiled){
            // Scalar with j-tiling (Bj)
            #pragma omp parallel for schedule(static)
            for (ptrdiff_t i=0;i<(ptrdiff_t)N;++i){
                double xi=x[i], yi=y[i], zi=z[i];
                double axi=0,ayi=0,azi=0;
                for (std::size_t jb=0; jb<N; jb+=Bj){
                    std::size_t jn = std::min(jb+Bj, N);
                    for (std::size_t j=jb; j<jn; ++j){
                        double rx=x[j]-xi, ry=y[j]-yi, rz=z[j]-zi;
                        double r2=rx*rx+ry*ry+rz*rz+eps2;
                        double inv=1.0/std::sqrt(r2);
                        double inv3=inv*inv*inv;
                        double a=m[j]*inv3;
                        axi += a*rx; ayi += a*ry; azi += a*rz;
                    }
                }
                ax[i]=axi; ay[i]=ayi; az[i]=azi;
            }
        } else {
            // Scalar naive (no tiling)
            #pragma omp parallel for schedule(static)
            for (ptrdiff_t i=0;i<(ptrdiff_t)N;++i){
                double xi=x[i], yi=y[i], zi=z[i];
                double axi=0,ayi=0,azi=0;
                for (std::size_t j=0; j<N; ++j){
                    double rx=x[j]-xi, ry=y[j]-yi, rz=z[j]-zi;
                    double r2=rx*rx+ry*ry+rz*rz+eps2;
                    double inv=1.0/std::sqrt(r2);
                    double inv3=inv*inv*inv;
                    double a=m[j]*inv3;
                    axi += a*rx; ayi += a*ry; azi += a*rz;
                }
                ax[i]=axi; ay[i]=ayi; az[i]=azi;
            }
        }
        return;
    }

    // ベクトル版を使う場合は、まず出力ベクトルをゼロクリア
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t i=0;i<(ptrdiff_t)N;++i){ ax[i]=ay[i]=az[i]=0.0; }

    const __m256d veps2 = _mm256_set1_pd(eps2);
    ptrdiff_t Ivec = (ptrdiff_t)(N & ~std::size_t(3));

    // i を4要素ずつ処理（__m256d は double×4）
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t ib=0; ib<Ivec; ib+=4){
        // 入力は std::vector 由来の一時領域もあり得るため、ロードは非整列版を使用
        __m256d xi = _mm256_loadu_pd(x+ib);
        __m256d yi = _mm256_loadu_pd(y+ib);
        __m256d zi = _mm256_loadu_pd(z+ib);
        __m256d axi = _mm256_setzero_pd();
        __m256d ayi = _mm256_setzero_pd();
        __m256d azi = _mm256_setzero_pd();

        for (std::size_t jb=0; jb<N; jb+=Bj){
            std::size_t jn = std::min(jb+Bj, N);
            std::size_t j  = jb;
            for (; j+3<jn; j+=4){
                // j をスカラで進めつつ、各 j のスカラ値を4レーンにブロードキャスト。
                // r2 = rx*rx + ry*ry + rz*rz + eps2 を FMA で評価し、
                // inv = rsqrt(r2) → inv3 = inv^3 → a = m[j]*inv3 として加算。
                // j+0
                {
                    __m256d xj = _mm256_broadcast_sd(x + j + 0);
                    __m256d yj = _mm256_broadcast_sd(y + j + 0);
                    __m256d zj = _mm256_broadcast_sd(z + j + 0);
                    __m256d mj = _mm256_broadcast_sd(m + j + 0);
                    __m256d rx = _mm256_sub_pd(xj, xi);
                    __m256d ry = _mm256_sub_pd(yj, yi);
                    __m256d rz = _mm256_sub_pd(zj, zi);
                    __m256d r2 = _mm256_fmadd_pd(rx, rx, veps2);
                    r2         = _mm256_fmadd_pd(ry, ry, r2);
                    r2         = _mm256_fmadd_pd(rz, rz, r2);
                    __m256d inv = rsqrt_dp_avx2(r2);
                    __m256d inv3= _mm256_mul_pd(_mm256_mul_pd(inv,inv), inv);
                    __m256d a   = _mm256_mul_pd(mj, inv3);
                    axi = _mm256_fmadd_pd(a, rx, axi);
                    ayi = _mm256_fmadd_pd(a, ry, ayi);
                    azi = _mm256_fmadd_pd(a, rz, azi);
                }
                // j+1
                {
                    __m256d xj = _mm256_broadcast_sd(x + j + 1);
                    __m256d yj = _mm256_broadcast_sd(y + j + 1);
                    __m256d zj = _mm256_broadcast_sd(z + j + 1);
                    __m256d mj = _mm256_broadcast_sd(m + j + 1);
                    __m256d rx = _mm256_sub_pd(xj, xi);
                    __m256d ry = _mm256_sub_pd(yj, yi);
                    __m256d rz = _mm256_sub_pd(zj, zi);
                    __m256d r2 = _mm256_fmadd_pd(rx, rx, veps2);
                    r2         = _mm256_fmadd_pd(ry, ry, r2);
                    r2         = _mm256_fmadd_pd(rz, rz, r2);
                    __m256d inv = rsqrt_dp_avx2(r2);
                    __m256d inv3= _mm256_mul_pd(_mm256_mul_pd(inv,inv), inv);
                    __m256d a   = _mm256_mul_pd(mj, inv3);
                    axi = _mm256_fmadd_pd(a, rx, axi);
                    ayi = _mm256_fmadd_pd(a, ry, ayi);
                    azi = _mm256_fmadd_pd(a, rz, azi);
                }
                // j+2
                {
                    __m256d xj = _mm256_broadcast_sd(x + j + 2);
                    __m256d yj = _mm256_broadcast_sd(y + j + 2);
                    __m256d zj = _mm256_broadcast_sd(z + j + 2);
                    __m256d mj = _mm256_broadcast_sd(m + j + 2);
                    __m256d rx = _mm256_sub_pd(xj, xi);
                    __m256d ry = _mm256_sub_pd(yj, yi);
                    __m256d rz = _mm256_sub_pd(zj, zi);
                    __m256d r2 = _mm256_fmadd_pd(rx, rx, veps2);
                    r2         = _mm256_fmadd_pd(ry, ry, r2);
                    r2         = _mm256_fmadd_pd(rz, rz, r2);
                    __m256d inv = rsqrt_dp_avx2(r2);
                    __m256d inv3= _mm256_mul_pd(_mm256_mul_pd(inv,inv), inv);
                    __m256d a   = _mm256_mul_pd(mj, inv3);
                    axi = _mm256_fmadd_pd(a, rx, axi);
                    ayi = _mm256_fmadd_pd(a, ry, ayi);
                    azi = _mm256_fmadd_pd(a, rz, azi);
                }
                // j+3
                {
                    __m256d xj = _mm256_broadcast_sd(x + j + 3);
                    __m256d yj = _mm256_broadcast_sd(y + j + 3);
                    __m256d zj = _mm256_broadcast_sd(z + j + 3);
                    __m256d mj = _mm256_broadcast_sd(m + j + 3);
                    __m256d rx = _mm256_sub_pd(xj, xi);
                    __m256d ry = _mm256_sub_pd(yj, yi);
                    __m256d rz = _mm256_sub_pd(zj, zi);
                    __m256d r2 = _mm256_fmadd_pd(rx, rx, veps2);
                    r2         = _mm256_fmadd_pd(ry, ry, r2);
                    r2         = _mm256_fmadd_pd(rz, rz, r2);
                    __m256d inv = rsqrt_dp_avx2(r2);
                    __m256d inv3= _mm256_mul_pd(_mm256_mul_pd(inv,inv), inv);
                    __m256d a   = _mm256_mul_pd(mj, inv3);
                    axi = _mm256_fmadd_pd(a, rx, axi);
                    ayi = _mm256_fmadd_pd(a, ry, ayi);
                    azi = _mm256_fmadd_pd(a, rz, azi);
                }
            }
            for (; j<jn; ++j){
                // 端数の j を1つずつ処理
                __m256d xj = _mm256_broadcast_sd(x + j);
                __m256d yj = _mm256_broadcast_sd(y + j);
                __m256d zj = _mm256_broadcast_sd(z + j);
                __m256d mj = _mm256_broadcast_sd(m + j);
                __m256d rx = _mm256_sub_pd(xj, xi);
                __m256d ry = _mm256_sub_pd(yj, yi);
                __m256d rz = _mm256_sub_pd(zj, zi);
                __m256d r2 = _mm256_fmadd_pd(rx, rx, veps2);
                r2         = _mm256_fmadd_pd(ry, ry, r2);
                r2         = _mm256_fmadd_pd(rz, rz, r2);
                __m256d inv = rsqrt_dp_avx2(r2);
                __m256d inv3= _mm256_mul_pd(_mm256_mul_pd(inv,inv), inv);
                __m256d a   = _mm256_mul_pd(mj, inv3);
                axi = _mm256_fmadd_pd(a, rx, axi);
                ayi = _mm256_fmadd_pd(a, ry, ayi);
                azi = _mm256_fmadd_pd(a, rz, azi);
            }
        }
        // 出力は整列していることが多いが、安全のため非整列ストア
        _mm256_storeu_pd(ax+ib, axi);
        _mm256_storeu_pd(ay+ib, ayi);
        _mm256_storeu_pd(az+ib, azi);
    }

    for (std::size_t i=Ivec; i<N; ++i){
        // 端数の i（ベクトル幅に満たない分）はスカラで処理
        double xi=x[i], yi=y[i], zi=z[i];
        double axi=0,ayi=0,azi=0;
        for (std::size_t jb=0; jb<N; jb+=Bj){
            std::size_t jn = std::min(jb+Bj, N);
            for (std::size_t j=jb; j<jn; ++j){
                double rx=x[j]-xi, ry=y[j]-yi, rz=z[j]-zi;
                double r2=rx*rx+ry*ry+rz*rz+eps2;
                double inv=1.0/std::sqrt(r2);
                double inv3=inv*inv*inv;
                double a=m[j]*inv3;
                axi += a*rx; ayi += a*ry; azi += a*rz;
            }
        }
        ax[i]=axi; ay[i]=ayi; az[i]=azi;
    }
}

// Fused kick（速度の更新）: v += tau * a(x)
// ----------------------------------------
// 役割
//  - 時間刻み tau の分だけ、位置 x によって決まる加速度 a(x) を速度 v に加算します。
//  - つまり a(x) を明示配列に書かず、その場で計算→tau を掛け→v に加算する「融合」カーネル。
//
// 引数
//  - x,y,z: 粒子の位置（長さ N の double 配列）
//  - m    : 粒子の質量（長さ N）
//  - vx,vy,vz: 速度ベクトル（更新先）
//  - N    : 粒子数
//  - eps2 : ソフトニング項 ε^2（近接特異性の緩和）
//  - tau  : 時間刻み Δt
//  - Bj   : j 方向タイル幅（キャッシュ効率向上のため、j を Bj ごとのブロックで処理）
//
// 計算内容（物理）
//  - r_ij = r_j - r_i,  r2 = |r_ij|^2 + eps2
//  - inv = 1/sqrt(r2)（rsqrt の近似 + Newton で高速）
//  - inv3 = inv^3,  a_contrib = m_j * inv3 * r_ij
//  - v_i += tau * a_contrib を全 j に対して加算
//  - 自分自身(j==i)の項は r_ij=0 なので寄与は 0 となり、特別扱いは不要
//
// 実装（AVX2/FMA ベクトル化）
//  - i を 4 粒子ずつ（__m256d = double×4）処理。
//  - j はスカラで進め、各 j のスカラ値を _mm256_broadcast_sd で 4 レーンに複製。
//  - r2 は FMA（_mm256_fmadd_pd）で rx*rx + ry*ry + rz*rz + eps2 を 3 命令で評価。
//  - inv = rsqrt_dp_avx2(r2) は 1/sqrt(r2) を倍精度で近似（Newton 2 回）。
//  - a = m_j * inv^3、tau_a = tau * a として、v += tau_a * r を FMA で加算。
//
// 性能・数値の注意
//  - Bj による j タイリングで、j ブロックの x,y,z,m を L1/L2 に載せやすくします。
//  - _mm256_load_pd は 32B 整列を前提とします（本コードは整列ロード）。
//    整列が保証できない環境では loadu/storeu への置換が必要です。
//  - FMA は丸めが一回で済むため、速度だけでなく数値誤差の抑制にも寄与します。
//  - eps2 により近接相互作用の発散を抑制。tau の大きさは積分器側で管理します。
void kick_accumulate_dp_avx2(const double* x, const double* y, const double* z, const double* m,
                             double* vx, double* vy, double* vz,
                             std::size_t N, double eps2, double tau, std::size_t Bj)
{
    /*
      実装詳細（レジスタ配置・依存鎖・スケジューリングの観点）
      ------------------------------------------------------------------
      レジスタ/値の役割（ベクトル化パスの i ブロック内で不変）
        - xi, yi, zi:   i 側位置（__m256d, double×4）
        - vxi, vyi, vzi: i 側速度の蓄積（__m256d, double×4）
        - veps2, vtau:  スカラーを4レーンへ複製した定数（__m256d）

      j ループ内（各 j または j+k）での一時レジスタ
        - xj, yj, zj, mj:   j 側スカラー値をブロードキャスト（__m256d）
        - rx, ry, rz:       相対位置 r_ij（__m256d）
        - r2:               |r|^2 + eps2（FMA 3発で評価）（__m256d）
        - inv:              1/sqrt(r2)（rsqrt_dp_avx2）（__m256d）
        - inv3:             inv*inv*inv（__m256d）
        - a:                m_j * inv3（__m256d）
        - tau_a:            tau * a（__m256d）

      依存鎖（概略）
        rx/ry/rz → r2(FMA×3) → inv(rsqrt+Newton) → inv3 → a → tau_a →
        v更新(FMA×3)
      ここがクリティカルパス。rsqrt+Newton が比較的長いため、ある程度の ILP を
      露出する工夫が有効です（本実装では j を 4 回ループ内でまとめて処理）。

      ILP とアンロール
        - 現状: 外側で j を 4 つずつ (j, j+1, j+2, j+3) 回すことで、
          それぞれの r2/inv/inv3 の計算が（vxi/vyi/vzi の蓄積を除けば）
          相互に独立になり、FMA/乗算/ブロードキャストをポート分散しやすい。
        - vxi/vyi/vzi の連続 FMA は依存鎖を作るため、必要に応じて
          4 回分を一時変数へ部分和として溜めてから合算する「木型リダクション」も選択肢。
          （本コードは可読性のため逐次蓄積）

      マイクロアーキテクチャの一般的な目安（参考）
        - FMA のスループットは多くのマイクロアーキテクチャで 2/サイクル（port0/1）だが、
          レイテンシは 4～5 サイクル程度。3 発連続の r2 評価は FMA の依存鎖になる。
        - _mm256_broadcast_sd は L1 ヒット時 1/cycle 程度で供給可能（実装依存）。
        - rsqrt_dp_avx2 内の Newton 2 回は mul/fnmadd が主で、
          1 呼び出しあたり ~6 FLOPs + 変換を含む依存鎖が加わる。
        - これらを隠すには、（本実装のように）j を小さくアンロールして独立作業を並べるのが有効。

      メモリ局所性と Bj
        - j ブロックの作業集合は 4 配列 × 8B × Bj = 32B×Bj。
          L1(32KB) を強く意識するなら Bj ≲ 1000 が目安（他のデータやインフライト命令もある）。
          実際には L2/L3 とのトレードオフで 256～2048 程度をスイープして最適化するのが良い。
        - i 側はベクトルレジスタに常駐（xi/yi/zi, vxi/vyi/vzi）し、j 側はブロードキャストで 8B 読み。
          帯域要求は比較的軽く、計算密度が高いカーネル。

      整列ロード/ストアと移植性
        - _mm256_load_pd/_mm256_store_pd は 32B 整列が前提。配列は 32B 境界に配置すること。
          不明な場合や外部入力では *_loadu_pd に切り替える。
    */
    // Scalar fallback
    if (!g_vector_mode){
        if (g_scalar_tiled){
            // Scalar with j-tiling (Bj)
            #pragma omp parallel for schedule(static)
            for (ptrdiff_t i=0;i<(ptrdiff_t)N;++i){
                double xi_=x[i], yi_=y[i], zi_=z[i];
                double vxi=vx[i], vyi=vy[i], vzi=vz[i];
                for (std::size_t jb=0; jb<N; jb+=Bj){
                    std::size_t jn = std::min(jb+Bj, N);
                    for (std::size_t j=jb; j<jn; ++j){
                        double rx=x[j]-xi_, ry=y[j]-yi_, rz=z[j]-zi_;
                        double r2=rx*rx+ry*ry+rz*rz+eps2;
                        double inv=1.0/std::sqrt(r2);
                        double inv3=inv*inv*inv;
                        double a=m[j]*inv3, tau_a=tau*a;
                        vxi += tau_a*rx; vyi += tau_a*ry; vzi += tau_a*rz;
                    }
                }
                vx[i]=vxi; vy[i]=vyi; vz[i]=vzi;
            }
        } else {
            // Scalar naive (no tiling)
            #pragma omp parallel for schedule(static)
            for (ptrdiff_t i=0;i<(ptrdiff_t)N;++i){
                double xi_=x[i], yi_=y[i], zi_=z[i];
                double vxi=vx[i], vyi=vy[i], vzi=vz[i];
                for (std::size_t j=0; j<N; ++j){
                    double rx=x[j]-xi_, ry=y[j]-yi_, rz=z[j]-zi_;
                    double r2=rx*rx+ry*ry+rz*rz+eps2;
                    double inv=1.0/std::sqrt(r2);
                    double inv3=inv*inv*inv;
                    double a=m[j]*inv3, tau_a=tau*a;
                    vxi += tau_a*rx; vyi += tau_a*ry; vzi += tau_a*rz;
                }
                vx[i]=vxi; vy[i]=vyi; vz[i]=vzi;
            }
        }
        return;
    }

    // ベクトル化パス
    const __m256d veps2 = _mm256_set1_pd(eps2); // 全レーンに eps2
    const __m256d vtau  = _mm256_set1_pd(tau);  // 全レーンに tau
    ptrdiff_t Ivec = (ptrdiff_t)(N & ~std::size_t(3));

    // i を 4 要素ずつ処理。各 i ブロックに対し、Bj 幅の j ブロックを順に当てる。
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t ib=0; ib<Ivec; ib+=4){
        // i 側の位置と速度をロード（double×4）。ロードは整列前提。
        __m256d xi = _mm256_load_pd(x+ib);
        __m256d yi = _mm256_load_pd(y+ib);
        __m256d zi = _mm256_load_pd(z+ib);
        __m256d vxi = _mm256_load_pd(vx+ib);
        __m256d vyi = _mm256_load_pd(vy+ib);
        __m256d vzi = _mm256_load_pd(vz+ib);

        for (std::size_t jb=0; jb<N; jb+=Bj){
            std::size_t jn = std::min(jb+Bj, N);
            std::size_t j  = jb;
            for (; j+3<jn; j+=4){
                // j+0..3 をひとまとめに処理（同じパターンを4回繰り返す）
                for (int k=0;k<4;++k){
                    // j 側のスカラ値を4レーンにブロードキャスト
                    __m256d xj = _mm256_broadcast_sd(x + j + k);
                    __m256d yj = _mm256_broadcast_sd(y + j + k);
                    __m256d zj = _mm256_broadcast_sd(z + j + k);
                    __m256d mj = _mm256_broadcast_sd(m + j + k);
                    // 相対位置 r_ij
                    __m256d rx = _mm256_sub_pd(xj, xi);
                    __m256d ry = _mm256_sub_pd(yj, yi);
                    __m256d rz = _mm256_sub_pd(zj, zi);
                    // r2 = rx^2 + ry^2 + rz^2 + eps2（FMAで依存連鎖を短縮）
                    __m256d r2 = _mm256_fmadd_pd(rx, rx, veps2);
                    r2         = _mm256_fmadd_pd(ry, ry, r2);
                    r2         = _mm256_fmadd_pd(rz, rz, r2);
                    // 1/sqrt(r2) を倍精度で近似（Newton 2回）
                    __m256d inv = rsqrt_dp_avx2(r2);
                    // inv^3, a = m_j * inv^3
                    __m256d inv3= _mm256_mul_pd(_mm256_mul_pd(inv,inv), inv);
                    __m256d a   = _mm256_mul_pd(mj, inv3);
                    // tau_a = tau * a
                    __m256d tau_a = _mm256_mul_pd(vtau, a);
                    // v += tau_a * r を FMA で加算
                    vxi = _mm256_fmadd_pd(tau_a, rx, vxi);
                    vyi = _mm256_fmadd_pd(tau_a, ry, vyi);
                    vzi = _mm256_fmadd_pd(tau_a, rz, vzi);
                }
            }
            for (; j<jn; ++j){
                // 端数の j を1つずつ処理
                __m256d xj = _mm256_broadcast_sd(x + j);
                __m256d yj = _mm256_broadcast_sd(y + j);
                __m256d zj = _mm256_broadcast_sd(z + j);
                __m256d mj = _mm256_broadcast_sd(m + j);
                __m256d rx = _mm256_sub_pd(xj, xi);
                __m256d ry = _mm256_sub_pd(yj, yi);
                __m256d rz = _mm256_sub_pd(zj, zi);
                __m256d r2 = _mm256_fmadd_pd(rx, rx, veps2);
                r2         = _mm256_fmadd_pd(ry, ry, r2);
                r2         = _mm256_fmadd_pd(rz, rz, r2);
                __m256d inv = rsqrt_dp_avx2(r2);
                __m256d inv3= _mm256_mul_pd(_mm256_mul_pd(inv,inv), inv);
                __m256d a   = _mm256_mul_pd(mj, inv3);
                __m256d tau_a = _mm256_mul_pd(vtau, a);
                vxi = _mm256_fmadd_pd(tau_a, rx, vxi);
                vyi = _mm256_fmadd_pd(tau_a, ry, vyi);
                vzi = _mm256_fmadd_pd(tau_a, rz, vzi);
            }
        }
        // 4 粒子分の速度をストア（整列ストア）
        _mm256_store_pd(vx+ib, vxi);
        _mm256_store_pd(vy+ib, vyi);
        _mm256_store_pd(vz+ib, vzi);
    }

    for (std::size_t i=Ivec; i<N; ++i){
        // 端数の i はスカラで処理（式はベクトル版と同じ）
        double xi_=x[i], yi_=y[i], zi_=z[i];
        double vxi=vx[i], vyi=vy[i], vzi=vz[i];
        for (std::size_t jb=0; jb<N; jb+=Bj){
            std::size_t jn = std::min(jb+Bj, N);
            for (std::size_t j=jb; j<jn; ++j){
                double rx=x[j]-xi_, ry=y[j]-yi_, rz=z[j]-zi_;
                double r2=rx*rx+ry*ry+rz*rz+eps2;
                double inv=1.0/std::sqrt(r2);
                double inv3=inv*inv*inv;
                double a=m[j]*inv3, tau_a=tau*a;
                vxi += tau_a*rx; vyi += tau_a*ry; vzi += tau_a*rz;
            }
        }
        vx[i]=vxi; vy[i]=vyi; vz[i]=vzi;
    }
}
