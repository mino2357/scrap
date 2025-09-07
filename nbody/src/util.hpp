#pragma once
#include <immintrin.h>
#include <cstddef>
#include <cstdlib>
#include <cstdint>
#include <chrono>
#include <random>
#include <cstring>

/*
  実行時ユーティリティ
  ----------------------------------------
  - aligned_malloc/aligned_free: AVX2 向けに 32B/64B 整列メモリを確保/解放。
  - enable_FTZ_DAZ: x86 の FTZ/DAZ を有効化（subnormal/denormal の性能ペナルティを回避）。
  - Timer: 単純な経過時間測定（steady_clock）。
  - RNG: 簡易な一様乱数 [a,b] 生成（mt19937_64）。
*/

inline void* aligned_malloc(std::size_t bytes, std::size_t align=64) {
#if defined(_MSC_VER)
  return _aligned_malloc(bytes, align);
#else
  void* p = nullptr;
  if (posix_memalign(&p, align, bytes) != 0) return nullptr;
  return p;
#endif
}
inline void aligned_free(void* p){
#if defined(_MSC_VER)
  _aligned_free(p);
#else
  free(p);
#endif
}

inline void enable_FTZ_DAZ(){
  unsigned mxcsr = _mm_getcsr();
  mxcsr |= 0x8040; // DAZ | FTZ （Denormals-Are-Zero, Flush-To-Zero）
  _mm_setcsr(mxcsr);
}

struct Timer {
  using clock = std::chrono::steady_clock;
  clock::time_point t0;
  void tic(){ t0 = clock::now(); }
  double toc() const {
    auto t1 = clock::now();
    return std::chrono::duration<double>(t1 - t0).count();
  }
};

struct RNG {
  std::mt19937_64 eng;
  RNG(uint64_t seed):eng(seed){}
  double uniform(double a, double b){
    std::uniform_real_distribution<double> d(a,b);
    return d(eng);
  }
};
