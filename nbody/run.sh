#!/usr/bin/env bash
# Unified runner for nbody: build, autosweep, sweep, heatmap, selftest
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
cd "$here"

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
export OMP_PLACES=${OMP_PLACES:-cores}
export OMP_PROC_BIND=${OMP_PROC_BIND:-close}

BUILD_DIR=${BUILD_DIR:-build}
BIN=./${BUILD_DIR}/nbody

usage() {
  cat <<USAGE
Usage: ./run.sh <command> [options]

Commands
  build [--clean] [--cc <cxx>] [--no-lld] [--opt default|strong|weak] [--dir DIR]
  run [nbody options...] [--bytes-per-pair B] [--fpp F] [--peak-gflops G] [--mem-gbs M]
  compare [common nbody options...]   # runs vector vs scalar and summarizes speedup
  compare_scalar [common nbody options...]   # runs scalar vs scalar_tiled
  compare_opt [common nbody options...]   # builds weak/strong in separate dirs and compares
  selftest [--steps N] [--dt DT] [--eps EPS] [--method M] [--threads T]
  autosweep [--N N] [--dt DT] [--eps EPS] [--threads T] [--list CSV]
  sweep [--N N] [--steps S] [--dt DT] [--eps EPS] [--threads T] [--kernel K]
        [--coarse-min A] [--coarse-max B] [--coarse-step S] 
        [--fine-radius R] [--fine-step S] [--out FILE]
  heatmap [--Ns CSV] [--bj-min A] [--bj-max B] [--bj-step S]
          [--steps S] [--dt DT] [--eps EPS] [--threads T] [--kernel K]
          [--out FILE]
  live [nbody options...]   # forwards args; adds --plot live if absent

Examples
  ./run.sh build --clean
  ./run.sh run --N 4096 --steps 200 --dt 1e-3 --eps 1e-3 --method yoshida4 --Bj 512 --threads 4 --kernel fused --check 1
  ./run.sh compare --N 4096 --steps 200 --dt 1e-3 --eps 1e-3 --method yoshida4 --Bj 512 --threads 4 --kernel fused
  ./run.sh autosweep --N 4096 --list 256,320,384,448,512,576,640,704,768,896
  ./run.sh sweep --N 4096 --steps 200 --kernel fused --out sweep_bj_200.csv
  ./run.sh heatmap --Ns 1024,2048,4096,8192 --out bjN_heatmap.csv
  ./run.sh live --N 2048 --steps 2000 --dt 1e-3 --eps 1e-3 --Bj 512 --plot-every 20 --plot-axes xy
USAGE
}

build() {
  local clean=0 cxx=${CXX:-clang++} use_lld=1 opt=${OPT:-} dir=${DIR:-$BUILD_DIR}
  while [[ $# -gt 0 ]]; do
    case $1 in
      --clean) clean=1; shift ;;
      --cc) cxx=$2; shift 2 ;;
      --no-lld) use_lld=0; shift ;;
      --opt) opt=$2; shift 2 ;;
      --dir) dir=$2; shift 2 ;;
      *) echo "[build] unknown option: $1"; usage; exit 2 ;;
    esac
  done
  [[ $clean -eq 1 ]] && rm -rf "$dir"
  local ldflags=""; [[ $use_lld -eq 1 ]] && ldflags="-DCMAKE_EXE_LINKER_FLAGS=-fuse-ld=lld"
  local optflag=""; [[ -n "$opt" ]] && optflag="-DNBODY_OPT_PRESET=$opt"
  cmake -S . -B "$dir" -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER="$cxx" $ldflags $optflag
  cmake --build "$dir" -j
  BUILD_DIR="$dir"; BIN="./$dir/nbody"
}

ensure_build() {
  if [[ ! -x "$BIN" ]]; then
    echo "[info] building nbody..."
    build "$@"
  fi
}

selftest() {
  ensure_build
  local steps=${STEPS:-20000} dt=${DT:-1e-4} eps=${EPS:-1e-6} method=${METHOD:-yoshida4} threads=${THREADS:-$OMP_NUM_THREADS}
  "$BIN" --selftest --method "$method" --steps "$steps" --dt "$dt" --eps "$eps" --threads "$threads"
}

run_cmd() {
  ensure_build
  # Defaults match the typical example if no args are provided
  local fpp=32 bpp=8 peak_gflops="" mem_gbs=""
  local pass_args=()
  if [[ $# -eq 0 ]]; then
    pass_args+=(--N 4096 --steps 200 --dt 1e-3 --eps 1e-3 --method yoshida4 --Bj 512 --threads 4 --kernel fused --mode vector --check 1)
  else
    while [[ $# -gt 0 ]]; do
      case $1 in
        --fpp) fpp=$2; shift 2 ;;
        --bytes-per-pair|--bpp) bpp=$2; shift 2 ;;
        --peak-gflops) peak_gflops=$2; shift 2 ;;
        --mem-gbs) mem_gbs=$2; shift 2 ;;
        *) pass_args+=("$1"); shift ;;
      esac
    done
  fi

  # Normalize common hyphenated/alias flags
  for i in "${!pass_args[@]}"; do
    case "${pass_args[$i]}" in
      --plot-every) pass_args[$i]=--plot_every ;;
      --plot-limit) pass_args[$i]=--plot_limit ;;
      --plot-axes)  pass_args[$i]=--plot_axes  ;;
      --bj)         pass_args[$i]=--Bj         ;;
    esac
  done

  # Run and capture output
  local out
  out=$("$BIN" "${pass_args[@]}" 2>&1 | tee /dev/fd/3 3>&1)

  # Parse metrics
  local sec pairs_per_sec N steps method eps Bj threads kernel
  sec=$(printf "%s\n" "$out" | awk -F 'Time = ' '/Time = /{split($2,a,"s"); t=a[1]} END{if(t!="") print t}')
  pairs_per_sec=$(printf "%s\n" "$out" | awk -F 'Pairs/s = ' '/Pairs\/s/{split($2,a," "); p=a[1]} END{if(p!="") print p}')
  N=$(printf "%s\n" "$out" | awk -F 'Done: N=' '/Done: N=/{sub(/ steps.*/,"",$2); print $2}' | tail -n1)
  steps=$(printf "%s\n" "$out" | awk -F 'steps=' '/Done: N=/{sub(/ dt.*/,"",$2); print $2}' | tail -n1)
  method=$(printf "%s\n" "$out" | awk -F 'method=' '/Done: N=/{print $2}' | tail -n1)
  eps=$(printf "%s\n" "$out" | awk -F ' eps=' '/Done: N=/{sub(/ method.*/,"",$2); print $2}' | tail -n1)
  Bj=$(printf "%s\n" "$out" | awk -F ' Bj=' '/\[INFO\] kernel=/{sub(/  threads.*/,"",$2); print $2}' | tail -n1)
  threads=$(printf "%s\n" "$out" | awk -F ' threads=' '/\[INFO\] kernel=/{print $2}' | tail -n1)
  kernel=$(printf "%s\n" "$out" | awk -F 'kernel=' '/\[INFO\] kernel=/{sub(/  Bj.*/,"",$2); print $2}' | tail -n1)

  # Compute extended metrics
  if [[ -n "$pairs_per_sec" && -n "$sec" ]]; then
    awk -v pps="$pairs_per_sec" -v sec="$sec" -v fpp="$fpp" -v bpp="$bpp" -v peak="$peak_gflops" -v bw="$mem_gbs" -v N="$N" -v steps="$steps" -v Bj="$Bj" -v thr="$threads" -v ker="$kernel" -v method="$method" -v eps="$eps" '
      BEGIN {
        pairs_total = "";
        if (N != "" && steps != "") {
          # Beware of overflow; print as scientific
          nt = N+0; st = steps+0; pairs_total = sprintf("%.6e", nt*nt*st);
        }
        gflops = pps * fpp / 1e9;
        gbs = pps * bpp / 1e9;
        ai = fpp / bpp;
        printf("\n[SUMMARY]\n");
        if (method != "") printf("- config: method=%s kernel=%s Bj=%s threads=%s eps=%s\n", method, ker, Bj, thr, eps);
        if (N != "" && steps != "") printf("- size: N=%s steps=%s\n", N, steps);
        printf("- timing: elapsed=%.3fs  pairs/s=%.6e\n", sec, pps);
        if (pairs_total != "") printf("- pairs: total=%s\n", pairs_total);
        printf("- compute: fpp=%d  GFLOP/s=%.3f\n", fpp, gflops);
        printf("- memory: bytes/pair=%d  est GB/s=%.3f  AI=%.2f flop/byte\n", bpp, gbs, ai);
        if (peak != "" && bw != "") {
          roof = bw * ai; eff_peak = (roof < peak ? roof : peak);
          comp_util = (peak>0 ? 100.0*gflops/peak : 0);
          mem_util = (bw>0 ? 100.0*gbs/bw : 0);
          printf("- roofline: peak=%.1f GFLOP/s  mem=%.1f GB/s  limit=%.1f GFLOP/s\n", peak, bw, eff_peak);
          printf("- utilization: compute=%.1f%%  memory=%.1f%%\n", comp_util, mem_util);
        }
      }
    '
  fi
}

selftest() {
  ensure_build
  local steps=${STEPS:-20000} dt=${DT:-1e-4} eps=${EPS:-1e-6} method=${METHOD:-yoshida4} threads=${THREADS:-$OMP_NUM_THREADS}
  "$BIN" --selftest --method "$method" --steps "$steps" --dt "$dt" --eps "$eps" --threads "$threads"
}

autosweep() {
  ensure_build
  local N=${N:-4096} dt=${DT:-1e-3} eps=${EPS:-1e-3} threads=${THREADS:-$OMP_NUM_THREADS} list=${LIST:-256,320,384,448,512,576,640,704,768,896}
  while [[ $# -gt 0 ]]; do
    case $1 in
      --N) N=$2; shift 2 ;;
      --dt) dt=$2; shift 2 ;;
      --eps) eps=$2; shift 2 ;;
      --threads) threads=$2; shift 2 ;;
      --list) list=$2; shift 2 ;;
      *) shift ;;
    esac
  done
  "$BIN" --N "$N" --steps 1 --dt "$dt" --eps "$eps" \
          --method yoshida4 --Bj 512 --threads "$threads" \
          --autosweep "$list" --kernel fused
}

sweep() {
  ensure_build
  local N=${N:-4096} steps=${STEPS:-200} dt=${DT:-1e-3} eps=${EPS:-1e-3} threads=${THREADS:-$OMP_NUM_THREADS}
  local kernel=${KERNEL:-fused}
  local mode=${MODE:-vector}
  local cmin=${COARSE_MIN:-128} cmax=${COARSE_MAX:-2048} cstep=${COARSE_STEP:-128}
  local fr=${FINE_RADIUS:-256} fstep=${FINE_STEP:-16}
  local out=${OUT:-sweep_bj_200.csv}

  # Parse CLI overrides
  while [[ $# -gt 0 ]]; do
    case $1 in
      --N) N=$2; shift 2 ;;
      --steps) steps=$2; shift 2 ;;
      --dt) dt=$2; shift 2 ;;
      --eps) eps=$2; shift 2 ;;
      --threads) threads=$2; shift 2 ;;
      --kernel) kernel=$2; shift 2 ;;
      --mode) mode=$2; shift 2 ;;
      --coarse-min) cmin=$2; shift 2 ;;
      --coarse-max) cmax=$2; shift 2 ;;
      --coarse-step) cstep=$2; shift 2 ;;
      --fine-radius) fr=$2; shift 2 ;;
      --fine-step) fstep=$2; shift 2 ;;
      --out) out=$2; shift 2 ;;
      *) shift ;;
    esac
  done

  export OMP_NUM_THREADS="$threads"
  echo "# N=$N steps=$steps dt=$dt eps=$eps threads=$threads kernel=$kernel" | tee "$out"
  echo "#Bj  PairsPerSec" | tee -a "$out"

  run_one() {
    local BJ="$1" out pps
    out=$("$BIN" --N "$N" --steps "$steps" --dt "$dt" --eps "$eps" \
                  --method yoshida4 --Bj "$BJ" --threads "$threads" \
                  --kernel "$kernel" --mode "$mode" --check 0 2>&1)
    pps=$(printf "%s\n" "$out" | awk -F 'Pairs/s = ' '/Pairs\/s/ {split($2,a," "); print a[1]}' | tail -n1)
    [[ -z "$pps" ]] && { echo "ERROR: failed at Bj=$BJ"; exit 1; }
    echo "$BJ $pps"
  }

  local best_bj=0 best_pps=0
  for BJ in $(seq "$cmin" "$cstep" "$cmax"); do
    read bj pps < <(run_one "$BJ")
    printf "%d %.6e\n" "$bj" "$pps" | tee -a "$out"
    if awk -v p="$pps" -v q="$best_pps" 'BEGIN{exit !(p>q)}'; then
      best_pps="$pps"; best_bj="$bj"
    fi
  done
  echo "# coarse_best Bj=$best_bj Pairs/s=$best_pps" | tee -a "$out"

  local fmin=$(( best_bj - fr )); (( fmin < cmin )) && fmin="$cmin"
  local fmax=$(( best_bj + fr )); (( fmax > cmax )) && fmax="$cmax"
  for BJ in $(seq "$fmin" "$fstep" "$fmax"); do
    read bj pps < <(run_one "$BJ")
    printf "%d %.6e\n" "$bj" "$pps" | tee -a "$out"
    if awk -v p="$pps" -v q="$best_pps" 'BEGIN{exit !(p>q)}'; then
      best_pps="$pps"; best_bj="$bj"
    fi
  done
  echo "# BEST Bj=$best_bj  Pairs/s=$best_pps" | tee -a "$out"
}

# Convenience wrapper to sweep Bj for scalar_tiled
sweep_scalar_tiled() {
  MODE=scalar_tiled sweep "$@"
}

heatmap() {
  ensure_build
  local Ns_csv=${NS:-${NSs:-"1024,2048,4096,8192"}}
  local bj_min=${BJ_MIN:-128} bj_max=${BJ_MAX:-3072} bj_step=${BJ_STEP:-128}
  local steps=${STEPS:-200} dt=${DT:-1e-3} eps=${EPS:-1e-3}
  local threads=${THREADS:-1} kernel=${KERNEL:-fused}
  local out=${OUT:-bjN_heatmap.csv}

  while [[ $# -gt 0 ]]; do
    case $1 in
      --Ns|--ns) Ns_csv=$2; shift 2 ;;
      --bj-min) bj_min=$2; shift 2 ;;
      --bj-max) bj_max=$2; shift 2 ;;
      --bj-step) bj_step=$2; shift 2 ;;
      --steps) steps=$2; shift 2 ;;
      --dt) dt=$2; shift 2 ;;
      --eps) eps=$2; shift 2 ;;
      --threads) threads=$2; shift 2 ;;
      --kernel) kernel=$2; shift 2 ;;
      --out) out=$2; shift 2 ;;
      *) shift ;;
    esac
  done

  IFS=',' read -r -a Ns <<< "$Ns_csv"
  {
    for N in "${Ns[@]}"; do
      for BJ in $(seq "$bj_min" "$bj_step" "$bj_max"); do
        out_line=$("$BIN" --N "$N" --steps "$steps" --dt "$dt" --eps "$eps" \
                           --method yoshida4 --Bj "$BJ" --threads "$threads" \
                           --kernel "$kernel" --check 0)
        pps=$(echo "$out_line" | awk -F 'Pairs/s = ' '/Pairs\/s/ {split($2,a," "); print a[1]}' | tail -n1)
        printf "%d,%d,%.6e\n" "$N" "$BJ" "$pps"
      done
    done
  } > "$out"
  echo "[heatmap] wrote $out"
}

live() {
  ensure_build
  if ! command -v gnuplot >/dev/null 2>&1; then
    echo "[warn] gnuplot not found. Install with: sudo apt install -y gnuplot" >&2
  fi
  local args=()
  if [[ $# -eq 0 ]]; then
    # sensible defaults for interactive live plotting
    args=(--N 2048 --steps 2000 --dt 1e-3 --eps 1e-3 \
          --method yoshida4 --Bj 512 --threads "${THREADS:-$OMP_NUM_THREADS}" --kernel fused \
          --plot live --plot_every 20 --plot_limit 4096 --plot_axes xy)
  else
    while [[ $# -gt 0 ]]; do args+=("$1"); shift; done
    # add --plot live if not specified
    local has_plot=0
    for a in "${args[@]}"; do [[ "$a" == "--plot" ]] && has_plot=1; done
    if [[ $has_plot -eq 0 ]]; then args+=(--plot live); fi
    # normalize common hyphenated aliases to program's underscore options
    for i in "${!args[@]}"; do
      case "${args[$i]}" in
        --plot-every) args[$i]=--plot_every ;;
        --plot-limit) args[$i]=--plot_limit ;;
        --plot-axes)  args[$i]=--plot_axes  ;;
      esac
    done
  fi
  echo -n "[live] launching: $BIN"; printf ' %q' "${args[@]}"; echo
  "$BIN" "${args[@]}"
}

compare() {
  ensure_build
  local args=("$@")
  # Normalize aliases
  for i in "${!args[@]}"; do
    case "${args[$i]}" in
      --plot-every) args[$i]=--plot_every ;;
      --plot-limit) args[$i]=--plot_limit ;;
      --plot-axes)  args[$i]=--plot_axes  ;;
      --bj)         args[$i]=--Bj         ;;
    esac
  done
  # Run vector
  local out_vec out_sca pps_vec pps_sca
  out_vec=$("$BIN" "${args[@]}" --mode vector 2>&1)
  pps_vec=$(echo "$out_vec" | awk -F 'Pairs/s = ' '/Pairs\/s/{split($2,a," "); print a[1]}' | tail -n1)
  # Run scalar
  out_sca=$("$BIN" "${args[@]}" --mode scalar 2>&1)
  pps_sca=$(echo "$out_sca" | awk -F 'Pairs/s = ' '/Pairs\/s/{split($2,a," "); print a[1]}' | tail -n1)
  echo "[compare] vector pairs/s: $pps_vec"
  echo "[compare] scalar pairs/s: $pps_sca"
  if [[ -n "$pps_vec" && -n "$pps_sca" ]]; then
    awk -v v="$pps_vec" -v s="$pps_sca" 'BEGIN{ if (s>0) printf("[compare] speedup vector/scalar = %.2fx\n", v/s); else print "[compare] speedup: N/A" }'
  fi
}

compare_opt() {
  # build and run weak vs strong in separate dirs
  local args=()
  while [[ $# -gt 0 ]]; do args+=("$1"); shift; done
  # weak
  cmake -S . -B build-weak -G Ninja -DCMAKE_BUILD_TYPE=Release -DNBODY_OPT_PRESET=weak >/dev/null
  cmake --build build-weak -j >/dev/null
  local out_w=$(./build-weak/nbody "${args[@]}" 2>&1)
  local pps_w=$(echo "$out_w" | awk -F 'Pairs/s = ' '/Pairs\/s/{split($2,a," "); print a[1]}' | tail -n1)
  # strong
  cmake -S . -B build-strong -G Ninja -DCMAKE_BUILD_TYPE=Release -DNBODY_OPT_PRESET=strong >/dev/null
  cmake --build build-strong -j >/dev/null
  local out_s=$(./build-strong/nbody "${args[@]}" 2>&1)
  local pps_s=$(echo "$out_s" | awk -F 'Pairs/s = ' '/Pairs\/s/{split($2,a," "); print a[1]}' | tail -n1)
  echo "[compare_opt] weak pairs/s:   $pps_w"
  echo "[compare_opt] strong pairs/s: $pps_s"
  if [[ -n "$pps_w" && -n "$pps_s" ]]; then
    awk -v s="$pps_s" -v w="$pps_w" 'BEGIN{ if (w>0) printf("[compare_opt] speedup strong/weak = %.2fx\n", s/w); else print "[compare_opt] speedup: N/A" }'
  fi
}

compare_scalar() {
  ensure_build
  local args=("$@")
  for i in "${!args[@]}"; do
    case "${args[$i]}" in
      --plot-every) args[$i]=--plot_every ;;
      --plot-limit) args[$i]=--plot_limit ;;
      --plot-axes)  args[$i]=--plot_axes  ;;
      --bj)         args[$i]=--Bj         ;;
    esac
  done
  local out_na out_tl pps_na pps_tl
  out_na=$("$BIN" "${args[@]}" --mode scalar 2>&1)
  pps_na=$(echo "$out_na" | awk -F 'Pairs/s = ' '/Pairs\/s/{split($2,a," "); print a[1]}' | tail -n1)
  out_tl=$("$BIN" "${args[@]}" --mode scalar_tiled 2>&1)
  pps_tl=$(echo "$out_tl" | awk -F 'Pairs/s = ' '/Pairs\/s/{split($2,a," "); print a[1]}' | tail -n1)
  echo "[compare_scalar] scalar pairs/s:       $pps_na"
  echo "[compare_scalar] scalar_tiled pairs/s: $pps_tl"
  if [[ -n "$pps_na" && -n "$pps_tl" ]]; then
    awk -v t="$pps_tl" -v n="$pps_na" 'BEGIN{ if (n>0) printf("[compare_scalar] speedup tiled/naive = %.2fx\n", t/n); else print "[compare_scalar] speedup: N/A" }'
  fi
}

cmd=${1:-}
case "$cmd" in
  build) shift; build "$@" ;;
  run) shift; run_cmd "$@" ;;
  selftest) shift; selftest "$@" ;;
  autosweep) shift; autosweep "$@" ;;
  sweep) shift; sweep "$@" ;;
  sweep_scalar_tiled) shift; sweep_scalar_tiled "$@" ;;
  heatmap) shift; heatmap "$@" ;;
  live) shift; live "$@" ;;
  compare) shift; compare "$@" ;;
  compare_scalar) shift; compare_scalar "$@" ;;
  compare_opt) shift; compare_opt "$@" ;;
  ""|-h|--help|help) usage ;;
  *) echo "Unknown command: $cmd"; usage; exit 2 ;;
esac
