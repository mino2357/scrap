#!/usr/bin/env bash
# Unified runner for nbody: build, autosweep, sweep, heatmap, selftest
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
cd "$here"

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
export OMP_PLACES=${OMP_PLACES:-cores}
export OMP_PROC_BIND=${OMP_PROC_BIND:-close}

BIN=./build/nbody

usage() {
  cat <<USAGE
Usage: ./run.sh <command> [options]

Commands
  build [--clean] [--cc <cxx>] [--no-lld]
  run [nbody options...] [--bytes-per-pair B] [--fpp F] [--peak-gflops G] [--mem-gbs M]
  selftest [--steps N] [--dt DT] [--eps EPS] [--method M] [--threads T]
  autosweep [--N N] [--dt DT] [--eps EPS] [--threads T] [--list CSV]
  sweep [--N N] [--steps S] [--dt DT] [--eps EPS] [--threads T] [--kernel K]
        [--coarse-min A] [--coarse-max B] [--coarse-step S] 
        [--fine-radius R] [--fine-step S] [--out FILE]
  heatmap [--Ns CSV] [--bj-min A] [--bj-max B] [--bj-step S]
          [--steps S] [--dt DT] [--eps EPS] [--threads T] [--kernel K]
          [--out FILE]
  live [--N N] [--steps S] [--dt DT] [--eps EPS] [--method M]
       [--Bj B] [--threads T] [--kernel K]
       [--plot-every K] [--plot-limit M] [--plot-axes xy|xz|yz]

Examples
  ./run.sh build --clean
  ./run.sh run --N 4096 --steps 200 --dt 1e-3 --eps 1e-3 --method yoshida4 --Bj 512 --threads 4 --kernel fused --check 1
  ./run.sh autosweep --N 4096 --list 256,320,384,448,512,576,640,704,768,896
  ./run.sh sweep --N 4096 --steps 200 --kernel fused --out sweep_bj_200.csv
  ./run.sh heatmap --Ns 1024,2048,4096,8192 --out bjN_heatmap.csv
  ./run.sh live --N 2048 --steps 2000 --dt 1e-3 --eps 1e-3 --Bj 512 --plot-every 20 --plot-axes xy
USAGE
}

build() {
  local clean=0 cxx=${CXX:-clang++} use_lld=1
  while [[ $# -gt 0 ]]; do
    case $1 in
      --clean) clean=1; shift ;;
      --cc) cxx=$2; shift 2 ;;
      --no-lld) use_lld=0; shift ;;
      *) echo "[build] unknown option: $1"; usage; exit 2 ;;
    esac
  done
  [[ $clean -eq 1 ]] && rm -rf build
  local ldflags=""; [[ $use_lld -eq 1 ]] && ldflags="-DCMAKE_EXE_LINKER_FLAGS=-fuse-ld=lld"
  cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER="$cxx" $ldflags
  cmake --build build -j
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
    pass_args+=(--N 4096 --steps 200 --dt 1e-3 --eps 1e-3 --method yoshida4 --Bj 512 --threads 4 --kernel fused --check 1)
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
  "$BIN" --N "$N" --steps 1 --dt "$dt" --eps "$eps" \
          --method yoshida4 --Bj 512 --threads "$threads" \
          --autosweep "$list" --kernel fused
}

sweep() {
  ensure_build
  local N=${N:-4096} steps=${STEPS:-200} dt=${DT:-1e-3} eps=${EPS:-1e-3} threads=${THREADS:-$OMP_NUM_THREADS}
  local kernel=${KERNEL:-fused}
  local cmin=${COARSE_MIN:-128} cmax=${COARSE_MAX:-2048} cstep=${COARSE_STEP:-128}
  local fr=${FINE_RADIUS:-256} fstep=${FINE_STEP:-16}
  local out=${OUT:-sweep_bj_200.csv}

  export OMP_NUM_THREADS="$threads"
  echo "# N=$N steps=$steps dt=$dt eps=$eps threads=$threads kernel=$kernel" | tee "$out"
  echo "#Bj  PairsPerSec" | tee -a "$out"

  run_one() {
    local BJ="$1" out pps
    out=$("$BIN" --N "$N" --steps "$steps" --dt "$dt" --eps "$eps" \
                  --method yoshida4 --Bj "$BJ" --threads "$threads" \
                  --kernel "$kernel" --check 0 2>&1)
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

heatmap() {
  ensure_build
  local Ns_csv=${NS:-${NSs:-"1024,2048,4096,8192"}}
  local bj_min=${BJ_MIN:-128} bj_max=${BJ_MAX:-3072} bj_step=${BJ_STEP:-128}
  local steps=${STEPS:-200} dt=${DT:-1e-3} eps=${EPS:-1e-3}
  local threads=${THREADS:-1} kernel=${KERNEL:-fused}
  local out=${OUT:-bjN_heatmap.csv}

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
  local N=${N:-2048} steps=${STEPS:-2000} dt=${DT:-1e-3} eps=${EPS:-1e-3}
  local method=${METHOD:-yoshida4} Bj=${BJ:-512} threads=${THREADS:-$OMP_NUM_THREADS}
  local kernel=${KERNEL:-fused}
  local plot_every=${PLOT_EVERY:-20} plot_limit=${PLOT_LIMIT:-4096} plot_axes=${PLOT_AXES:-xy}

  export OMP_NUM_THREADS="$threads"
  echo "[live] N=$N steps=$steps dt=$dt eps=$eps Bj=$Bj threads=$threads kernel=$kernel plot=$plot_axes/$plot_every limit=$plot_limit"
  "$BIN" --N "$N" --steps "$steps" --dt "$dt" --eps "$eps" \
         --method "$method" --Bj "$Bj" --threads "$threads" --kernel "$kernel" \
         --plot live --plot_every "$plot_every" --plot_limit "$plot_limit" --plot_axes "$plot_axes"
}

cmd=${1:-}
case "$cmd" in
  build) shift; build "$@" ;;
  run) shift; run_cmd "$@" ;;
  selftest) shift; selftest "$@" ;;
  autosweep) shift; autosweep "$@" ;;
  sweep) shift; sweep "$@" ;;
  heatmap) shift; heatmap "$@" ;;
  live) shift; live "$@" ;;
  ""|-h|--help|help) usage ;;
  *) echo "Unknown command: $cmd"; usage; exit 2 ;;
esac
