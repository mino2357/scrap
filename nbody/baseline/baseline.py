#!/usr/bin/env python3
import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path


def run(cmd, cwd=None):
    proc = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")
    return proc.stdout


def default_cases():
    # Deterministic, CPU-agnostic: scalar mode, single thread
    common = {
        "N": 100,
        "steps": 4000,
        "dt": 5e-4,
        "eps": 1e-2,
        "Bj": 512,
        "threads": 1,
        "seed": 1337,
        "rand_vel": 0.05,
        "mode": "scalar",
        "check_every": 200,
    }
    methods = ["yoshida4", "verlet", "rk4"]
    kernels = ["two_pass", "fused"]
    cases = []
    for ker in kernels:
        for meth in methods:
            cid = f"{ker}_{meth}"
            c = {"id": cid, "kernel": ker, "method": meth}
            c.update(common)
            cases.append(c)
    return cases


def run_case(bin_path: Path, outdir: Path, case: dict) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / f"{case['id']}.csv"
    args = [
        str(bin_path),
        "--N", str(case["N"]),
        "--steps", str(case["steps"]),
        "--dt", f"{case['dt']}",
        "--eps", f"{case['eps']}",
        "--Bj", str(case["Bj"]),
        "--threads", str(case["threads"]),
        "--check", "1",
        "--check_every", str(case["check_every"]),
        "--seed", str(case["seed"]),
        "--rand_vel", f"{case['rand_vel']}",
        "--kernel", case["kernel"],
        "--mode", case["mode"],
        "--method", case["method"],
        "--out_csv", str(csv_path),
    ]
    run(args)

    # Parse CSV: step,t,H,Lnorm,rel_dE,rel_dL
    H0 = None
    L0 = None
    max_rel_dE = 0.0
    max_rel_dL = 0.0
    Hend = None
    Lend = None
    with open(csv_path, newline="") as f:
        it = csv.reader(f)
        header = next(it)
        for row in it:
            step = int(row[0])
            H = float(row[2])
            L = float(row[3])
            rel_dE = float(row[4])
            rel_dL = float(row[5])
            if H0 is None:
                H0 = H
                L0 = L
            if rel_dE > max_rel_dE:
                max_rel_dE = rel_dE
            if rel_dL > max_rel_dL:
                max_rel_dL = rel_dL
            Hend = H
            Lend = L

    return {
        "id": case["id"],
        "kernel": case["kernel"],
        "method": case["method"],
        "params": {
            k: case[k] for k in ["N", "steps", "dt", "eps", "Bj", "threads", "seed", "rand_vel", "mode", "check_every"]
        },
        "metrics": {
            "H0": H0,
            "L0": L0,
            "H_end": Hend,
            "L_end": Lend,
            "max_rel_dE": max_rel_dE,
            "max_rel_dL": max_rel_dL,
        },
    }


def approx_equal(a, b, rtol=1e-10, atol=5e-13):
    if a is None or b is None:
        return False
    return abs(a - b) <= (atol + rtol * max(abs(a), abs(b)))


def compare_results(current: dict, baseline: dict, rtolE=1e-10, atolE=5e-13, rtolL=1e-10, atolL=5e-13) -> list:
    diffs = []
    for cid, cur in current.items():
        if cid not in baseline:
            diffs.append((cid, "missing_baseline", None, None))
            continue
        base = baseline[cid]
        # Compare key metrics
        for k, rtol, atol in [
            ("max_rel_dE", rtolE, atolE),
            ("max_rel_dL", rtolL, atolL),
            ("H0", rtolE, atolE),
            ("L0", rtolL, atolL),
            ("H_end", rtolE, atolE),
            ("L_end", rtolL, atolL),
        ]:
            a = cur["metrics"][k]
            b = base["metrics"][k]
            if not approx_equal(a, b, rtol=rtol, atol=atol):
                diffs.append((cid, k, a, b))
    return diffs


def main():
    parser = argparse.ArgumentParser(description="Run and check nbody baseline")
    parser.add_argument("action", choices=["run", "update", "check"], help="run: produce current.json, update: write baseline.json, check: compare")
    parser.add_argument("--bin", dest="bin_path", default=None, help="Path to nbody binary (default: nbody/build/nbody)")
    parser.add_argument("--out", dest="out_dir", default="nbody/out/baseline", help="Output directory for CSV and JSON")
    parser.add_argument("--cases", dest="cases_file", default=None, help="Optional cases.json file")
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[2]
    bin_path = Path(args.bin_path) if args.bin_path else (repo / "nbody" / "build" / "nbody")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load cases
    if args.cases_file:
        with open(args.cases_file) as f:
            cases = json.load(f)
    else:
        cases = default_cases()

    # Run all cases
    current = {}
    for c in cases:
        res = run_case(bin_path, out_dir, c)
        current[res["id"]] = res

    current_path = out_dir / "current.json"
    with open(current_path, "w") as f:
        json.dump(current, f, indent=2, sort_keys=True)

    if args.action == "run":
        print(f"Wrote {current_path}")
        return 0

    baseline_path = repo / "nbody" / "baseline" / "baseline.json"
    if args.action == "update":
        with open(baseline_path, "w") as f:
            json.dump(current, f, indent=2, sort_keys=True)
        print(f"Updated baseline at {baseline_path}")
        return 0

    # check
    with open(baseline_path) as f:
        baseline = json.load(f)
    diffs = compare_results(current, baseline)
    if not diffs:
        print("Baseline check passed: all metrics within tolerance")
        return 0
    print("Baseline check FAILED:")
    for cid, key, a, b in diffs:
        print(f"  case={cid} metric={key} current={a} baseline={b}")
    return 2


if __name__ == "__main__":
    sys.exit(main())

