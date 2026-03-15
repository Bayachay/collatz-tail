#!/usr/bin/env python3
"""
reproduce_cited_files.py

Reproduces the key CSV artifacts referenced in the paper from the two primary inputs:
  - dangerous_words.csv         (201 pump words with A, r, cL)
  - seeds_b24_depth3.csv        (325 depth-3 seeds with nL0 and t)

Outputs:
  - results_r4_repro.csv        (recomputed R4 and basic flags)
  - violations_N1e12_repro.csv  (empty if no hits)

Also prints the "certificate summary" reviewers can quote:
  rows, odd/even counts, min digits, min bitlen, and count(R4<=1e12).
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Tuple

import pandas as pd


def compute_R4(A: int, r: int, c: int, L: int = 64) -> int:
    a = pow(3, L)
    M = 1 << A
    inv_a = pow(a, -1, M)
    R = int(r)
    c = int(c)
    for k in range(1, 4):
        q = (a * R + c) // M
        denom = pow(M, k - 1)
        num = R - q
        if num % denom != 0:
            raise ValueError(f"Lift divisibility failed at k={k} (A={A})")
        rhs = (num // denom) % M
        s = (rhs * inv_a) % M
        R = R + pow(M, k) * s
    return R


def certificate_summary(r4_series, bound: int = 10**12) -> str:
    r4 = r4_series.astype(object)
    odd = int((r4 % 2 == 1).sum())
    even = int((r4 % 2 == 0).sum())
    rmin = int(r4.min())
    return "\n".join([
        f"rows: {len(r4)}",
        f"odd count: {odd}",
        f"even count: {even}",
        f"min R4: {rmin}",
        f"min digits: {len(str(rmin))}",
        f"min bitlen: {rmin.bit_length()}",
        f"any R4 <= 1e12: {bool((r4 <= bound).any())}",
        f"count R4 <= 1e12: {int((r4 <= bound).sum())}",
    ])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dangerous", default="dangerous_words.csv")
    ap.add_argument("--seeds", default="seeds_b24_depth3.csv")
    ap.add_argument("--L", type=int, default=64)
    ap.add_argument("--B", type=int, default=24)
    ap.add_argument("--Nmax", type=int, default=10**12)
    ap.add_argument("--out_r4", default="results_r4_repro.csv")
    ap.add_argument("--out_viol", default="violations_N1e12_repro.csv")
    args = ap.parse_args()

    dw = pd.read_csv(args.dangerous)
    if not {"A","r","cL"}.issubset(set(dw.columns)):
        raise ValueError("dangerous_words.csv must contain columns: A, r, cL")

    # Compute R4
    R4_list = []
    mod_list = []
    for _, row in dw.iterrows():
        A = int(row["A"])
        R4 = compute_R4(A, int(row["r"]), int(row["cL"]), L=args.L)
        R4_list.append(R4)
        mod_list.append(1 << (4*A))

    out = dw.copy()
    out["R4_decimal"] = [str(x) for x in R4_list]
    out["R4"] = out["R4_decimal"]
    out["mod_4A"] = [str(x) for x in mod_list]
    out["in_range"] = False  # recomputed below
    out["impossible"] = True # preserved as a legacy field; not used for conclusions

    # Range flag at N=2e7 for the original "dataset checkmate" statement
    N0 = 2 * 10**7
    out["in_range"] = [x <= N0 for x in R4_list]

    out_path = Path(args.out_r4)
    out.to_csv(out_path, index=False)

    print("Wrote:", out_path)
    print("\nresults_r4 certificate summary:")
    print(certificate_summary(pd.Series(R4_list)))

    # Reproduce 65,325 lift congruence test in n-space
    seeds = pd.read_csv(args.seeds)
    if not {"nL0","t"}.issubset(set(seeds.columns)) and not {"nL0","t0"}.issubset(set(seeds.columns)):
        raise ValueError("seeds_b24_depth3.csv must contain columns: nL0 and t (or t0)")
    tcol = "t" if "t" in seeds.columns else "t0"

    K = 2 * pow(3, args.L)
    B = args.B
    Nmax = int(args.Nmax)

    viol_rows = []
    for si, srow in seeds.iterrows():
        nL0 = int(srow["nL0"])
        t0 = int(srow[tcol])
        base = nL0 + K * t0
        for wi, (A, R4, mod) in enumerate(zip(dw["A"].astype(int), R4_list, mod_list)):
            g = 1 << (B + 1)        # gcd(step,mod)
            rhs = (R4 - base) % mod
            if rhs % g != 0:
                continue
            mod2 = mod // g
            inv_3L = pow(pow(3, args.L), -1, mod2)
            s0 = ((rhs // g) * inv_3L) % mod2
            t_hit = t0 + (1 << B) * s0
            n_hit = nL0 + K * t_hit
            if n_hit <= Nmax:
                viol_rows.append((si, wi, nL0, t0, int(A), int(R4), int(s0), int(t_hit), int(n_hit), Nmax))

    outv = Path(args.out_viol)
    with outv.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["seed_idx","word_idx","nL0","t0","A","R4","s0","t_hit","n_hit","Nmax"])
        w.writerows(viol_rows)

    print("Wrote:", outv, f"(hits={len(viol_rows)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
