#!/usr/bin/env python3
"""
holographic_checkmate_65325_v4.py

Correct "65,325 congruence" holographic test in n-space.

Given:
  - depth-3 seed parameters (nL0, t0) at precision B (default 24), with the node parameterization
        n = nL0 + K * t,   where K = 2 * 3^L  (default L=64)
  - pump words with (A, r, cL) so that the length-L block map is
        F_w(n) = (3^L n + cL) / 2^A
    and cylinder membership is n ≡ r (mod 2^A)

We compute the depth-4 repeat residue R4(w) modulo 2^(4A) via the lift in §2.1 of the paper.

For each (seed, word) pair we test whether there exists s >= 0 such that
    nL0 + K*(t0 + 2^B s) ≡ R4(w)   (mod 2^(4A))
and, if so, whether the minimal solution produces an integer n_hit <= Nmax.

Because K has v2(K)=1, this is a power-of-two linear congruence. Let:
    mod = 2^(4A)
    step = K * 2^B = (2*3^L)*2^B = 2^(B+1) * 3^L
    rhs  = (R4 - (nL0 + K*t0)) mod mod

Then solutions exist iff rhs ≡ 0 (mod 2^(B+1)).
When solvable:
    s0 ≡ (rhs / 2^(B+1)) * (3^L)^{-1}   (mod 2^(4A-B-1))

We take the minimal nonnegative s0 and form:
    t_hit = t0 + 2^B * s0
    n_hit = nL0 + K * t_hit

If n_hit <= Nmax, we record a violation.

Inputs:
  seeds CSV: columns nL0 and t (or t0)
  words CSV: columns A, r, cL (and optional node_id/word_id)

Output:
  violations CSV: header only if no hits.

Example:
  python holographic_checkmate_65325_v4.py --seeds seeds_b24_depth3.csv --words dangerous_words.csv --out violations_N1e12.csv --B 24 --L 64 --Nmax 1000000000000
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List


def parse_int(x: str) -> int:
    x = x.strip()
    return int(x, 0)


@dataclass(frozen=True)
class Seed:
    nL0: int
    t0: int


@dataclass(frozen=True)
class Word:
    A: int
    r: int
    c: int


def read_rows(path: Path) -> List[dict]:
    with path.open("r", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def load_seeds(path: Path) -> List[Seed]:
    rows = read_rows(path)
    if not rows:
        raise ValueError(f"Seeds CSV empty: {path}")
    if "nL0" not in rows[0]:
        raise ValueError("Seeds CSV missing column nL0")
    tcol = "t0" if "t0" in rows[0] else ("t" if "t" in rows[0] else None)
    if tcol is None:
        raise ValueError("Seeds CSV missing column t0 (or alias t)")
    out = []
    for r in rows:
        out.append(Seed(nL0=parse_int(r["nL0"]), t0=parse_int(r[tcol])))
    return out


def load_words(path: Path) -> List[Word]:
    rows = read_rows(path)
    if not rows:
        raise ValueError(f"Words CSV empty: {path}")
    for col in ("A", "r", "cL"):
        if col not in rows[0]:
            raise ValueError(f"Words CSV missing column {col}")
    out = []
    for r in rows:
        out.append(Word(A=int(r["A"]), r=parse_int(r["r"]), c=parse_int(r["cL"])))
    return out


def compute_R4(A: int, r: int, c: int, L: int = 64) -> int:
    """Compute the minimal representative R4 modulo 2^(4A) using the §2.1 lift."""
    a = pow(3, L)
    M = 1 << A
    inv_a = pow(a, -1, M)
    R = r
    for k in range(1, 4):
        q = (a * R + c) // M
        denom = pow(M, k - 1)
        num = R - q
        if num % denom != 0:
            raise ValueError(f"Lift divisibility failed at k={k} for A={A}")
        rhs = (num // denom) % M
        s = (rhs * inv_a) % M
        R = R + pow(M, k) * s
    return R


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", required=True)
    ap.add_argument("--words", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--B", type=int, default=24)
    ap.add_argument("--L", type=int, default=64)
    ap.add_argument("--Nmax", type=int, default=10**12)
    args = ap.parse_args()

    seeds = load_seeds(Path(args.seeds))
    words = load_words(Path(args.words))

    B = args.B
    L = args.L
    Nmax = args.Nmax
    K = 2 * pow(3, L)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Precompute per-word constants
    word_cache = []
    for w in words:
        mod = 1 << (4 * w.A)
        R4 = compute_R4(w.A, w.r, w.c, L=L)
        # gcd(step, mod) = 2^(B+1) since step = 2^(B+1)*3^L and 4A >= 328
        g = 1 << (B + 1)
        mod2 = mod // g
        inv_3L = pow(pow(3, L), -1, mod2)  # odd inverse modulo power of two
        word_cache.append((w.A, mod, R4, g, mod2, inv_3L))

    with out_path.open("w", newline="", encoding="utf-8") as f:
        wtr = csv.writer(f)
        wtr.writerow(["seed_idx", "word_idx", "nL0", "t0", "A", "R4", "s0", "t_hit", "n_hit", "Nmax"])

        hits = 0
        for si, s in enumerate(seeds):
            base = s.nL0 + K * s.t0
            for wi, (A, mod, R4, g, mod2, inv_3L) in enumerate(word_cache):
                rhs = (R4 - base) % mod
                if rhs % g != 0:
                    continue
                s0 = ((rhs // g) * inv_3L) % mod2
                t_hit = s.t0 + (1 << B) * s0
                n_hit = s.nL0 + K * t_hit
                if n_hit <= Nmax:
                    hits += 1
                    wtr.writerow([si, wi, s.nL0, s.t0, A, R4, s0, t_hit, n_hit, Nmax])

    print(f"seeds={len(seeds)} words={len(words)} pairs={len(seeds)*len(words)} hits={hits}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
