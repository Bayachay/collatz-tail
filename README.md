#!/usr/bin/env python3
"""mc_pump_survival_split.py

Monte Carlo (multilevel splitting) estimate of pump-block survival probabilities.

Model:
- Odd-only Collatz map: n -> (3n+1) / 2^v2(3n+1)
- A 'block' = L odd-steps (default L=64). Block valuation A = sum v2(3n_i+1) over the L odd steps.
- A block is 'pump' if A <= Amax (default 101).

Goal:
Estimate S(k) = P(first k blocks are pump) under uniform sampling of odd residues mod 2^B.

Notes:
- Multilevel splitting is used because S(k) gets tiny fast.
- This is an estimator, not a proof. For proof-style results, replace sampling with exact dyadic counting on cylinders.

Outputs:
- CSV with mean and 5–95% quantiles across independent runs.
"""

from __future__ import annotations
import argparse, random, math
import numpy as np
import pandas as pd

def v2(x:int)->int:
    return (x & -x).bit_length()-1

def odd_step(n:int):
    m = 3*n + 1
    a = v2(m)
    return m >> a, a

def block_A(n:int, L:int):
    A=0
    for _ in range(L):
        n,a = odd_step(n)
        A += a
    return n, A

def splitting_once(N:int, B:int, L:int, Amax:int, K:int, seed:int):
    rng = random.Random(seed)
    pop = [ (rng.getrandbits(B) | 1) or 1 for _ in range(N) ]
    logS = 0.0
    S = [1.0]
    for level in range(1, K+1):
        survivors=[]
        for n in pop:
            n2,A = block_A(n, L)
            if A <= Amax:
                survivors.append(n2)
        m=len(survivors)
        if m == 0:
            S.extend([0.0]*(K-level+1))
            break
        logS += math.log(m/N)
        S.append(math.exp(logS))
        pop = [ survivors[rng.randrange(m)] for _ in range(N) ]
    return np.array(S)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--B', type=int, default=24)
    ap.add_argument('--L', type=int, default=64)
    ap.add_argument('--Amax', type=int, default=101)
    ap.add_argument('--K', type=int, default=10)
    ap.add_argument('--N', type=int, default=8000)
    ap.add_argument('--reps', type=int, default=80)
    ap.add_argument('--seed', type=int, default=1000)
    ap.add_argument('--out', default='mc_pump_survival_estimates_B24_L64_Amax101.csv')
    args = ap.parse_args()

    Ss=[]
    for r in range(args.reps):
        Ss.append(splitting_once(args.N, args.B, args.L, args.Amax, args.K, args.seed+r))
    Ss=np.stack(Ss)
    mean=Ss.mean(axis=0)
    q05=np.quantile(Ss,0.05,axis=0)
    q95=np.quantile(Ss,0.95,axis=0)
    df=pd.DataFrame({
        'k_blocks': np.arange(Ss.shape[1]),
        'S_mean': mean,
        'S_q05': q05,
        'S_q95': q95
    })
    df.to_csv(args.out, index=False)
    print('wrote', args.out)
    print(df)

if __name__=='__main__':
    main()
