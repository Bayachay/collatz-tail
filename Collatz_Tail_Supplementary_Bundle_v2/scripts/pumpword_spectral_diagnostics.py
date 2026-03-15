#!/usr/bin/env python3
"""
Pump-word spectral diagnostics (length-64 valuation words).

Inputs:
  - dangerous_words.csv with columns: node_id,count,A,r,cL,word
    where word is dot-separated valuations (64 integers).

Outputs:
  - pumpwords_spectral_signatures.csv
  - pumpwords_spectral_scatter.png

Metrics (computed on the parity-encoding b[i]=+1 if a[i] odd else -1):
  - max_abs_circular_autocorr: max_{k=1..63} | (1/64) sum_i b[i]*b[i+k] |
  - power_varratio: var(|FFT(b)|^2) / mean(|FFT(b)|^2)^2
  - flat_score: 1/(1 + power_varratio)

Also prints a short summary and the Rudin–Shapiro (length 64) reference values.
"""
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def parse_word(s):
    return [int(x) for x in str(s).strip().split('.') if x != ""]

def b_parity(a):
    return np.array([1 if (x % 2 == 1) else -1 for x in a], dtype=np.float64)

def maxabs_circ_autocorr(b):
    n=b.size
    F=np.fft.fft(b)
    ac=np.fft.ifft(F*np.conjugate(F)).real / n
    return float(np.max(np.abs(ac[1:])))

def power_varratio(b):
    F=np.fft.fft(b)
    power=(np.abs(F)**2)
    m=power.mean()
    v=power.var()
    return float(v/(m*m + 1e-30))

def flat_score_from_varratio(vr):
    return float(1.0/(1.0+vr))

def rudin_shapiro(n):
    N=1<<n
    a=np.empty(N, dtype=np.float64)
    for k in range(N):
        x=k
        prev=0
        cnt=0
        while x:
            bit=x & 1
            if bit==1 and prev==1:
                cnt += 1
            prev=bit
            x >>= 1
        a[k]=1.0 if (cnt % 2 == 0) else -1.0
    return a

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--words", default="dangerous_words.csv")
    ap.add_argument("--out_csv", default="pumpwords_spectral_signatures.csv")
    ap.add_argument("--out_plot", default="pumpwords_spectral_scatter.png")
    args=ap.parse_args()

    df=pd.read_csv(args.words)
    rows=[]
    for _,row in df.iterrows():
        a=parse_word(row["word"])
        if len(a)!=64:
            raise ValueError("Expected length 64 word, got %d" % len(a))
        b=b_parity(a)
        ma=maxabs_circ_autocorr(b)
        vr=power_varratio(b)
        fs=flat_score_from_varratio(vr)
        odd_frac=sum(x%2 for x in a)/64.0
        rows.append({
            "node_id": int(row["node_id"]),
            "A": int(row["A"]),
            "odd_frac": odd_frac,
            "max_abs_circular_autocorr": ma,
            "power_varratio": vr,
            "flat_score": fs,
        })

    out=pd.DataFrame(rows)
    out.to_csv(args.out_csv, index=False)

    rs=rudin_shapiro(6)
    rs_ma=maxabs_circ_autocorr(rs)
    rs_vr=power_varratio(rs)
    rs_fs=flat_score_from_varratio(rs_vr)

    # plot
    plt.figure(figsize=(6,4))
    plt.scatter(out["max_abs_circular_autocorr"], out["flat_score"], alpha=0.7, s=20)
    plt.scatter([rs_ma],[rs_fs], marker="x", s=80)
    plt.xlabel("max |circular autocorr| (parity encoding)")
    plt.ylabel("spectral flatness score")
    plt.title("Pump words: autocorrelation vs spectral flatness (length 64)")
    plt.tight_layout()
    plt.savefig(args.out_plot, dpi=200)
    plt.close()

    # print summary
    print("rows =", len(out))
    print("odd_frac: mean=%.3f min=%.3f max=%.3f" % (out["odd_frac"].mean(), out["odd_frac"].min(), out["odd_frac"].max()))
    print("max_abs_circular_autocorr: mean=%.3f min=%.3f max=%.3f" % (out["max_abs_circular_autocorr"].mean(), out["max_abs_circular_autocorr"].min(), out["max_abs_circular_autocorr"].max()))
    print("flat_score: mean=%.3f min=%.3f max=%.3f" % (out["flat_score"].mean(), out["flat_score"].min(), out["flat_score"].max()))
    print("Rudin–Shapiro length 64: max_abs_circular_autocorr=%.3f flat_score=%.3f" % (rs_ma, rs_fs))

if __name__=="__main__":
    main()
