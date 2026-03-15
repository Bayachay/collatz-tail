Collatz Tail Supplementary Bundle (v2)

This bundle contains the minimal tables and scripts needed to reproduce the paper's finite-depth certificates.

Contents

data/
  dangerous_words.csv
  seeds_b24_depth3.csv
  results_r4.csv
  violations_N1e12.csv

scripts/
  reproduce_cited_files.py
  holographic_checkmate_65325_v4.py

Optional diagnostics (not needed for any theorem)
  data/pumpwords_spectral_signatures.csv
  scripts/pumpword_spectral_diagnostics.py
  figures/fig_pumpwords_spectral_scatter.png

Quick start

  python scripts/reproduce_cited_files.py --dangerous data/dangerous_words.csv --seeds data/seeds_b24_depth3.csv --out_r4 data/results_r4_repro.csv --out_viol data/violations_N1e12_repro.csv

Expected headline check (from results_r4.csv)
  rows=201, odd=201, min_digits(R4)=99, min_bitlen(R4)=328, count(R4 <= 1e12)=0

  python scripts/holographic_checkmate_65325_v4.py --seeds data/seeds_b24_depth3.csv --words data/dangerous_words.csv --out data/violations_N1e12_repro.csv --B 24 --L 64 --Nmax 1000000000000

Integrity

See SHA256SUMS.txt for file checksums.
