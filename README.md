# Collatz Tail One-Click Repro

Pure Python

This wraps the supplied supplementary bundle into a single-command entry point.

## Run

```bash
python run_all.py
```

The script will:

1. verify the core input files exist
2. regenerate the `R4` certificate table from `dangerous_words.csv`
3. print the worked example for `node_id = 1`
4. run the `Nmax = 10^12` violation check
5. write regenerated outputs under `outputs/`

## Main data + scripts

Inside `Collatz_Tail_Supplementary_Bundle_v2/`:

- `data/dangerous_words.csv`
- `data/seeds_b24_depth3.csv`
- `data/results_r4.csv`
- `data/violations_N1e12.csv`
- `scripts/reproduce_cited_files.py`
- `scripts/holographic_checkmate_65325_v4.py`

## Expected headline check

The one-click run should confirm:

- rows = 201
- odd count = 201
- min bitlen(R4) = 328
- count(R4 <= 10^12) = 0

and print the worked example for `node_id = 1`, including the 99-digit `R4`.

## Zenodo

doi:10.5281/zenodo.18125856
