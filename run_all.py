from pathlib import Path
import subprocess
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parent
BUNDLE = ROOT / "Collatz_Tail_Supplementary_Bundle_v2"
DATA = BUNDLE / "data"
SCRIPTS = BUNDLE / "scripts"
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)

def info(msg):
    print(f"[collatz-tail] {msg}")

def run_py(script, args):
    cmd = [sys.executable, str(script)] + args
    info("running: " + " ".join(cmd))
    return subprocess.run(cmd, cwd=BUNDLE, check=False)

def print_worked_example(results_path):
    df = pd.read_csv(results_path, dtype={"R4_decimal": str})
    row = None
    for key in ("node_id", "word_id"):
        if key in df.columns:
            mask = df[key].astype(str) == "1"
            if mask.any():
                row = df[mask].iloc[0]
                break
    if row is None:
        row = df.iloc[0]
    r4 = int(row["R4_decimal"])
    print()
    print("worked example")
    if "node_id" in row.index:
        print("node_id =", row["node_id"])
    if "A" in row.index:
        print("A =", row["A"])
    for maybe in ("r", "r_decimal"):
        if maybe in row.index:
            print("r =", row[maybe])
            break
    print("R4 =", r4)
    print("digits(R4) =", len(str(r4)))
    print("bitlen(R4) =", r4.bit_length())
    print()

def print_certificate_summary(results_path):
    df = pd.read_csv(results_path, dtype={"R4_decimal": str})
    r4 = df["R4_decimal"].map(int)
    print()
    print("certificate summary")
    print("rows =", len(df))
    print("odd count =", int((r4 % 2 == 1).sum()))
    print("even count =", int((r4 % 2 == 0).sum()))
    print("min digits(R4) =", len(str(r4.min())))
    print("min bitlen(R4) =", r4.min().bit_length())
    print("count(R4 <= 10^12) =", int((r4 <= 10**12).sum()))
    print()

def main():
    need = [
        DATA / "dangerous_words.csv",
        DATA / "seeds_b24_depth3.csv",
        SCRIPTS / "reproduce_cited_files.py",
        SCRIPTS / "holographic_checkmate_65325_v4.py",
    ]
    missing = [p.name for p in need if not p.exists()]
    if missing:
        info("missing required files: " + ", ".join(missing))
        sys.exit(1)

    repro_results = OUT / "results_r4_repro.csv"
    repro_viol = OUT / "violations_N1e12_repro.csv"

    run_py(
        SCRIPTS / "reproduce_cited_files.py",
        [
            "--dangerous", str(DATA / "dangerous_words.csv"),
            "--seeds", str(DATA / "seeds_b24_depth3.csv"),
            "--out_r4", str(repro_results),
            "--out_viol", str(repro_viol),
        ],
    )

    results_for_summary = repro_results if repro_results.exists() else (DATA / "results_r4.csv")
    print_certificate_summary(results_for_summary)
    print_worked_example(results_for_summary)

    run_py(
        SCRIPTS / "holographic_checkmate_65325_v4.py",
        [
            "--seeds", str(DATA / "seeds_b24_depth3.csv"),
            "--words", str(DATA / "dangerous_words.csv"),
            "--out", str(OUT / "violations_N1e12_from_holo.csv"),
            "--B", "24",
            "--L", "64",
            "--Nmax", "1000000000000",
        ],
    )

    info("done, outputs are in ./outputs")

if __name__ == "__main__":
    main()
