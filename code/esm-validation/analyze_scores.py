#!/usr/bin/env python3
"""
Statistical analysis of ESM2 pseudo-perplexity scores.
Wilcoxon rank-sum tests + Cohen's d, matching paper's statistical framework.
"""

import csv
import math
from pathlib import Path
from scipy import stats
import numpy as np

RESULTS_PATH = Path(__file__).resolve().parent / "results" / "esm2_pseudoppl_scores.csv"

# Load data
rows = []
with open(RESULTS_PATH) as f:
    for row in csv.DictReader(f):
        row["pseudo_perplexity"] = float(row["pseudo_perplexity"])
        row["mean_log_likelihood"] = float(row["mean_log_likelihood"])
        rows.append(row)

families = ["RRM", "SH3", "WW", "Kunitz", "zf-C2H2", "PDZ", "Pkinase", "Defensin_beta"]

print(f"{'Family':<16} {'SA gen ppl':>12} {'Stored ppl':>12} {'Wilcoxon p':>12} {'Cohen d':>10} {'Interpretation'}")
print("-" * 80)

for fam in families:
    sa = [r["pseudo_perplexity"] for r in rows if r["family"] == fam and r["source"] == "SA_generation"]
    st = [r["pseudo_perplexity"] for r in rows if r["family"] == fam and r["source"] == "Stored"]

    if not sa or not st:
        print(f"{fam:<16} MISSING DATA")
        continue

    sa_arr = np.array(sa)
    st_arr = np.array(st)

    # Wilcoxon rank-sum (two-sided)
    stat, p_val = stats.mannwhitneyu(sa_arr, st_arr, alternative="two-sided")

    # Cohen's d (SA - Stored), negative means SA is lower (better)
    pooled_std = math.sqrt(((len(sa_arr) - 1) * sa_arr.std(ddof=1)**2 +
                            (len(st_arr) - 1) * st_arr.std(ddof=1)**2) /
                           (len(sa_arr) + len(st_arr) - 2))
    d = (sa_arr.mean() - st_arr.mean()) / pooled_std if pooled_std > 0 else 0

    # Significance label
    if p_val < 0.001:
        sig = "***"
    elif p_val < 0.01:
        sig = "**"
    elif p_val < 0.05:
        sig = "*"
    else:
        sig = "n.s."

    # Interpretation
    if sig == "n.s.":
        interp = "Indistinguishable"
    elif d < 0:
        interp = f"SA BETTER ({sig})"
    else:
        interp = f"Stored better ({sig})"

    print(f"{fam:<16} {sa_arr.mean():8.2f}+/-{sa_arr.std():.2f} {st_arr.mean():8.2f}+/-{st_arr.std():.2f} {p_val:12.4e} {d:10.3f} {interp}")

# Also do mean log-likelihood
print(f"\n{'='*80}")
print(f"Mean log-likelihood (higher = more protein-like)")
print(f"{'='*80}")
print(f"{'Family':<16} {'SA gen mll':>12} {'Stored mll':>12} {'Wilcoxon p':>12} {'Cohen d':>10} {'Interpretation'}")
print("-" * 80)

for fam in families:
    sa = [r["mean_log_likelihood"] for r in rows if r["family"] == fam and r["source"] == "SA_generation"]
    st = [r["mean_log_likelihood"] for r in rows if r["family"] == fam and r["source"] == "Stored"]

    if not sa or not st:
        continue

    sa_arr = np.array(sa)
    st_arr = np.array(st)

    stat, p_val = stats.mannwhitneyu(sa_arr, st_arr, alternative="two-sided")

    pooled_std = math.sqrt(((len(sa_arr) - 1) * sa_arr.std(ddof=1)**2 +
                            (len(st_arr) - 1) * st_arr.std(ddof=1)**2) /
                           (len(sa_arr) + len(st_arr) - 2))
    d = (sa_arr.mean() - st_arr.mean()) / pooled_std if pooled_std > 0 else 0

    if p_val < 0.001:
        sig = "***"
    elif p_val < 0.01:
        sig = "**"
    elif p_val < 0.05:
        sig = "*"
    else:
        sig = "n.s."

    # For log-likelihood, higher is better, so positive d means SA is better
    if sig == "n.s.":
        interp = "Indistinguishable"
    elif d > 0:
        interp = f"SA BETTER ({sig})"
    else:
        interp = f"Stored better ({sig})"

    print(f"{fam:<16} {sa_arr.mean():8.3f}+/-{sa_arr.std():.3f} {st_arr.mean():8.3f}+/-{st_arr.std():.3f} {p_val:12.4e} {d:10.3f} {interp}")

# Overall summary
print(f"\n{'='*80}")
print("SUMMARY")
print(f"{'='*80}")
all_sa_ppl = [r["pseudo_perplexity"] for r in rows if r["source"] == "SA_generation"]
all_st_ppl = [r["pseudo_perplexity"] for r in rows if r["source"] == "Stored"]
print(f"Overall SA gen pseudo-perplexity: {np.mean(all_sa_ppl):.2f} +/- {np.std(all_sa_ppl):.2f} (n={len(all_sa_ppl)})")
print(f"Overall Stored pseudo-perplexity: {np.mean(all_st_ppl):.2f} +/- {np.std(all_st_ppl):.2f} (n={len(all_st_ppl)})")
