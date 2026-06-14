#!/usr/bin/env python3
"""
Position-matched (and conservation-matched) DMS null for the SA substitution analysis.

Reviewer concern: the original null compares the SA-generated substitution
tolerance rate against the *overall* tolerance rate of all single mutations in
each DMS dataset. Because SA preferentially substitutes at variable (tolerant)
positions, that comparison conflates *where* SA mutates with *which* residue it
chooses. The enrichment could therefore be an artifact of position selection.

This script adds the appropriate controls:

  (1) POSITION-MATCHED null. For each position SA actually substituted, the null
      tolerance probability is the fraction of *all tested mutant amino acids at
      that same position* that are tolerated. Summing these per-event
      probabilities gives the expected number of tolerated SA substitutions if
      SA picked a random tested residue at each position it chose. This holds the
      position distribution fixed and isolates residue choice. A Poisson-binomial
      z-test and a Monte-Carlo permutation test give one-sided p-values.

  (2) CONSERVATION-MATCHED null (coarser, robustness check). Positions are binned
      by alignment conservation (max per-column amino-acid frequency in the stored
      seed alignment); each SA substitution is matched to a random tested DMS
      substitution drawn from positions in the same conservation quintile.

Reuses the exact matched-substitution extraction logic of analyze_dms.py so the
matched set is identical to the published table.
"""

import csv
import math
from collections import defaultdict, Counter
from pathlib import Path
import numpy as np


def norm_sf(z):
    """One-sided upper-tail probability of the standard normal (no scipy dependency)."""
    return 0.5 * math.erfc(z / math.sqrt(2.0))

BASE = Path(__file__).resolve().parent.parent
DMS_DATA = Path(__file__).resolve().parent / "data"
PFAM_DATA = BASE / "pfam-families" / "data"
OUT_DIR = Path(__file__).resolve().parent / "results"
OUT_DIR.mkdir(exist_ok=True)

RNG = np.random.default_rng(0)
N_PERM = 20000
N_CONS_BINS = 5

# ── Helpers (copied from analyze_dms.py so the matched set is identical) ─────────

def read_fasta(path):
    seqs, name, buf = [], None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(buf)))
                name, buf = line[1:], []
            else:
                buf.append(line)
    if name is not None:
        seqs.append((name, "".join(buf)))
    return seqs


def build_consensus(seqs):
    L = len(seqs[0][1])
    consensus = []
    for pos in range(L):
        counts = Counter()
        for _, seq in seqs:
            aa = seq[pos]
            if aa not in (".", "-"):
                counts[aa] += 1
        consensus.append(counts.most_common(1)[0][0] if counts else "-")
    return "".join(consensus)


def column_conservation(seqs):
    """Max per-column amino-acid frequency (ignoring gaps) for each alignment column."""
    L = len(seqs[0][1])
    cons = []
    for pos in range(L):
        counts = Counter()
        for _, seq in seqs:
            aa = seq[pos]
            if aa not in (".", "-"):
                counts[aa] += 1
        total = sum(counts.values())
        cons.append(counts.most_common(1)[0][1] / total if total else 0.0)
    return cons


def load_single_mutant_dms(csv_path):
    dms = defaultdict(dict)
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            m = row["mutant"]
            if ":" in m:
                continue
            pos = int(m[1:-1])
            mut_aa = m[-1]
            dms[pos][mut_aa] = (float(row["DMS_score"]), int(row["DMS_score_bin"]), m[0])
    return dms


def find_domain_in_stored(stored_seqs, protein_prefix):
    for name, seq in stored_seqs:
        if protein_prefix in name:
            parts = name.split("/")
            if len(parts) == 2:
                start, _ = map(int, parts[1].split("-"))
                aln_to_prot, prot_pos = {}, start
                for aln_pos, char in enumerate(seq):
                    if char not in (".", "-"):
                        aln_to_prot[aln_pos] = prot_pos
                        prot_pos += 1
                return name, seq, aln_to_prot
    return None, None, None


FAMILIES = [
    {"name": "PDZ", "pfam_id": "PF00595", "dms_file": "DLG4_RAT_McLaughlin_2012.csv",
     "protein_prefix": "DLG4_RAT/313"},
    {"name": "SH3", "pfam_id": "PF00018", "dms_file": "GRB2_HUMAN_Faure_2021.csv",
     "protein_prefix": "GRB2_CHICK"},
]

all_results = []

for fam in FAMILIES:
    print(f"\n{'='*70}\nFamily: {fam['name']} ({fam['pfam_id']})\n{'='*70}")

    dms = load_single_mutant_dms(DMS_DATA / fam["dms_file"])

    # Per-position tolerance rate over all tested mutant AAs (the position-matched null prob)
    pos_tol_rate = {p: np.mean([v[1] for v in d.values()]) for p, d in dms.items()}
    pos_tested_bins = {p: [v[1] for v in d.values()] for p, d in dms.items()}

    all_bins = [v[1] for d in dms.values() for v in d.values()]
    overall_tol_rate = float(np.mean(all_bins))

    stored = read_fasta(PFAM_DATA / fam["pfam_id"] / "stored_sequences.fasta")
    consensus = build_consensus(stored)
    cons_by_col = column_conservation(stored)

    _, _, aln_to_prot = find_domain_in_stored(stored, fam["protein_prefix"])
    aln_to_dms = {a: p for a, p in aln_to_prot.items() if p in dms}

    sa_seqs = read_fasta(PFAM_DATA / fam["pfam_id"] / "sa_generation_sequences.fasta")

    # Collect each matched SA substitution event: (dms_pos, conservation, observed_bin)
    events = []   # (dms_pos, observed_bin)
    cons_of_event = []
    for _, seq in sa_seqs:
        clean = seq.replace(".", "").replace("-", "")
        seq_idx = 0
        for aln_pos in range(len(consensus)):
            cons_aa = consensus[aln_pos]
            if cons_aa in (".", "-"):
                continue
            if seq_idx >= len(clean):
                break
            gen_aa = clean[seq_idx]
            seq_idx += 1
            if gen_aa != cons_aa and aln_pos in aln_to_dms:
                dms_pos = aln_to_dms[aln_pos]
                if gen_aa in dms[dms_pos]:
                    events.append((dms_pos, dms[dms_pos][gen_aa][1]))
                    cons_of_event.append(cons_by_col[aln_pos])

    n = len(events)
    obs_tol = sum(b for _, b in events)
    obs_rate = obs_tol / n

    # ── Position-matched null ──────────────────────────────────────────────────
    p_i = np.array([pos_tol_rate[p] for p, _ in events])      # per-event null prob
    exp_tol = p_i.sum()
    exp_rate = float(p_i.mean())
    var = float((p_i * (1 - p_i)).sum())                      # Poisson-binomial variance
    z = (obs_tol - exp_tol) / np.sqrt(var) if var > 0 else np.nan
    z_p = norm_sf(z)                                           # one-sided (SA > null)

    # Monte-Carlo permutation: random tested AA at each chosen position. Drawing a
    # uniform tested residue at a position and reading its bin is Bernoulli(p_i) with
    # p_i the position tolerance rate, so the null count is Poisson-binomial(p_i) and
    # can be sampled vectorized exactly.
    perm_counts = (RNG.random((N_PERM, n)) < p_i[None, :]).sum(axis=1)
    perm_rate = perm_counts / n
    perm_p = (np.sum(perm_counts >= obs_tol) + 1) / (N_PERM + 1)
    enr_pos = obs_rate / exp_rate

    # ── Conservation-matched null (quintiles of column conservation) ────────────
    # Build pool of tested bins per conservation bin, using the DMS-covered columns.
    covered_cons = [cons_by_col[a] for a in aln_to_dms]
    edges = np.quantile(covered_cons, np.linspace(0, 1, N_CONS_BINS + 1))
    edges[0] -= 1e-9
    def cbin(c):
        return int(np.clip(np.searchsorted(edges, c, side="left") - 1, 0, N_CONS_BINS - 1))
    pool = defaultdict(list)
    for a, p in aln_to_dms.items():
        cb = cbin(cons_by_col[a])
        pool[cb].extend(pos_tested_bins[p])
    pool = {k: np.array(v) for k, v in pool.items()}
    event_cbins = [cbin(c) for c in cons_of_event]
    cons_p_i = np.array([pool[cb].mean() for cb in event_cbins])  # per-event null prob
    cons_perm = (RNG.random((N_PERM, n)) < cons_p_i[None, :]).sum(axis=1)
    cons_exp_rate = float(cons_p_i.mean())
    cons_perm_p = (np.sum(cons_perm >= obs_tol) + 1) / (N_PERM + 1)
    enr_cons = obs_rate / cons_exp_rate

    print(f"  matched substitutions n           = {n}")
    print(f"  observed SA tolerance rate         = {obs_rate:.4f}")
    print(f"  overall DMS null rate              = {overall_tol_rate:.4f}  (enrichment {obs_rate/overall_tol_rate:.3f}x)")
    print(f"  POSITION-matched null rate         = {exp_rate:.4f}  (enrichment {enr_pos:.3f}x)")
    print(f"    Poisson-binomial z = {z:.2f}, p = {z_p:.3e}")
    print(f"    permutation p ({N_PERM} draws) = {perm_p:.3e}  [null rate {perm_rate.mean():.4f} +/- {perm_rate.std():.4f}]")
    print(f"  CONSERVATION-matched null rate     = {cons_exp_rate:.4f}  (enrichment {enr_cons:.3f}x)")
    print(f"    permutation p = {cons_perm_p:.3e}  [null rate {(cons_perm/n).mean():.4f} +/- {(cons_perm/n).std():.4f}]")

    all_results.append({
        "family": fam["name"], "pfam_id": fam["pfam_id"], "n_matched": n,
        "obs_tol_rate": round(obs_rate, 4),
        "overall_null_rate": round(overall_tol_rate, 4),
        "overall_enrichment": round(obs_rate / overall_tol_rate, 3),
        "posmatched_null_rate": round(exp_rate, 4),
        "posmatched_enrichment": round(enr_pos, 3),
        "posmatched_z": round(float(z), 3),
        "posmatched_z_p": z_p,
        "posmatched_perm_p": perm_p,
        "consmatched_null_rate": round(cons_exp_rate, 4),
        "consmatched_enrichment": round(enr_cons, 3),
        "consmatched_perm_p": cons_perm_p,
    })

out_path = OUT_DIR / "dms_position_matched_null.csv"
with open(out_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=all_results[0].keys())
    w.writeheader()
    w.writerows(all_results)
print(f"\nSaved {out_path}")
