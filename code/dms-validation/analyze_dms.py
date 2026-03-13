#!/usr/bin/env python3
"""
DMS validation: check whether SA-generated substitutions are enriched
at experimentally tolerated positions vs. deleterious ones.

For each family with DMS data (PDZ, SH3, RRM):
1. Extract single-mutant DMS scores (position -> {AA: score, bin})
2. Build consensus from stored alignment
3. For each SA-generated sequence, identify substitutions vs consensus
4. Look up each substitution in the DMS data
5. Compare: fraction of SA substitutions that are tolerated (DMS_score_bin=1)
   vs. the overall fraction of tolerated mutations in the DMS dataset
"""

import csv
import os
import sys
from collections import defaultdict, Counter
from pathlib import Path
import numpy as np
from scipy import stats

BASE = Path(__file__).resolve().parent.parent
DMS_DATA = Path(__file__).resolve().parent / "data"
PFAM_DATA = BASE / "pfam-families" / "data"
OUT_DIR = Path(__file__).resolve().parent / "results"
OUT_DIR.mkdir(exist_ok=True)

# ── Helpers ────────────────────────────────────────────────────────────────────

def read_fasta(path):
    seqs = []
    name, buf = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(buf)))
                name = line[1:]
                buf = []
            else:
                buf.append(line)
    if name is not None:
        seqs.append((name, "".join(buf)))
    return seqs


def build_consensus(seqs):
    """Build position-wise consensus (most frequent AA, ignoring gaps)."""
    L = len(seqs[0][1])
    consensus = []
    for pos in range(L):
        counts = Counter()
        for _, seq in seqs:
            aa = seq[pos]
            if aa not in (".", "-"):
                counts[aa] += 1
        if counts:
            consensus.append(counts.most_common(1)[0][0])
        else:
            consensus.append("-")
    return "".join(consensus)


def load_single_mutant_dms(csv_path):
    """Load single-mutant DMS scores. Returns dict: pos -> {mut_aa: (score, bin)}"""
    dms = defaultdict(dict)
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            m = row["mutant"]
            if ":" in m:
                continue  # skip multi-mutants
            wt_aa = m[0]
            pos = int(m[1:-1])
            mut_aa = m[-1]
            score = float(row["DMS_score"])
            score_bin = int(row["DMS_score_bin"])
            dms[pos][mut_aa] = (score, score_bin, wt_aa)
    return dms


def find_domain_in_stored(stored_seqs, protein_prefix, full_protein_seq, dms_start, dms_end):
    """
    Find the stored sequence matching the DMS protein and return
    the mapping from alignment position -> DMS position.
    """
    # Find the stored sequence header matching the protein
    for name, seq in stored_seqs:
        if protein_prefix in name:
            # Parse the range from header, e.g., "DLG4_RAT/313-391"
            parts = name.split("/")
            if len(parts) == 2:
                range_str = parts[1]
                start, end = map(int, range_str.split("-"))
                # This stored sequence spans positions start-end in the full protein
                # The alignment may have gaps (dots/dashes)
                clean_seq = seq.replace(".", "").replace("-", "")
                print(f"  Found {name}: protein positions {start}-{end}, "
                      f"aligned length {len(seq)}, clean length {len(clean_seq)}")

                # Build alignment_pos -> protein_pos mapping
                # Walk through the aligned sequence, incrementing protein_pos
                # only for non-gap characters
                aln_to_prot = {}
                prot_pos = start
                for aln_pos, char in enumerate(seq):
                    if char not in (".", "-"):
                        aln_to_prot[aln_pos] = prot_pos
                        prot_pos += 1

                return name, seq, aln_to_prot
    return None, None, None


# ── Analysis per family ────────────────────────────────────────────────────────

FAMILIES = [
    {
        "name": "PDZ",
        "pfam_id": "PF00595",
        "dms_file": "DLG4_RAT_McLaughlin_2012.csv",
        "protein_prefix": "DLG4_RAT/313",  # matches the 313-391 range
        "dms_start": 311,
        "dms_end": 393,
    },
    {
        "name": "SH3",
        "pfam_id": "PF00018",
        "dms_file": "GRB2_HUMAN_Faure_2021.csv",
        "protein_prefix": "GRB2_CHICK",  # closest homolog in our alignment
        "dms_start": 159,
        "dms_end": 214,
    },
    # RRM excluded: PABP_YEAST DMS (positions 126-200) covers a different RRM domain
    # than the PABP sequences in our Pfam seed alignment (positions 263-332, 92-161, etc.)
    # No direct position mapping is possible without full sequence alignment.
]

all_results = []

for fam in FAMILIES:
    print(f"\n{'='*70}")
    print(f"Family: {fam['name']} ({fam['pfam_id']})")
    print(f"DMS: {fam['dms_file']}")
    print(f"{'='*70}")

    # Load DMS data
    dms = load_single_mutant_dms(DMS_DATA / fam["dms_file"])
    dms_positions = sorted(dms.keys())
    print(f"  DMS: {len(dms)} positions, {sum(len(v) for v in dms.values())} single mutations")

    # Overall DMS stats
    all_bins = [v[1] for pos_dict in dms.values() for v in pos_dict.values()]
    overall_tol_rate = sum(all_bins) / len(all_bins)
    print(f"  Overall DMS tolerance rate: {overall_tol_rate:.3f} ({sum(all_bins)}/{len(all_bins)})")

    # Load stored sequences and build consensus
    stored = read_fasta(PFAM_DATA / fam["pfam_id"] / "stored_sequences.fasta")
    consensus = build_consensus(stored)
    clean_consensus = consensus.replace(".", "").replace("-", "")
    print(f"  Stored: {len(stored)} sequences, alignment L={len(consensus)}, clean L={len(clean_consensus)}")

    # Find the DMS protein in stored sequences and get position mapping
    ref_name, ref_seq, aln_to_prot = find_domain_in_stored(
        stored, fam["protein_prefix"], None, fam["dms_start"], fam["dms_end"]
    )
    if aln_to_prot is None:
        print(f"  WARNING: Could not find {fam['protein_prefix']} in stored sequences!")
        continue

    # Build reverse mapping: for each alignment position, what is the DMS protein position?
    # Only positions covered by the DMS data
    aln_to_dms = {}
    for aln_pos, prot_pos in aln_to_prot.items():
        if prot_pos in dms:
            aln_to_dms[aln_pos] = prot_pos
    print(f"  Mapped {len(aln_to_dms)} alignment positions to DMS positions")

    # Load SA-generated sequences
    sa_seqs = read_fasta(PFAM_DATA / fam["pfam_id"] / "sa_generation_sequences.fasta")
    print(f"  SA-generated: {len(sa_seqs)} sequences")

    # For each SA-generated sequence, find substitutions vs consensus
    # and look up DMS scores
    sa_tolerated = 0
    sa_deleterious = 0
    sa_not_in_dms = 0  # substitution at a position covered by DMS but AA not tested
    sa_no_dms_pos = 0  # substitution at a position not in DMS
    sa_total_subs = 0
    sa_scores = []

    for seq_name, seq in sa_seqs:
        clean_seq = seq.replace(".", "").replace("-", "")
        # Walk alignment positions
        seq_idx = 0
        for aln_pos in range(len(consensus)):
            cons_aa = consensus[aln_pos]
            if cons_aa in (".", "-"):
                continue

            if seq_idx >= len(clean_seq):
                break

            gen_aa = clean_seq[seq_idx]
            seq_idx += 1

            if gen_aa != cons_aa:
                sa_total_subs += 1
                if aln_pos in aln_to_dms:
                    dms_pos = aln_to_dms[aln_pos]
                    if gen_aa in dms[dms_pos]:
                        score, score_bin, wt_aa = dms[dms_pos][gen_aa]
                        sa_scores.append(score)
                        if score_bin == 1:
                            sa_tolerated += 1
                        else:
                            sa_deleterious += 1
                    else:
                        sa_not_in_dms += 1
                else:
                    sa_no_dms_pos += 1

    sa_matched = sa_tolerated + sa_deleterious
    sa_tol_rate = sa_tolerated / sa_matched if sa_matched > 0 else 0

    print(f"\n  SA-generated substitution analysis:")
    print(f"    Total substitutions vs consensus: {sa_total_subs}")
    print(f"    Matched to DMS data: {sa_matched}")
    print(f"      Tolerated (DMS_bin=1): {sa_tolerated} ({sa_tol_rate:.1%})")
    print(f"      Deleterious (DMS_bin=0): {sa_deleterious} ({1-sa_tol_rate:.1%})")
    print(f"    Not in DMS (AA not tested): {sa_not_in_dms}")
    print(f"    Position not in DMS: {sa_no_dms_pos}")
    print(f"    Overall DMS tolerance rate: {overall_tol_rate:.1%}")
    print(f"    SA enrichment ratio: {sa_tol_rate/overall_tol_rate:.2f}x")

    # Statistical test: binomial test
    # H0: SA substitutions are drawn from the overall DMS distribution
    if sa_matched == 0:
        print("    No DMS matches found - skipping statistical tests")
        binom_p = None
    else:
        binom_result = stats.binomtest(sa_tolerated, sa_matched, overall_tol_rate,
                                        alternative="greater")
        binom_p = binom_result.pvalue
        print(f"    Binomial test (SA tol rate > overall): p = {binom_p:.4e}")

    # Also compare DMS scores of SA substitutions vs all DMS mutations
    all_scores = [v[0] for pos_dict in dms.values() for v in pos_dict.values()]
    mw_p = None
    d = None
    if sa_scores:
        stat, mw_p = stats.mannwhitneyu(sa_scores, all_scores, alternative="greater")
        d = (np.mean(sa_scores) - np.mean(all_scores)) / np.sqrt(
            ((len(sa_scores)-1)*np.std(sa_scores, ddof=1)**2 +
             (len(all_scores)-1)*np.std(all_scores, ddof=1)**2) /
            (len(sa_scores) + len(all_scores) - 2)
        )
        print(f"    SA DMS score: {np.mean(sa_scores):.3f} +/- {np.std(sa_scores):.3f}")
        print(f"    All DMS score: {np.mean(all_scores):.3f} +/- {np.std(all_scores):.3f}")
        print(f"    Mann-Whitney (SA > All): p = {mw_p:.4e}, Cohen's d = {d:.3f}")

    all_results.append({
        "family": fam["name"],
        "pfam_id": fam["pfam_id"],
        "n_sa_subs": sa_total_subs,
        "n_matched": sa_matched,
        "n_tolerated": sa_tolerated,
        "n_deleterious": sa_deleterious,
        "sa_tol_rate": sa_tol_rate,
        "overall_tol_rate": overall_tol_rate,
        "enrichment": sa_tol_rate / overall_tol_rate if overall_tol_rate > 0 else 0,
        "binom_p": binom_p,
        "sa_mean_dms": np.mean(sa_scores) if sa_scores else None,
        "all_mean_dms": np.mean(all_scores),
        "mw_p": mw_p if sa_scores else None,
        "cohens_d": d if sa_scores else None,
    })

# Save results
out_path = OUT_DIR / "dms_validation_results.csv"
with open(out_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=all_results[0].keys())
    writer.writeheader()
    writer.writerows(all_results)
print(f"\nResults saved to {out_path}")

# Summary table
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")
print(f"{'Family':<8} {'Matched':>8} {'SA tol%':>8} {'DMS tol%':>9} {'Enrich':>7} {'Binom p':>12}")
print("-" * 60)
for r in all_results:
    print(f"{r['family']:<8} {r['n_matched']:>8} {r['sa_tol_rate']:>7.1%} {r['overall_tol_rate']:>8.1%} {r['enrichment']:>6.2f}x {r['binom_p']:>12.2e}")
