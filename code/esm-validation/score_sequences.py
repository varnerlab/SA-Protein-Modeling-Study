#!/usr/bin/env python3
"""
ESM2 pseudo-perplexity scoring of SA-generated vs stored (natural) sequences.

For each family, scores:
  - sa_generation_sequences.fasta  (150 SA-generated)
  - stored_sequences.fasta         (natural family members)

Outputs per-sequence scores to CSV and prints summary statistics.
"""

import os, sys, csv, math, time
from pathlib import Path

import torch
import esm

# ── Configuration ──────────────────────────────────────────────────────────────
FAMILIES = [
    ("PF00076", "RRM"),
    ("PF00018", "SH3"),
    ("PF00397", "WW"),
    ("PF00014", "Kunitz"),
    ("PF00096", "zf-C2H2"),
    ("PF00595", "PDZ"),
    ("PF00069", "Pkinase"),
    ("PF00711", "Defensin_beta"),
]

DATA_ROOT = Path(__file__).resolve().parent.parent / "pfam-families" / "data"
OUT_DIR = Path(__file__).resolve().parent / "results"
OUT_DIR.mkdir(exist_ok=True)

DEVICE = "mps" if torch.backends.mps.is_available() else "cpu"
BATCH_SIZE = 8  # sequences per batch for pseudo-ppl


# ── Helpers ────────────────────────────────────────────────────────────────────
def read_fasta(path):
    """Read FASTA, strip gaps (dots and dashes), return list of (name, seq)."""
    seqs = []
    name, buf = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(buf).replace(".", "").replace("-", "")))
                name = line[1:]
                buf = []
            else:
                buf.append(line)
    if name is not None:
        seqs.append((name, "".join(buf).replace(".", "").replace("-", "")))
    return seqs


def pseudo_perplexity(model, alphabet, sequences, device, batch_size=8):
    """
    Compute pseudo-perplexity for each sequence using masked marginal scoring.

    For each position i, mask it, run a forward pass, and record
    log P(true_aa | context). Sum over positions to get log-likelihood.
    Pseudo-perplexity = exp(-mean_log_prob).

    Returns list of dicts with per-sequence scores.
    """
    batch_converter = alphabet.get_batch_converter()
    results = []

    for idx, (name, seq) in enumerate(sequences):
        L = len(seq)
        # Build masked versions: one per position
        # Process in batches to avoid OOM
        log_probs = []

        for start in range(0, L, batch_size):
            end = min(start + batch_size, L)
            batch_data = []
            for pos in range(start, end):
                masked_seq = seq[:pos] + "<mask>" + seq[pos + 1:]
                batch_data.append((f"pos_{pos}", masked_seq))

            _, _, batch_tokens = batch_converter(batch_data)
            batch_tokens = batch_tokens.to(device)

            with torch.no_grad():
                logits = model(batch_tokens)["logits"]  # (B, L+2, vocab)

            for i, pos in enumerate(start + j for j in range(end - start)):
                # token index is pos+1 (accounting for BOS token)
                token_idx = alphabet.get_idx(seq[pos])
                log_p = torch.log_softmax(logits[i, pos + 1, :], dim=-1)
                log_probs.append(log_p[token_idx].item())

        mean_log_p = sum(log_probs) / L
        ppl = math.exp(-mean_log_p)
        total_log_p = sum(log_probs)

        results.append({
            "name": name,
            "length": L,
            "total_log_likelihood": total_log_p,
            "mean_log_likelihood": mean_log_p,
            "pseudo_perplexity": ppl,
        })

        if (idx + 1) % 10 == 0 or idx == 0:
            print(f"    [{idx+1}/{len(sequences)}] {name[:40]:40s}  ppl={ppl:.2f}  mean_llk={mean_log_p:.3f}")

    return results


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    print(f"Device: {DEVICE}")
    print("Loading ESM2-650M...")
    t0 = time.time()
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.eval().to(DEVICE)
    print(f"Model loaded in {time.time()-t0:.1f}s")

    all_rows = []

    for pfam_id, family_name in FAMILIES:
        print(f"\n{'='*60}")
        print(f"Family: {family_name} ({pfam_id})")
        print(f"{'='*60}")

        family_dir = DATA_ROOT / pfam_id

        # Score SA-generated sequences
        gen_path = family_dir / "sa_generation_sequences.fasta"
        stored_path = family_dir / "stored_sequences.fasta"

        for source, fasta_path in [("SA_generation", gen_path), ("Stored", stored_path)]:
            if not fasta_path.exists():
                print(f"  SKIP {source}: {fasta_path} not found")
                continue

            seqs = read_fasta(fasta_path)
            # For stored sequences, cap at 50 to match structure validation
            if source == "Stored" and len(seqs) > 50:
                seqs = seqs[:50]

            print(f"\n  Scoring {source}: {len(seqs)} sequences (L={len(seqs[0][1])})")
            t1 = time.time()
            scores = pseudo_perplexity(model, alphabet, seqs, DEVICE, BATCH_SIZE)
            elapsed = time.time() - t1
            print(f"  Done in {elapsed:.1f}s ({elapsed/len(seqs):.2f}s/seq)")

            # Summary stats
            ppls = [s["pseudo_perplexity"] for s in scores]
            mlls = [s["mean_log_likelihood"] for s in scores]
            print(f"  Pseudo-perplexity: {sum(ppls)/len(ppls):.2f} +/- {(sum((p-sum(ppls)/len(ppls))**2 for p in ppls)/len(ppls))**0.5:.2f}")
            print(f"  Mean log-likelihood: {sum(mlls)/len(mlls):.3f} +/- {(sum((m-sum(mlls)/len(mlls))**2 for m in mlls)/len(mlls))**0.5:.3f}")

            for s in scores:
                all_rows.append({
                    "family": family_name,
                    "pfam_id": pfam_id,
                    "source": source,
                    "name": s["name"],
                    "length": s["length"],
                    "total_log_likelihood": s["total_log_likelihood"],
                    "mean_log_likelihood": s["mean_log_likelihood"],
                    "pseudo_perplexity": s["pseudo_perplexity"],
                })

    # Save all results
    out_path = OUT_DIR / "esm2_pseudoppl_scores.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_rows[0].keys())
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"\nAll results saved to {out_path}")
    print(f"Total sequences scored: {len(all_rows)}")


if __name__ == "__main__":
    main()
