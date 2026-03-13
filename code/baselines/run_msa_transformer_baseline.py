#!/usr/bin/env python3
"""
MSA Transformer Baseline: Gibbs-like MLM sampling from an MSA context.

The MSA Transformer (Rao et al., 2021) is a masked language model trained on
protein MSAs. We use it for generation via iterative masking and infilling:
  1. Append a fully masked query row to the family MSA
  2. Forward pass → sample from softmax at each masked position
  3. Re-mask a fraction of positions, re-sample (Gibbs)
  4. After N_GIBBS iterations, the generated row is the output

Install:
    pip install fair-esm torch

Usage:
    python run_msa_transformer_baseline.py
    python run_msa_transformer_baseline.py --device cpu --n_gibbs 50
"""

import os
import sys
import time
import argparse
import warnings
import numpy as np

warnings.filterwarnings("ignore", category=UserWarning, message=".*pkg_resources.*")

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

N_SEQUENCES = 150
N_GIBBS_ITER = 50
MASK_FRACTION = 0.15
MAX_MSA_DEPTH = 128  # MSA Transformer max depth
PFAM_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pfam-families", "data")
BASELINE_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

AA_SET = set("ACDEFGHIKLMNPQRSTVWY")


def read_stockholm(filepath):
    """Parse a Stockholm alignment file."""
    sequences = {}
    seq_order = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or line.startswith("//") or not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            name, seq = parts[0], parts[1].upper()
            if name in sequences:
                sequences[name] += seq
            else:
                sequences[name] = seq
                seq_order.append(name)
    return [(name, sequences[name]) for name in seq_order]


def clean_alignment(raw_seqs, max_gap_col=0.5, max_gap_seq=0.3):
    """Clean alignment: remove high-gap columns and sequences."""
    seqs = [s for _, s in raw_seqs]
    names = [n for n, _ in raw_seqs]
    L = len(seqs[0])

    mat = np.array([[c for c in seq] for seq in seqs])

    is_gap = np.isin(mat, ['.', '-', '~'])
    col_gap_frac = is_gap.mean(axis=0)
    keep_cols = col_gap_frac <= max_gap_col
    mat = mat[:, keep_cols]

    is_gap = np.isin(mat, ['.', '-', '~'])
    seq_gap_frac = is_gap.mean(axis=1)
    keep_seqs = seq_gap_frac <= max_gap_seq
    mat = mat[keep_seqs]
    names = [n for n, keep in zip(names, keep_seqs) if keep]

    return mat, names


def write_fasta(sequences, filepath):
    """Write sequences as FASTA."""
    with open(filepath, "w") as f:
        for name, seq in sequences:
            f.write(f">{name}\n{seq}\n")


def main():
    parser = argparse.ArgumentParser(description="MSA Transformer MLM sampling baseline")
    parser.add_argument("--device", default="cpu", help="Device: cpu or cuda")
    parser.add_argument("--n_gibbs", type=int, default=N_GIBBS_ITER,
                        help=f"Gibbs iterations per sequence (default: {N_GIBBS_ITER})")
    parser.add_argument("--families", nargs="*", default=None,
                        help="Subset of family IDs to run (default: all)")
    args = parser.parse_args()

    try:
        import torch
        import esm
    except ImportError as e:
        print(f"ERROR: {e}")
        print("Install with: pip install fair-esm torch")
        sys.exit(1)

    device = args.device
    print(f"Loading MSA Transformer (device={device}) …")
    model, alphabet = esm.pretrained.esm_msa1b_t12_100M_UR50S()
    model = model.to(device)
    model.eval()
    batch_converter = alphabet.get_batch_converter()
    mask_idx = alphabet.mask_idx
    print(f"Model loaded. mask_idx={mask_idx}")

    # Build a reverse mapping from token index to character
    tok_to_char = {}
    for i in range(len(alphabet)):
        tok = alphabet.get_tok(i)
        tok_to_char[i] = tok

    families_to_run = FAMILIES
    if args.families:
        families_to_run = [(fid, fn) for fid, fn in FAMILIES if fid in args.families]

    for pfam_id, fam_name in families_to_run:
        print(f"\n{'='*60}")
        print(f"  {fam_name} ({pfam_id})")
        print(f"{'='*60}")

        sto_file = os.path.join(PFAM_DATA_DIR, pfam_id, f"{pfam_id}_seed.sto")
        if not os.path.isfile(sto_file):
            print(f"  WARNING: Missing {sto_file}, skipping")
            continue

        raw_seqs = read_stockholm(sto_file)
        char_mat, names = clean_alignment(raw_seqs)
        K, L = char_mat.shape
        print(f"  Alignment: {K} sequences × {L} positions")

        # Convert to list of strings (AAs + gaps as '-')
        msa_strings = []
        for i in range(K):
            seq = "".join(char_mat[i])
            seq = "".join(c if c in AA_SET else "-" for c in seq.upper())
            msa_strings.append(seq)

        # Subsample MSA if too deep
        rng = np.random.RandomState(42)
        if K > MAX_MSA_DEPTH - 1:
            idx = rng.choice(K, MAX_MSA_DEPTH - 1, replace=False)
            context_seqs = [msa_strings[i] for i in sorted(idx)]
        else:
            context_seqs = msa_strings.copy()

        print(f"  Context MSA: {len(context_seqs)} sequences")

        out_dir = os.path.join(BASELINE_DATA_DIR, pfam_id)
        os.makedirs(out_dir, exist_ok=True)

        # Check for existing output (resume support)
        out_fasta = os.path.join(out_dir, "msat_sequences.fasta")
        existing = []
        if os.path.isfile(out_fasta):
            with open(out_fasta) as f:
                for line in f:
                    if line.startswith(">"):
                        existing.append(line.strip())
            print(f"  Found {len(existing)} existing sequences, resuming …")

        generated = []
        t_start = time.time()
        start_idx = len(existing)

        import torch

        for seq_i in range(start_idx, N_SEQUENCES):
            if (seq_i + 1) % 5 == 0:
                elapsed = time.time() - t_start
                done = seq_i - start_idx + 1
                rate = done / elapsed if elapsed > 0 else 0
                remaining = N_SEQUENCES - seq_i - 1
                eta = remaining / rate if rate > 0 else 0
                print(f"  Sequence {seq_i+1}/{N_SEQUENCES} "
                      f"({elapsed:.0f}s, ~{eta:.0f}s remaining)")

            try:
                # Build MSA: query row (will be masked) + context sequences
                # ESM batch_converter expects equal-length sequences,
                # so use a dummy query of correct length, then mask after tokenization
                dummy_query = context_seqs[0]  # same length as alignment
                msa_input = [("query", dummy_query)]
                for j, seq in enumerate(context_seqs):
                    msa_input.append((f"ctx_{j}", seq))

                # Tokenize
                _, _, batch_tokens = batch_converter([msa_input])
                batch_tokens = batch_tokens.to(device)
                # shape: (1, n_seqs, L+2) — includes <cls> at 0 and <eos> at L+1

                # Mask the query row (positions 1..L, skip cls at 0)
                batch_tokens[0, 0, 1:L+1] = mask_idx

                # Initialize current query tokens
                current_query = batch_tokens[0, 0, :].clone()

                # Gibbs sampling
                for g_iter in range(args.n_gibbs):
                    # Set query row in batch
                    batch_tokens[0, 0, :] = current_query

                    if g_iter > 0:
                        # Re-mask a fraction of positions (positions 1..L)
                        n_mask = max(1, int(MASK_FRACTION * L))
                        mask_positions = rng.choice(L, n_mask, replace=False) + 1  # +1 for cls offset
                        for pos in mask_positions:
                            batch_tokens[0, 0, pos] = mask_idx
                    else:
                        # First iteration: all positions are already masked
                        mask_positions = np.arange(1, L + 1)

                    with torch.no_grad():
                        results = model(batch_tokens)
                        logits = results["logits"]  # (1, n_seqs, L+2, vocab)

                    # Sample from softmax at masked positions in query row
                    for pos in mask_positions:
                        pos_logits = logits[0, 0, pos, :]
                        probs = torch.softmax(pos_logits, dim=-1)
                        sampled_tok = torch.multinomial(probs, 1).item()
                        current_query[pos] = sampled_tok

                # Decode final query tokens to sequence
                gen_seq = ""
                for pos in range(1, L + 1):  # skip cls token
                    tok_idx = current_query[pos].item()
                    c = tok_to_char.get(tok_idx, "")
                    if c in AA_SET:
                        gen_seq += c
                    # skip gaps, special tokens

                if len(gen_seq) > 0:
                    generated.append((f"msat_{seq_i+1}", gen_seq))
                    # Write incrementally so we don't lose progress
                    mode = "a" if (existing or seq_i > start_idx) else "w"
                    with open(out_fasta, mode) as f:
                        f.write(f">msat_{seq_i+1}\n{gen_seq}\n")

            except Exception as e:
                print(f"  WARNING: Failed for seq {seq_i+1}: {e}")
                import traceback
                traceback.print_exc()
                continue

        t_elapsed = time.time() - t_start
        n_new = len(generated)
        total = len(existing) + n_new
        print(f"  Generated {n_new} new sequences in {t_elapsed:.1f}s")
        print(f"  Total: {total} sequences in {out_fasta}")

        if total == 0:
            print(f"  WARNING: No sequences generated for {fam_name}")

    print("\nDone. Run evaluate_baselines.jl to compute metrics.")


if __name__ == "__main__":
    main()
