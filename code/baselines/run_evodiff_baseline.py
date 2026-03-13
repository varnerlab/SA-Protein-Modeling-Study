#!/usr/bin/env python3
"""
EvoDiff Baseline: MSA-conditioned protein sequence generation.

Uses the EvoDiff MSA_OA_DM_RANDSUB model (order-agnostic discrete diffusion,
MSA-conditioned) to generate sequences from family alignments.

Install:
    pip install evodiff torch "setuptools<81"

Usage:
    python run_evodiff_baseline.py
    python run_evodiff_baseline.py --families PF00076 PF00018

References:
    Alamdari et al. (2023) "Protein generation with evolutionary diffusion"
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


def write_aligned_fasta(char_mat, names, filepath):
    """Write a cleaned alignment as FASTA (uppercase AAs, '-' for gaps)."""
    with open(filepath, "w") as f:
        for i, name in enumerate(names):
            seq = "".join(char_mat[i])
            # standardize: uppercase AAs, '-' for anything else
            seq = "".join(c if c in AA_SET else "-" for c in seq.upper())
            f.write(f">{name}\n{seq}\n")


def write_fasta(sequences, filepath):
    """Write sequences as FASTA."""
    with open(filepath, "w") as f:
        for name, seq in sequences:
            f.write(f">{name}\n{seq}\n")


def main():
    parser = argparse.ArgumentParser(description="EvoDiff MSA-conditioned baseline")
    parser.add_argument("--device", default="cpu", help="Device: cpu or cuda")
    parser.add_argument("--families", nargs="*", default=None,
                        help="Subset of family IDs to run (default: all)")
    args = parser.parse_args()

    try:
        import torch
        from evodiff.pretrained import MSA_OA_DM_RANDSUB
        from evodiff.generate_msa import generate_query_oadm_msa_simple
    except ImportError as e:
        print(f"ERROR: {e}")
        print("Install with: pip install evodiff torch 'setuptools<81'")
        sys.exit(1)

    device = args.device
    print(f"Loading EvoDiff MSA_OA_DM_RANDSUB model (device={device}) …")
    checkpoint = MSA_OA_DM_RANDSUB()
    model, collater, tokenizer, scheme = checkpoint
    model = model.to(device)
    model.eval()
    print("Model loaded.")

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

        # Write cleaned alignment as FASTA for EvoDiff
        out_dir = os.path.join(BASELINE_DATA_DIR, pfam_id)
        os.makedirs(out_dir, exist_ok=True)
        aligned_fasta = os.path.join(out_dir, "aligned_for_evodiff.fasta")
        write_aligned_fasta(char_mat, names, aligned_fasta)

        # EvoDiff requires n_sequences (MSA depth) and seq_length
        # For small families, n_sequences = K (use all); model pads if needed
        n_msa_seqs = min(K, 64)  # EvoDiff default max depth

        # Check if MSA is too small for EvoDiff
        if K < 2:
            print(f"  WARNING: MSA too small (K={K}), skipping")
            continue

        generated = []
        t_start = time.time()

        # Generate sequences one at a time
        for i in range(N_SEQUENCES):
            if (i + 1) % 10 == 0:
                elapsed = time.time() - t_start
                rate = (i + 1) / elapsed if elapsed > 0 else 0
                eta = (N_SEQUENCES - i - 1) / rate if rate > 0 else 0
                print(f"  Sequence {i+1}/{N_SEQUENCES} "
                      f"({elapsed:.0f}s elapsed, ~{eta:.0f}s remaining)")

            try:
                # generate_query_oadm_msa_simple returns (generated_samples, query_sequences)
                gen_result = generate_query_oadm_msa_simple(
                    path_to_msa=aligned_fasta,
                    model=model,
                    tokenizer=tokenizer,
                    n_sequences=n_msa_seqs,
                    seq_length=L,
                    batch_size=1,
                    device=device,
                    selection_type='MaxHamming',
                )
                # gen_result is (sample_tensor, query_list)
                # sample tensor shape: (1, n_sequences, seq_length)
                sample_tensor = gen_result[0]
                # Row 0 is the generated query sequence
                gen_tokens = sample_tensor[0, 0, :].cpu().numpy()

                # Decode tokens back to amino acids
                gen_seq = ""
                for tok in gen_tokens:
                    tok = int(tok)
                    if tok < len(tokenizer.alphabet):
                        c = tokenizer.alphabet[tok]
                        if c in AA_SET:
                            gen_seq += c
                    # skip gaps, pads, special tokens

                if len(gen_seq) > 0:
                    generated.append((f"evodiff_{i+1}", gen_seq))

            except Exception as e:
                print(f"  WARNING: Generation failed for seq {i+1}: {e}")
                # If MSA is too small, EvoDiff may raise an exception
                if "n_sequences" in str(e).lower() or "msa" in str(e).lower():
                    print(f"  MSA too small for EvoDiff (K={K}). Stopping this family.")
                    break
                continue

        t_elapsed = time.time() - t_start
        print(f"  Generated {len(generated)} sequences in {t_elapsed:.1f}s")

        if generated:
            out_fasta = os.path.join(out_dir, "evodiff_sequences.fasta")
            write_fasta(generated, out_fasta)
            print(f"  Saved: {out_fasta}")
        else:
            print(f"  WARNING: No sequences generated for {fam_name}")

    print("\nDone. Run evaluate_baselines.jl to compute metrics.")


if __name__ == "__main__":
    main()
