#!/usr/bin/env python3
"""
Potts model (plmDCA) baseline: fit pseudo-likelihood DCA and generate sequences
via Gibbs sampling.

Implements pseudo-likelihood maximization for Potts model parameters (fields h
and couplings J), then samples sequences using Gibbs sweeps.

Reference: Ekeberg et al. (2013), Trinquier et al. (2021)
"""

import numpy as np
from scipy.optimize import minimize
from pathlib import Path
import time

# ── Configuration ──────────────────────────────────────────────────────────
FAMILIES = [
    {"id": "PF00076", "name": "RRM"},
    {"id": "PF00018", "name": "SH3"},
    {"id": "PF00397", "name": "WW"},
    {"id": "PF00014", "name": "Kunitz"},
    {"id": "PF00096", "name": "zf-C2H2"},
    {"id": "PF00595", "name": "PDZ"},
    {"id": "PF00069", "name": "Pkinase"},
    {"id": "PF00711", "name": "Defensin_beta"},
]

N_GEN = 150          # sequences to generate (match SA output)
LAMBDA_J = 0.01      # L2 regularization on couplings (standard plmDCA value)
LAMBDA_H = 0.01      # L2 regularization on fields
GIBBS_SWEEPS = 500   # number of full Gibbs sweeps per sequence
N_CHAINS = 150       # one chain per output sequence

# Amino acid alphabet (20 standard + gap)
ALPHABET = "ACDEFGHIKLMNPQRSTVWY-"
Q = len(ALPHABET)    # 21
AA_TO_IDX = {aa: i for i, aa in enumerate(ALPHABET)}

BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / "data"


def read_aligned_fasta(fasta_path):
    """Read aligned FASTA file, return list of sequences (strings)."""
    sequences = []
    current_seq = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_seq:
                    sequences.append("".join(current_seq))
                current_seq = []
            else:
                current_seq.append(line)
    if current_seq:
        sequences.append("".join(current_seq))
    return sequences


def encode_alignment(sequences):
    """Convert alignment to integer matrix (K x L)."""
    K = len(sequences)
    L = len(sequences[0])
    msa = np.zeros((K, L), dtype=np.int32)
    gap_idx = AA_TO_IDX["-"]
    for k, seq in enumerate(sequences):
        for i, c in enumerate(seq):
            msa[k, i] = AA_TO_IDX.get(c, gap_idx)
    return msa


def compute_weights(msa, threshold=0.8):
    """Compute sequence weights (inverse of number of similar sequences)."""
    K, L = msa.shape
    weights = np.ones(K)
    for k1 in range(K):
        for k2 in range(k1 + 1, K):
            identity = np.mean(msa[k1] == msa[k2])
            if identity >= threshold:
                weights[k1] += 1
                weights[k2] += 1
    weights = 1.0 / weights
    return weights / weights.sum()


def fit_plmdca(msa, weights, lambda_J=0.01, lambda_h=0.01):
    """Fit Potts model via pseudo-likelihood maximization.

    For each position i, optimize:
        max sum_k w_k log P(x_i^k | x_{-i}^k; h_i, J_i) - regularization

    Parameters per position i:
        h_i: (q,) local fields
        J_i: (q, n_other, q) couplings, where J_i[a, j_idx, b] is the
             coupling between amino acid a at position i and amino acid b
             at position other_positions[j_idx].

    Logits: logits[k, a] = h_i[a] + sum_{j!=i} J_i[a, j_idx, x_j^k]

    Returns:
        h: (L, q) local fields
        J: (L, L, q, q) couplings (symmetrized)
    """
    K, L = msa.shape
    q = Q

    h = np.zeros((L, q))
    J = np.zeros((L, L, q, q))

    for i in range(L):
        if i % 5 == 0:
            print(f"  Fitting position {i+1}/{L}...", flush=True)

        target = msa[:, i]  # (K,)
        other_positions = [j for j in range(L) if j != i]
        n_other = len(other_positions)

        # One-hot features for all other positions: (K, n_other * q)
        features = np.zeros((K, n_other * q))
        for idx, j in enumerate(other_positions):
            for k in range(K):
                features[k, idx * q + msa[k, j]] = 1.0

        # Total parameters: q (fields) + q * n_other * q (couplings)
        n_params = q + q * n_other * q

        def objective_and_grad(params):
            hi = params[:q]
            Ji = params[q:].reshape(q, n_other * q)  # (q, n_other*q)

            # logits[k, a] = h_i[a] + features[k] @ Ji[a]
            logits = hi[np.newaxis, :] + features @ Ji.T  # (K, q)

            # Stable softmax
            logits -= logits.max(axis=1, keepdims=True)
            exp_logits = np.exp(logits)
            probs = exp_logits / exp_logits.sum(axis=1, keepdims=True)

            # Weighted negative log-likelihood
            log_probs = np.log(probs[np.arange(K), target] + 1e-30)
            nll = -np.sum(weights * log_probs)
            reg = (lambda_h / 2) * np.sum(hi**2) + (lambda_J / 2) * np.sum(params[q:]**2)
            loss = nll + reg

            # Gradient
            indicator = np.zeros((K, q))
            indicator[np.arange(K), target] = 1.0
            w_residual = (probs - indicator) * weights[:, np.newaxis]  # (K, q)

            grad_h = w_residual.sum(axis=0) + lambda_h * hi
            grad_Ji = w_residual.T @ features + lambda_J * Ji  # (q, n_other*q)

            grad = np.concatenate([grad_h, grad_Ji.ravel()])
            return loss, grad

        x0 = np.zeros(n_params)
        result = minimize(
            objective_and_grad,
            x0,
            jac=True,
            method="L-BFGS-B",
            options={"maxiter": 200, "ftol": 1e-6}
        )

        hi_fit = result.x[:q]
        Ji_fit = result.x[q:].reshape(q, n_other, q)

        h[i] = hi_fit
        for idx, j in enumerate(other_positions):
            J[i, j] = Ji_fit[:, idx, :]

    # Symmetrize: J_ij(a,b) = (J_ij(a,b) + J_ji(b,a)) / 2
    for i in range(L):
        for j in range(i + 1, L):
            J_sym = (J[i, j] + J[j, i].T) / 2
            J[i, j] = J_sym
            J[j, i] = J_sym.T

    return h, J


def gibbs_sample(h, J, msa, n_sequences=150, n_sweeps=500, rng=None):
    """Generate sequences by Gibbs sampling from fitted Potts model.

    Each chain is initialized from a random stored sequence.
    """
    if rng is None:
        rng = np.random.default_rng(42)

    L, q = h.shape
    K = msa.shape[0]
    samples = []

    for s in range(n_sequences):
        if (s + 1) % 25 == 0:
            print(f"  Generating sequence {s+1}/{n_sequences}...", flush=True)

        # Initialize from random stored sequence
        x = msa[rng.integers(0, K)].copy()

        # Gibbs sweeps
        for sweep in range(n_sweeps):
            order = rng.permutation(L)
            for i in order:
                # P(x_i = a | x_{-i}) = softmax(h_i(a) + sum_{j!=i} J_ij(a, x_j))
                logits = h[i].copy()
                for j in range(L):
                    if j != i:
                        logits += J[i, j, :, x[j]]
                logits -= logits.max()
                probs = np.exp(logits)
                probs /= probs.sum()
                x[i] = rng.choice(q, p=probs)

        samples.append(x.copy())

    return np.array(samples)


def write_fasta(sequences, filepath, prefix="potts"):
    """Write sequences in FASTA format."""
    with open(filepath, "w") as f:
        for i, seq_arr in enumerate(sequences):
            seq = "".join(ALPHABET[a] for a in seq_arr)
            clean_seq = seq.rstrip("-")
            f.write(f">{prefix}_seq_{i+1}\n{clean_seq}\n")


def main():
    for family in FAMILIES:
        fam_id = family["id"]
        fam_name = family["name"]

        print(f"\n{'='*70}")
        print(f"Processing {fam_name} ({fam_id})")
        print(f"{'='*70}")

        fasta_path = DATA_DIR / fam_id / "aligned_for_evodiff.fasta"
        if not fasta_path.exists():
            print(f"  Alignment not found at {fasta_path}, skipping.")
            continue

        sequences = read_aligned_fasta(str(fasta_path))
        msa = encode_alignment(sequences)
        K, L = msa.shape
        n_couplings = L * (L - 1) // 2 * Q * Q
        n_fields = L * Q
        print(f"  Alignment: K={K} sequences, L={L} positions, q={Q} states")
        print(f"  Parameters: {n_fields} fields + {n_couplings} couplings = {n_fields + n_couplings}")
        print(f"  Data/parameter ratio: {K / (n_fields + n_couplings):.6f}")

        # Compute sequence weights
        print(f"\n  Computing sequence weights...")
        weights = compute_weights(msa)
        print(f"  K_eff = {1/np.sum(weights**2):.1f}")

        # Fit plmDCA
        print(f"\n  Fitting pseudo-likelihood DCA...")
        t0 = time.time()
        h, J = fit_plmdca(msa, weights, lambda_J=LAMBDA_J, lambda_h=LAMBDA_H)
        fit_time = time.time() - t0
        print(f"  plmDCA fit completed in {fit_time:.1f}s")
        print(f"  Max |h|: {np.abs(h).max():.3f}, Max |J|: {np.abs(J).max():.3f}")

        # Generate sequences
        print(f"\n  Generating {N_GEN} sequences via Gibbs sampling ({GIBBS_SWEEPS} sweeps)...")
        t0 = time.time()
        samples = gibbs_sample(h, J, msa, n_sequences=N_GEN, n_sweeps=GIBBS_SWEEPS)
        gen_time = time.time() - t0
        print(f"  Generation completed in {gen_time:.1f}s")

        # Save
        out_path = DATA_DIR / fam_id / "potts_sequences.fasta"
        write_fasta(samples, out_path, prefix=f"{fam_name}_potts")
        print(f"  Saved to {out_path}")

        # Quick diagnostics
        print(f"\n  Diagnostics:")

        # Nearest sequence identity
        identities = []
        for gen_seq in samples:
            max_id = max(np.mean(gen_seq == msa[k]) for k in range(K))
            identities.append(max_id)
        print(f"  Nearest seq identity: {np.mean(identities):.3f} +/- {np.std(identities):.3f}")

        # AA composition KL
        gen_counts = np.bincount(samples.ravel(), minlength=Q).astype(float)
        stored_counts = np.bincount(msa.ravel(), minlength=Q).astype(float)
        gen_freq = gen_counts / gen_counts.sum()
        stored_freq = stored_counts / stored_counts.sum()
        mask = (gen_freq > 0) & (stored_freq > 0)
        kl = np.sum(gen_freq[mask] * np.log(gen_freq[mask] / stored_freq[mask]))
        print(f"  AA composition KL: {kl:.4f}")

        # Gap fraction
        gap_idx = AA_TO_IDX["-"]
        print(f"  Gap fraction - generated: {np.mean(samples == gap_idx):.3f}, "
              f"stored: {np.mean(msa == gap_idx):.3f}")

        print(f"\n  Total time: {fit_time + gen_time:.1f}s")


if __name__ == "__main__":
    main()
