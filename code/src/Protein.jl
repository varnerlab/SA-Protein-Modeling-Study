# ──────────────────────────────────────────────────────────────────────────────
# Protein.jl — Shared data pipeline for protein sequence experiments
#
# Provides: MSA download/parsing, one-hot encoding, PCA helpers, sequence
# metrics, and phase-transition analysis. All experiment folders import this.
# ──────────────────────────────────────────────────────────────────────────────

# ══════════════════════════════════════════════════════════════════════════════
# Amino acid alphabet and encoding constants
# ══════════════════════════════════════════════════════════════════════════════

const AA_ALPHABET = collect("ACDEFGHIKLMNPQRSTVWY")  # 20 standard amino acids
const AA_TO_IDX = Dict(aa => i for (i, aa) in enumerate(AA_ALPHABET))
const N_AA = length(AA_ALPHABET)  # 20

# ══════════════════════════════════════════════════════════════════════════════
# Data I/O: download, parse, clean
# ══════════════════════════════════════════════════════════════════════════════

"""
    download_pfam_seed(pfam_id; cache_dir) -> String

Download the Pfam seed alignment in Stockholm format from InterPro.
Returns the path to the cached file. Skips download if already cached.
"""
function download_pfam_seed(pfam_id::String; cache_dir::String=".")
    mkpath(cache_dir)
    cache_file = joinpath(cache_dir, "$(pfam_id)_seed.sto")

    if isfile(cache_file)
        @info "  Using cached alignment: $cache_file"
        return cache_file
    end

    url = "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/$(pfam_id)/?annotation=alignment:seed"
    gz_file = cache_file * ".gz"
    @info "  Downloading seed alignment from InterPro …"
    @info "  URL: $url"

    try
        Downloads.download(url, gz_file)
        run(`gunzip -f $gz_file`)
        @info "  Saved to $cache_file"
    catch e
        isfile(gz_file) && rm(gz_file)
        @warn "  Download failed: $e"
        @warn "  Please download the $(pfam_id) seed alignment manually and place it at:"
        @warn "    $cache_file"
        error("Cannot proceed without alignment data.")
    end

    return cache_file
end

"""
    parse_stockholm(filepath) -> Vector{Tuple{String, String}}

Parse a Stockholm-format multiple sequence alignment.
Returns vector of (name, aligned_sequence) tuples.
Handles multi-block Stockholm files by concatenating sequence fragments.
"""
function parse_stockholm(filepath::String)
    sequences = Dict{String, String}()
    seq_order = String[]

    for line in eachline(filepath)
        startswith(line, "#") && continue
        startswith(line, "//") && continue
        stripped = strip(line)
        isempty(stripped) && continue

        parts = split(stripped)
        length(parts) >= 2 || continue

        name = parts[1]
        seq = uppercase(parts[2])

        if haskey(sequences, name)
            sequences[name] *= seq  # append for multi-block Stockholm
        else
            sequences[name] = seq
            push!(seq_order, name)
        end
    end

    return [(name, sequences[name]) for name in seq_order]
end

"""
    parse_fasta(filepath) -> Vector{Tuple{String, String}}

Parse a FASTA-format file. Returns vector of (name, sequence) tuples.
"""
function parse_fasta(filepath::String)
    sequences = Tuple{String,String}[]
    current_name = ""
    current_seq = IOBuffer()

    for line in eachline(filepath)
        if startswith(line, ">")
            if !isempty(current_name)
                push!(sequences, (current_name, String(take!(current_seq))))
            end
            current_name = strip(line[2:end])
            current_seq = IOBuffer()
        else
            write(current_seq, uppercase(strip(line)))
        end
    end
    if !isempty(current_name)
        push!(sequences, (current_name, String(take!(current_seq))))
    end

    return sequences
end

"""
    clean_alignment(raw_seqs; max_gap_frac_col, max_gap_frac_seq) -> (Matrix{Char}, Vector{String})

Process a raw alignment: remove high-gap columns and high-gap sequences.
Returns a character matrix (K × L) and the retained sequence names.
"""
function clean_alignment(raw_seqs::Vector{Tuple{String,String}};
                         max_gap_frac_col::Float64=0.5,
                         max_gap_frac_seq::Float64=0.3)

    names = [s[1] for s in raw_seqs]
    seqs  = [s[2] for s in raw_seqs]

    L_raw = length(seqs[1])
    K_raw = length(seqs)
    char_mat = fill('.', K_raw, L_raw)
    for (i, seq) in enumerate(seqs)
        for (j, c) in enumerate(seq)
            j <= L_raw && (char_mat[i, j] = c)
        end
    end

    is_gap(c) = c in ('.', '-', '~')

    # remove columns with >max_gap_frac_col gaps
    col_gap_frac = [count(is_gap, char_mat[:, j]) / K_raw for j in 1:L_raw]
    keep_cols = findall(f -> f <= max_gap_frac_col, col_gap_frac)
    char_mat = char_mat[:, keep_cols]
    L = length(keep_cols)

    # remove sequences with >max_gap_frac_seq gaps (in remaining columns)
    seq_gap_frac = [count(is_gap, char_mat[i, :]) / L for i in 1:K_raw]
    keep_seqs = findall(f -> f <= max_gap_frac_seq, seq_gap_frac)
    char_mat = char_mat[keep_seqs, :]
    names = names[keep_seqs]

    @info "  Alignment: $(K_raw) seqs × $(L_raw) cols → $(length(keep_seqs)) seqs × $L cols (after gap filtering)"

    return char_mat, names
end

# ══════════════════════════════════════════════════════════════════════════════
# One-hot encoding / decoding
# ══════════════════════════════════════════════════════════════════════════════

"""
    onehot_encode(char_mat) -> Matrix{Float64}

One-hot encode a character matrix (K × L) into a matrix (20L × K).
Gap positions and non-standard amino acids map to all-zeros (20-dim zero vector).
"""
function onehot_encode(char_mat::Matrix{Char})
    K, L = size(char_mat)
    d_full = N_AA * L  # 20L
    X = zeros(Float64, d_full, K)

    for k in 1:K
        for pos in 1:L
            aa = char_mat[k, pos]
            idx = get(AA_TO_IDX, aa, 0)
            if idx > 0
                X[(pos-1)*N_AA + idx, k] = 1.0
            end
        end
    end

    return X
end

"""
    decode_onehot(x, L) -> String

Decode a continuous vector (20L-dim) back to an amino acid sequence.
At each position, take the argmax over the 20 amino acid channels.
Positions where all channels are near-zero are decoded as gaps ('-').
"""
function decode_onehot(x::Vector{Float64}, L::Int)
    seq = Char[]
    for pos in 1:L
        start_idx = (pos - 1) * N_AA + 1
        end_idx = pos * N_AA
        block = x[start_idx:end_idx]
        best_idx = argmax(block)
        if maximum(block) < 1e-10
            push!(seq, '-')
        else
            push!(seq, AA_ALPHABET[best_idx])
        end
    end
    return String(seq)
end

# ══════════════════════════════════════════════════════════════════════════════
# PCA pipeline: encode → reduce → normalize → memory matrix
# ══════════════════════════════════════════════════════════════════════════════

"""
    build_memory_matrix(char_mat; pratio=0.95) -> (X̂, pca_model, L, d_full)

Full pipeline: one-hot encode → PCA → unit-norm columns.

Returns:
- `X̂`: Memory matrix (d_pca × K), unit-norm columns in PCA space
- `pca_model`: Fitted PCA model (for inverse transform / decoding)
- `L`: Alignment length (number of positions)
- `d_full`: Full one-hot dimensionality (20L)
"""
function build_memory_matrix(char_mat::Matrix{Char}; pratio::Float64=0.95)
    K, L = size(char_mat)

    # one-hot encode
    X_onehot = onehot_encode(char_mat)  # (20L × K)
    d_full = size(X_onehot, 1)
    @info "  One-hot: $d_full × $K"

    # PCA
    pca_model = MultivariateStats.fit(PCA, X_onehot; pratio=pratio)
    d_pca = MultivariateStats.outdim(pca_model)
    Z = MultivariateStats.transform(pca_model, X_onehot)  # d_pca × K
    var_retained = round(100 * sum(MultivariateStats.principalvars(pca_model)) /
                         MultivariateStats.tvar(pca_model), digits=1)
    @info "  PCA: $d_full → $d_pca dimensions ($var_retained% variance)"

    # normalize to unit norm
    X̂ = copy(Z)
    for k in 1:K
        nk = norm(X̂[:, k])
        X̂[:, k] ./= (nk + 1e-12)
    end
    @info "  Memory matrix X̂: $d_pca × $K (unit-norm columns)"

    return X̂, pca_model, L, d_full
end

"""
    decode_sample(ξ_pca, pca_model, L) -> String

Decode a PCA-space vector back to an amino acid sequence.
Maps through inverse PCA, then argmax decoding at each position.
"""
function decode_sample(ξ_pca::Vector{Float64}, pca_model, L::Int)
    x_onehot = vec(MultivariateStats.reconstruct(pca_model, ξ_pca))
    return decode_onehot(x_onehot, L)
end

# ══════════════════════════════════════════════════════════════════════════════
# Protein-specific evaluation metrics
# ══════════════════════════════════════════════════════════════════════════════

"""
    sequence_identity(seq1, seq2) -> Float64

Fraction of positions where two sequences have the same amino acid.
Gaps are excluded from the comparison.
"""
function sequence_identity(seq1::String, seq2::String)
    L = min(length(seq1), length(seq2))
    matches = 0
    compared = 0
    for i in 1:L
        (seq1[i] == '-' || seq2[i] == '-') && continue
        compared += 1
        seq1[i] == seq2[i] && (matches += 1)
    end
    return compared > 0 ? matches / compared : 0.0
end

"""
    nearest_sequence_identity(gen_seq, stored_seqs) -> Float64

Maximum sequence identity between a generated sequence and any stored sequence.
"""
function nearest_sequence_identity(gen_seq::String, stored_seqs::Vector{String})
    return maximum(sequence_identity(gen_seq, s) for s in stored_seqs)
end

"""
    valid_residue_fraction(seq) -> Float64

Fraction of non-gap positions that are standard amino acids.
"""
function valid_residue_fraction(seq::String)
    non_gap = count(c -> c != '-', seq)
    non_gap == 0 && return 0.0
    valid = count(c -> c in AA_ALPHABET, seq)
    return valid / non_gap
end

"""
    aa_composition_kl(gen_seqs, stored_seqs) -> Float64

KL divergence of amino acid composition: D_KL(stored || generated).
Measures how well generated sequences preserve the family's amino acid frequencies.
"""
function aa_composition_kl(gen_seqs::Vector{String}, stored_seqs::Vector{String})
    function aa_freqs(seqs)
        counts = zeros(N_AA)
        for seq in seqs, c in seq
            idx = get(AA_TO_IDX, c, 0)
            idx > 0 && (counts[idx] += 1)
        end
        total = sum(counts)
        total > 0 ? counts ./ total : ones(N_AA) ./ N_AA
    end

    p = aa_freqs(stored_seqs)  # reference
    q = aa_freqs(gen_seqs)     # generated

    eps = 1e-10
    p .+= eps; p ./= sum(p)
    q .+= eps; q ./= sum(q)

    return sum(p[i] * log(p[i] / q[i]) for i in 1:N_AA)
end

"""
    aa_freq_matrix(seqs, L) -> Matrix{Float64}

Compute per-position amino acid frequency matrix (N_AA × L).
Each column sums to 1 (or 0 if no residues at that position).
"""
function aa_freq_matrix(seqs::Vector{String}, L::Int)
    freq = zeros(N_AA, L)
    for seq in seqs
        for pos in 1:min(L, length(seq))
            idx = get(AA_TO_IDX, seq[pos], 0)
            idx > 0 && (freq[idx, pos] += 1)
        end
    end
    for j in 1:L
        s = sum(freq[:, j])
        s > 0 && (freq[:, j] ./= s)
    end
    return freq
end

# ══════════════════════════════════════════════════════════════════════════════
# Phase transition analysis
# ══════════════════════════════════════════════════════════════════════════════

"""
    find_entropy_inflection(X̂; α, n_betas, β_range) -> NamedTuple

Compute the entropy inflection point β* for the memory matrix X̂.
Returns β*, SNR*, theoretical predictions, and the full β/H curves.
"""
function find_entropy_inflection(X̂::Matrix{Float64};
                                  α::Float64=0.01,
                                  n_betas::Int=50,
                                  β_range::Tuple{Float64,Float64}=(0.1, 500.0))

    d, K = size(X̂)
    βs = 10 .^ range(log10(β_range[1]), log10(β_range[2]), length=n_betas)

    # compute mean entropy at each β using a few random probes
    n_probes = min(K, 20)
    Hs = zeros(n_betas)
    for (bi, β) in enumerate(βs)
        H_sum = 0.0
        for k in 1:n_probes
            H_sum += attention_entropy(X̂[:, k], X̂, β)
        end
        Hs[bi] = H_sum / n_probes
    end

    # numerical second derivative to find inflection (in log-β space)
    log_βs = log.(βs)
    dH = diff(Hs) ./ diff(log_βs)
    d2H = diff(dH) ./ diff(log_βs[1:end-1])

    # inflection: where d2H is most negative
    inflection_idx = 1
    min_d2H = Inf
    for i in 1:length(d2H)
        if d2H[i] < min_d2H
            min_d2H = d2H[i]
            inflection_idx = i + 1
        end
    end

    β_star = βs[inflection_idx]
    snr_star = sqrt(α * β_star / (2 * d))

    # theoretical prediction for random unit-norm patterns
    β_star_theory = sqrt(d)
    snr_star_theory = sqrt(α / (2 * sqrt(d)))

    @info "  Phase transition analysis (d=$d, K=$K):"
    @info "    Empirical inflection:  β* = $(round(β_star, digits=2)),  SNR* = $(round(snr_star, digits=4))"
    @info "    Theoretical (√d):      β* = $(round(β_star_theory, digits=2)),  SNR* = $(round(snr_star_theory, digits=4))"

    return (β_star=β_star, snr_star=snr_star,
            β_star_theory=β_star_theory, snr_star_theory=snr_star_theory,
            βs=βs, Hs=Hs)
end

# ══════════════════════════════════════════════════════════════════════════════
# MSA statistics (for β* prediction — Section 2b of paper)
# ══════════════════════════════════════════════════════════════════════════════

"""
    msa_column_entropy(char_mat) -> Vector{Float64}

Compute Shannon entropy at each column of the alignment.
Returns a vector of length L (one entropy per position).
"""
function msa_column_entropy(char_mat::Matrix{Char})
    K, L = size(char_mat)
    entropies = zeros(L)
    for j in 1:L
        counts = zeros(N_AA)
        total = 0
        for i in 1:K
            idx = get(AA_TO_IDX, char_mat[i, j], 0)
            if idx > 0
                counts[idx] += 1
                total += 1
            end
        end
        total == 0 && continue
        for c in counts
            p = c / total
            p > 0 && (entropies[j] -= p * log(p))
        end
    end
    return entropies
end

"""
    effective_num_sequences(char_mat; threshold=0.8) -> Float64

Compute K_eff: the effective number of sequences after reweighting
by pairwise sequence identity. Sequences with identity > threshold
are downweighted. Standard practice in DCA/coevolution analysis.
"""
function effective_num_sequences(char_mat::Matrix{Char}; threshold::Float64=0.8)
    K, L = size(char_mat)
    weights = ones(K)
    for i in 1:K
        n_similar = 0
        for j in 1:K
            i == j && continue
            matches = 0
            compared = 0
            for pos in 1:L
                ci, cj = char_mat[i, pos], char_mat[j, pos]
                (ci in ('.', '-', '~') || cj in ('.', '-', '~')) && continue
                compared += 1
                ci == cj && (matches += 1)
            end
            if compared > 0 && (matches / compared) >= threshold
                n_similar += 1
            end
        end
        weights[i] = 1.0 / (1.0 + n_similar)
    end
    return sum(weights)
end
