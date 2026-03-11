#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Protein Sequence Experiment (Pfam RRM family, PF00076)
# Mirrors the protocol in run_multidigit_experiment.jl (MNIST digit 3)
# Runs all baselines + SA/MALA, computes metrics, saves figures.
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment …"
include(joinpath(@__DIR__, "Include-Sequence.jl"))
@info "Environment loaded."

# ══════════════════════════════════════════════════════════════════════════════
# Part 1: Constants and amino acid alphabet
# ══════════════════════════════════════════════════════════════════════════════
const AA_ALPHABET = collect("ACDEFGHIKLMNPQRSTVWY")  # 20 standard amino acids
const AA_TO_IDX = Dict(aa => i for (i, aa) in enumerate(AA_ALPHABET))
const N_AA = length(AA_ALPHABET)  # 20

# ══════════════════════════════════════════════════════════════════════════════
# Part 2: Data loading — download and parse Pfam seed alignment
# ══════════════════════════════════════════════════════════════════════════════

"""
    download_pfam_seed(pfam_id; cache_dir) -> String

Download the Pfam seed alignment in Stockholm format from InterPro.
Returns the path to the cached file.
"""
function download_pfam_seed(pfam_id::String; cache_dir::String=_PATH_TO_DATA)
    mkpath(cache_dir)
    cache_file = joinpath(cache_dir, "$(pfam_id)_seed.sto")

    if isfile(cache_file)
        @info "  Using cached alignment: $cache_file"
        return cache_file
    end

    # InterPro API returns gzip-compressed Stockholm format
    url = "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/$(pfam_id)/?annotation=alignment:seed"
    gz_file = cache_file * ".gz"
    @info "  Downloading seed alignment from InterPro …"
    @info "  URL: $url"

    try
        Downloads.download(url, gz_file)
        run(`gunzip -f $gz_file`)
        @info "  Saved to $cache_file"
    catch e
        # clean up partial downloads
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
"""
function parse_stockholm(filepath::String)
    sequences = Dict{String, String}()
    seq_order = String[]

    for line in eachline(filepath)
        # skip comments, metadata, and end markers
        startswith(line, "#") && continue
        startswith(line, "//") && continue
        stripped = strip(line)
        isempty(stripped) && continue

        # sequence lines: "name/start-end   aligned_sequence"
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
Returns a character matrix (K × L) and sequence names.
"""
function clean_alignment(raw_seqs::Vector{Tuple{String,String}};
                         max_gap_frac_col::Float64=0.5,
                         max_gap_frac_seq::Float64=0.3)

    names = [s[1] for s in raw_seqs]
    seqs  = [s[2] for s in raw_seqs]

    # convert to character matrix
    L_raw = length(seqs[1])
    K_raw = length(seqs)
    char_mat = fill('.', K_raw, L_raw)
    for (i, seq) in enumerate(seqs)
        for (j, c) in enumerate(seq)
            j <= L_raw && (char_mat[i, j] = c)
        end
    end

    # identify gap characters
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
# Part 3: One-hot encoding and PCA
# ══════════════════════════════════════════════════════════════════════════════

"""
    onehot_encode(char_mat) -> Matrix{Float64}

One-hot encode a character matrix (K × L) into a matrix (20L × K).
Gap positions map to all-zeros (20-dim zero vector).
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
            # gaps and non-standard AAs → zero vector (already 0)
        end
    end

    return X
end

"""
    decode_onehot(x, L) -> String

Decode a continuous vector (20L-dim) back to an amino acid sequence.
At each position, take the argmax over the 20 amino acid channels.
"""
function decode_onehot(x::Vector{Float64}, L::Int)
    seq = Char[]
    for pos in 1:L
        start_idx = (pos - 1) * N_AA + 1
        end_idx = pos * N_AA
        block = x[start_idx:end_idx]
        best_idx = argmax(block)
        if maximum(block) < 1e-10
            push!(seq, '-')  # gap if all near-zero
        else
            push!(seq, AA_ALPHABET[best_idx])
        end
    end
    return String(seq)
end

# ══════════════════════════════════════════════════════════════════════════════
# Part 4: Protein-specific evaluation metrics
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

KL divergence of amino acid composition: generated vs stored sequences.
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

    # add small epsilon to avoid log(0)
    eps = 1e-10
    p .+= eps; p ./= sum(p)
    q .+= eps; q ./= sum(q)

    return sum(p[i] * log(p[i] / q[i]) for i in 1:N_AA)
end

# ══════════════════════════════════════════════════════════════════════════════
# Part 5: Phase transition analysis (entropy inflection)
# ══════════════════════════════════════════════════════════════════════════════

"""
    find_entropy_inflection(X̂; α, n_betas, β_range) -> (β_star, snr_star, βs, Hs)

Compute the entropy inflection point β* for the memory matrix X̂.
Returns β*, SNR*, and the full β/H curves for plotting.
"""
function find_entropy_inflection(X̂::Matrix{Float64};
                                  α::Float64=0.01,
                                  n_betas::Int=50,
                                  β_range::Tuple{Float64,Float64}=(0.1, 500.0))

    d, K = size(X̂)
    βs = 10 .^ range(log10(β_range[1]), log10(β_range[2]), length=n_betas)

    # compute mean entropy at each β using a few random states
    n_probes = min(K, 20)
    Hs = zeros(n_betas)
    for (bi, β) in enumerate(βs)
        H_sum = 0.0
        for k in 1:n_probes
            H_sum += attention_entropy(X̂[:, k], X̂, β)
        end
        Hs[bi] = H_sum / n_probes
    end

    # numerical second derivative to find inflection
    # use finite differences on log-spaced β
    log_βs = log.(βs)
    dH = diff(Hs) ./ diff(log_βs)         # first derivative in log-β space
    d2H = diff(dH) ./ diff(log_βs[1:end-1])  # second derivative

    # inflection: where d2H changes sign (most negative → crossing zero)
    inflection_idx = 1
    min_d2H = Inf
    for i in 1:length(d2H)
        if d2H[i] < min_d2H
            min_d2H = d2H[i]
            inflection_idx = i + 1  # offset for double diff
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
# Part 6: GMM-PCA helpers (same as MNIST experiment)
# ══════════════════════════════════════════════════════════════════════════════
# ── VAE struct (defined once at top level) ────────────────────────────────
struct _ProteinVAE
    enc_shared; enc_μ; enc_logσ²; dec
end
Flux.@layer _ProteinVAE

const R_pca_gmm = 50   # PCA components for GMM baseline
const C_gmm     = 10   # GMM mixture components
const EM_iters  = 500
const EM_tol    = 1e-8

function _logsumexp_cols(A::Matrix{Float64})
    m = maximum(A, dims=1)
    return m .+ log.(sum(exp.(A .- m), dims=1))
end

function _fit_gmm(X::Matrix{Float64}, C::Int; n_iters=EM_iters, tol=EM_tol, seed=42)
    r, N = size(X)
    rng  = MersenneTwister(seed)
    # K-means++ init
    μ = zeros(r, C)
    μ[:, 1] = X[:, rand(rng, 1:N)]
    for k in 2:C
        dists = [minimum(sum((X[:, n] .- μ[:, j]).^2) for j in 1:(k-1)) for n in 1:N]
        μ[:, k] = X[:, StatsBase.sample(rng, 1:N, Weights(dists ./ sum(dists)))]
    end
    log_σ² = log.(repeat(max.(vec(var(X, dims=2)), 1e-6), 1, C))
    π_k    = ones(C) / C
    log_r  = zeros(C, N)
    log_lik_prev = -Inf
    for _ in 1:n_iters
        σ² = exp.(log_σ²)
        for k in 1:C
            diff = X .- μ[:, k]
            log_r[k, :] = log(π_k[k]) .-
                           0.5 .* vec(sum(diff.^2 ./ σ²[:, k], dims=1)) .-
                           0.5 .* sum(log.(2π .* σ²[:, k]))
        end
        lse     = _logsumexp_cols(log_r)
        log_lik = sum(lse)
        log_r .-= lse
        r_mat   = exp.(log_r)
        abs(log_lik - log_lik_prev) / (abs(log_lik_prev) + 1e-30) < tol && break
        log_lik_prev = log_lik
        N_k = vec(sum(r_mat, dims=2))
        π_k = max.(N_k / N, 1e-8);  π_k ./= sum(π_k)
        for k in 1:C
            μ[:, k]      = X * r_mat[k, :] / N_k[k]
            diff         = X .- μ[:, k]
            log_σ²[:, k] = log.(max.(vec(sum(diff.^2 .* r_mat[k, :]', dims=2)) / N_k[k], 1e-6))
        end
    end
    return μ, log_σ², π_k
end

function _sample_gmm(μ, log_σ², π_k, n_samples; seed=9999)
    rng = MersenneTwister(seed)
    r   = size(μ, 1)
    wts = Weights(π_k)
    [let k = StatsBase.sample(rng, 1:length(π_k), wts)
         μ[:, k] .+ sqrt.(exp.(log_σ²[:, k])) .* randn(rng, r)
     end for _ in 1:n_samples]
end

# ══════════════════════════════════════════════════════════════════════════════
# Part 7: Main experiment
# ══════════════════════════════════════════════════════════════════════════════

# -- Experiment parameters (mirrors MNIST protocol) --
const PFAM_ID = "PF00076"  # RRM (RNA Recognition Motif)
const K_MAX = 100           # max number of sequences to use
const α_step = 0.01
const S = 150               # total samples to generate
const n_chains = 30
const T_per_chain = 5000
const T_burnin = 2000
const thin_interval = 100
const samples_per_chain = 5
const σ_init = 0.01
const R_pca_sa = 0.95       # fraction of variance to retain in PCA

function run_protein_experiment(; pfam_id=PFAM_ID, figpath=_PATH_TO_FIG)

    mkpath(figpath)
    mkpath(_PATH_TO_DATA)

    @info "═══  Protein Sequence Experiment: $pfam_id  ═══"

    # ── Step 1: Download and parse alignment ──────────────────────────────────
    @info "Step 1: Loading alignment data …"
    sto_file = download_pfam_seed(pfam_id)
    raw_seqs = parse_stockholm(sto_file)

    if isempty(raw_seqs)
        # try FASTA format as fallback
        @info "  No Stockholm sequences found, trying FASTA …"
        raw_seqs = parse_fasta(sto_file)
    end

    @info "  Parsed $(length(raw_seqs)) sequences"

    # ── Step 2: Clean alignment ───────────────────────────────────────────────
    char_mat, seq_names = clean_alignment(raw_seqs)
    K_total, L = size(char_mat)

    # subsample to K_MAX if needed
    if K_total > K_MAX
        Random.seed!(42)
        keep = StatsBase.sample(1:K_total, K_MAX, replace=false) |> sort
        char_mat = char_mat[keep, :]
        seq_names = seq_names[keep]
        @info "  Subsampled to $K_MAX sequences"
    end
    K = size(char_mat, 1)

    # store original sequences for evaluation
    stored_seqs = [String(char_mat[k, :]) for k in 1:K]

    # ── Step 3: One-hot encode ────────────────────────────────────────────────
    @info "Step 3: One-hot encoding ($K seqs × $L positions → $(N_AA*L) dimensions) …"
    X_onehot = onehot_encode(char_mat)  # (20L × K)
    d_full = size(X_onehot, 1)
    @info "  One-hot matrix: $d_full × $K"

    # ── Step 4: PCA reduction ─────────────────────────────────────────────────
    @info "Step 4: PCA reduction (retaining $(Int(R_pca_sa*100))% variance) …"
    pca_model = MultivariateStats.fit(PCA, X_onehot; pratio=R_pca_sa)
    d_pca = outdim(pca_model)  # number of PCA components retained
    Z = MultivariateStats.transform(pca_model, X_onehot)  # d_pca × K
    @info "  PCA: $d_full → $d_pca dimensions ($(round(100*sum(principalvars(pca_model))/tvar(pca_model), digits=1))% variance)"

    # ── Step 5: Normalize to unit norm (memory matrix) ────────────────────────
    ϵ = 1e-12
    X̂ = copy(Z)
    for k in 1:K
        nk = norm(X̂[:, k])
        X̂[:, k] ./= (nk + ϵ)
    end
    d = size(X̂, 1)
    @info "  Memory matrix X̂: $d × $K (unit-norm columns in PCA space)"

    # ── Step 6: Phase transition analysis ─────────────────────────────────────
    @info "Step 6: Phase transition analysis …"
    pt = find_entropy_inflection(X̂; α=α_step)
    β_retrieval = round(Int, 20 * pt.β_star)  # ~20× the transition for structured retrieval
    β_generation = round(Int, 2 * pt.β_star)  # ~2× the transition for generation
    @info "  Selected β_retrieval = $β_retrieval,  β_generation = $β_generation"

    # save entropy curve
    p_entropy = plot(pt.βs, pt.Hs ./ log(K),
        xlabel="β (inverse temperature)", ylabel="H(β) / log K",
        title="Attention entropy — $pfam_id (d=$d, K=$K)",
        xscale=:log10, label="H(β)/log K", lw=2, color=:coral,
        legend=:topright, size=(600, 400))
    vline!([pt.β_star], label="β* = $(round(pt.β_star, digits=1)) (empirical)", ls=:dash, color=:blue)
    vline!([pt.β_star_theory], label="β* = $(round(pt.β_star_theory, digits=1)) (√d theory)", ls=:dot, color=:green)
    savefig(p_entropy, joinpath(figpath, "Fig_protein_entropy_curve.pdf"))
    @info "  Saved entropy curve figure"

    # ── Step 7: SA multi-chain (retrieval regime) ─────────────────────────────
    β_inv_temp = Float64(β_retrieval)
    @info "Step 7: Running SA multi-chain (β=$β_retrieval) …"
    sa_samples = Vector{Vector{Float64}}()
    Random.seed!(42)
    pattern_indices = StatsBase.sample(1:K, n_chains, replace=(n_chains > K))
    for (c, k) in enumerate(pattern_indices)
        Random.seed!(12345 + c)
        sₒ = X̂[:, k] .+ σ_init .* randn(d)
        (_, Ξ) = sample(X̂, sₒ, T_per_chain; β=β_inv_temp, α=α_step, seed=12345+c)
        chain_pool = Vector{Vector{Float64}}()
        for tᵢ in (T_burnin+1):thin_interval:T_per_chain
            push!(chain_pool, Ξ[tᵢ, :])
        end
        n_avail = length(chain_pool)
        idxs = round.(Int, range(1, n_avail, length=min(samples_per_chain, n_avail)))
        for idx in idxs
            push!(sa_samples, chain_pool[idx])
        end
    end
    @info "  SA (retrieval): $(length(sa_samples)) samples"

    # ── Step 7b: SA multi-chain (generation regime) ───────────────────────────
    β_gen_temp = Float64(β_generation)
    @info "Step 7b: Running SA multi-chain GENERATION (β=$β_generation) …"
    sa_gen_samples = Vector{Vector{Float64}}()
    Random.seed!(42)
    pattern_indices_gen = StatsBase.sample(1:K, n_chains, replace=(n_chains > K))
    for (c, k) in enumerate(pattern_indices_gen)
        Random.seed!(12345 + c)
        sₒ = X̂[:, k] .+ σ_init .* randn(d)
        (_, Ξ) = sample(X̂, sₒ, T_per_chain; β=β_gen_temp, α=α_step, seed=12345+c)
        chain_pool = Vector{Vector{Float64}}()
        for tᵢ in (T_burnin+1):thin_interval:T_per_chain
            push!(chain_pool, Ξ[tᵢ, :])
        end
        n_avail = length(chain_pool)
        idxs = round.(Int, range(1, n_avail, length=min(samples_per_chain, n_avail)))
        for idx in idxs
            push!(sa_gen_samples, chain_pool[idx])
        end
    end
    @info "  SA (generation): $(length(sa_gen_samples)) samples"

    # ── Step 8: MALA multi-chain ──────────────────────────────────────────────
    @info "Step 8: Running MALA multi-chain (β=$β_retrieval) …"
    mala_samples = Vector{Vector{Float64}}()
    mala_accept_rates = Float64[]
    Random.seed!(42)
    pattern_indices = StatsBase.sample(1:K, n_chains, replace=(n_chains > K))
    for (c, k) in enumerate(pattern_indices)
        Random.seed!(12345 + c)
        sₒ = X̂[:, k] .+ σ_init .* randn(d)
        (_, Ξ, ar) = mala_sample(X̂, sₒ, T_per_chain; β=β_inv_temp, α=α_step, seed=12345+c)
        push!(mala_accept_rates, ar)
        chain_pool = Vector{Vector{Float64}}()
        for tᵢ in (T_burnin+1):thin_interval:T_per_chain
            push!(chain_pool, Ξ[tᵢ, :])
        end
        n_avail = length(chain_pool)
        idxs = round.(Int, range(1, n_avail, length=min(samples_per_chain, n_avail)))
        for idx in idxs
            push!(mala_samples, chain_pool[idx])
        end
    end
    mala_mean_ar = round(mean(mala_accept_rates), digits=4)
    @info "  MALA: $(length(mala_samples)) samples, accept rate = $mala_mean_ar"

    # ── Step 9: Baselines ─────────────────────────────────────────────────────
    @info "Step 9: Running baselines …"

    # Bootstrap (replay)
    bs_samples = Vector{Vector{Float64}}()
    Random.seed!(12345)
    for _ in 1:S
        k = rand(1:K)
        push!(bs_samples, copy(X̂[:, k]))
    end

    # Gaussian perturbation
    gp_samples = Vector{Vector{Float64}}()
    σ_noise = sqrt(2 * α_step / β_inv_temp)
    Random.seed!(12345)
    for _ in 1:S
        k = rand(1:K)
        ξ = X̂[:, k] .+ σ_noise .* randn(d)
        push!(gp_samples, ξ)
    end

    # Random convex combination
    rc_samples = Vector{Vector{Float64}}()
    dirichlet_dist = Dirichlet(K, 1.0)
    Random.seed!(12345)
    for _ in 1:S
        w = rand(dirichlet_dist)
        ξ = X̂ * w
        push!(rc_samples, ξ)
    end

    # GMM-PCA
    @info "  Fitting GMM-PCA (r=$(min(R_pca_gmm, d)), C=$C_gmm) …"
    R_gmm = min(R_pca_gmm, d - 1)  # can't have more PCA components than dimensions
    pca_gmm_model = MultivariateStats.fit(PCA, X̂; maxoutdim=R_gmm)
    Z_gmm = MultivariateStats.transform(pca_gmm_model, X̂)
    C_actual = min(C_gmm, K - 1)  # can't have more components than data points
    μ_gmm, log_σ²_gmm, π_gmm = _fit_gmm(Z_gmm, C_actual; seed=42)
    pca_samps = _sample_gmm(μ_gmm, log_σ²_gmm, π_gmm, S; seed=7777)
    gmm_samples = [vec(MultivariateStats.reconstruct(pca_gmm_model, z)) for z in pca_samps]
    @info "  GMM-PCA: $(length(gmm_samples)) samples"

    # ── VAE baseline ─────────────────────────────────────────────────────────
    @info "  Training VAE (latent=8, two-phase) …"
    VAE_LATENT = 8
    VAE_PHASE1 = 2000
    VAE_PHASE2 = 2000
    VAE_LR = 1f-3
    VAE_KL_FINAL = 0.0001f0

    # Architecture scaled to PCA dimension
    _h1 = max(d ÷ 2, 16)
    _h2 = max(d ÷ 4, 8)
    pvae = _ProteinVAE(
        Chain(Dense(d => _h1, relu), Dense(_h1 => _h2, relu)),
        Dense(_h2 => VAE_LATENT),
        Dense(_h2 => VAE_LATENT),
        Chain(Dense(VAE_LATENT => _h2, relu), Dense(_h2 => _h1, relu), Dense(_h1 => d))
    )

    X_train_vae = Float32.(X̂)
    opt_state_vae = Flux.setup(Adam(VAE_LR), pvae)

    # Phase 1: pure autoencoder
    for epoch in 1:VAE_PHASE1
        ε = randn(Float32, VAE_LATENT, K)
        loss, grads = Flux.withgradient(pvae) do m
            h = m.enc_shared(X_train_vae)
            μ = m.enc_μ(h); lσ² = m.enc_logσ²(h)
            z = μ .+ exp.(0.5f0 .* lσ²) .* ε
            o = m.dec(z)
            x̂ = o ./ (sqrt.(sum(o .^ 2; dims=1)) .+ 1f-8)
            mean(sum((X_train_vae .- x̂) .^ 2; dims=1))
        end
        Flux.update!(opt_state_vae, pvae, grads[1])
    end

    # Phase 2: VAE with KL warmup
    for epoch in 1:VAE_PHASE2
        kl_w = VAE_KL_FINAL * Float32(epoch) / Float32(VAE_PHASE2)
        ε = randn(Float32, VAE_LATENT, K)
        loss, grads = Flux.withgradient(pvae) do m
            h = m.enc_shared(X_train_vae)
            μ = m.enc_μ(h); lσ² = m.enc_logσ²(h)
            z = μ .+ exp.(0.5f0 .* lσ²) .* ε
            o = m.dec(z)
            x̂ = o ./ (sqrt.(sum(o .^ 2; dims=1)) .+ 1f-8)
            recon = mean(sum((X_train_vae .- x̂) .^ 2; dims=1))
            kl = -0.5f0 * mean(sum(1f0 .+ lσ² .- μ .^ 2 .- exp.(lσ²); dims=1))
            recon + kl_w * kl
        end
        Flux.update!(opt_state_vae, pvae, grads[1])
    end

    # Generate samples
    Random.seed!(9999)
    Z_vae = randn(Float32, VAE_LATENT, S)
    raw_vae = let o = pvae.dec(Z_vae); o ./ (sqrt.(sum(o .^ 2; dims=1)) .+ 1f-8) end
    vae_samples = [Float64.(raw_vae[:, i]) for i in 1:S]
    @info "  VAE: $(length(vae_samples)) samples"

    # ── Step 10: Decode all samples to amino acid sequences ───────────────────
    @info "Step 10: Decoding samples to amino acid sequences …"

    # helper: PCA-space vector → amino acid sequence
    function decode_sample(ξ_pca::Vector{Float64})
        # un-normalize: we need to map back through PCA
        x_onehot = vec(MultivariateStats.reconstruct(pca_model, ξ_pca))
        return decode_onehot(x_onehot, L)
    end

    # note: for PCA-space samples, we need to handle the normalization.
    # stored patterns were: onehot → PCA → normalize. To decode, we need to
    # scale back, then inverse-PCA. We use the PCA reconstruction directly
    # on the PCA-space samples (SA/MALA/baselines operate in normalized PCA space).

    sa_seqs     = [decode_sample(ξ) for ξ in sa_samples]
    sa_gen_seqs = [decode_sample(ξ) for ξ in sa_gen_samples]
    mala_seqs   = [decode_sample(ξ) for ξ in mala_samples]
    bs_seqs   = [decode_sample(ξ) for ξ in bs_samples]
    gp_seqs   = [decode_sample(ξ) for ξ in gp_samples]
    rc_seqs   = [decode_sample(ξ) for ξ in rc_samples]
    gmm_seqs  = [decode_sample(ξ) for ξ in gmm_samples]
    vae_seqs  = [decode_sample(ξ) for ξ in vae_samples]

    # ── Step 11: Compute metrics ──────────────────────────────────────────────
    @info "Step 11: Computing metrics …"

    # Helper: split S samples into n_chains groups for SE computation
    function chain_metric_se(samps, metric_fn)
        nc = n_chains
        spc = samples_per_chain
        group_vals = [metric_fn(samps[(i-1)*spc+1:i*spc]) for i in 1:nc]
        return std(group_vals) / sqrt(nc)
    end

    methods = [
        "Bootstrap (replay)"        => (bs_samples,   bs_seqs),
        "Gaussian perturbation"     => (gp_samples,   gp_seqs),
        "Random convex combination" => (rc_samples,   rc_seqs),
        "GMM-PCA"                   => (gmm_samples,  gmm_seqs),
        "VAE (latent=8)"            => (vae_samples,  vae_seqs),
        "MALA (β=$β_retrieval)"     => (mala_samples, mala_seqs),
        "SA (β=$β_retrieval, ret)"  => (sa_samples,   sa_seqs),
        "SA (β=$β_generation, gen)" => (sa_gen_samples, sa_gen_seqs),
    ]

    rows = Vector{NamedTuple}()
    for (name, (samps, seqs)) in methods
        # standard metrics (in PCA space)
        novelty_vals  = [sample_novelty(ξ, X̂) for ξ in samps]
        energy_vals   = [hopfield_energy(ξ, X̂, β_inv_temp) for ξ in samps]
        novelty_mean  = mean(novelty_vals)
        diversity_mean = sample_diversity(samps)
        energy_mean   = mean(energy_vals)
        novelty_se    = chain_metric_se(samps, g -> mean(sample_novelty(ξ, X̂) for ξ in g))
        diversity_se  = chain_metric_se(samps, g -> sample_diversity(g))
        energy_se     = chain_metric_se(samps, g -> mean(hopfield_energy(ξ, X̂, β_inv_temp) for ξ in g))

        # protein-specific metrics
        seq_id_vals   = [nearest_sequence_identity(s, stored_seqs) for s in seqs]
        seq_id_mean   = mean(seq_id_vals)
        seq_id_se     = chain_metric_se(seqs, g -> mean(nearest_sequence_identity(s, stored_seqs) for s in g))
        valid_frac    = mean(valid_residue_fraction(s) for s in seqs)
        kl_div        = aa_composition_kl(seqs, stored_seqs)

        push!(rows, (Method=name,
                     Novelty=novelty_mean,   Novelty_SE=novelty_se,
                     Diversity=diversity_mean, Diversity_SE=diversity_se,
                     Energy=energy_mean,      Energy_SE=energy_se,
                     SeqID=seq_id_mean,       SeqID_SE=seq_id_se,
                     ValidAA=valid_frac,      KL_AA=kl_div))
    end
    df = DataFrame(rows)

    # ── Step 12: Print results ────────────────────────────────────────────────
    println("\n══════════════════════════════════════════════════════")
    println("PROTEIN EXPERIMENT RESULTS — $pfam_id")
    println("Memory: K=$K sequences, L=$L positions, d_PCA=$d")
    println("β_retrieval=$β_retrieval (SNR=$(round(sqrt(α_step*β_inv_temp/(2*d)), digits=4)))")
    println("β_generation=$β_generation (SNR=$(round(sqrt(α_step*β_gen_temp/(2*d)), digits=4)))")
    println("β* (entropy inflection) = $(round(pt.β_star, digits=2))")
    println("MALA acceptance rate: $mala_mean_ar")
    println("══════════════════════════════════════════════════════")

    # formatted table
    fmt(x, digits) = let v = round(x; digits); abs(v) == 0.0 ? abs(v) : v end
    fmtpm(v, se, dv, dse) = "\$$(fmt(v,dv)) \\pm $(fmt(se,dse))\$"

    println("\n% --- LaTeX table rows (copy into paper) ---")
    for row in eachrow(df)
        name = row.Method
        if name == "Stochastic attention"
            name = "\\textbf{SA (β=$β_retrieval, retrieval)}"
        end
        nv  = fmtpm(row.Novelty,   row.Novelty_SE,   3, 3)
        div = fmtpm(row.Diversity, row.Diversity_SE, 3, 3)
        en  = fmtpm(row.Energy,    row.Energy_SE,    3, 3)
        sid = fmtpm(row.SeqID,     row.SeqID_SE,     3, 3)
        println("$(name) & $(nv) & $(div) & $(en) & $(sid) \\\\")
    end

    # plain-text summary
    println("\n% --- Summary ---")
    println(df[:, [:Method, :Novelty, :Diversity, :Energy, :SeqID, :ValidAA, :KL_AA]])

    # ── Step 13: Save sequence samples ────────────────────────────────────────
    @info "Step 13: Saving results …"

    # save generated sequences to FASTA
    open(joinpath(_PATH_TO_DATA, "sa_generated_sequences.fasta"), "w") do io
        for (i, seq) in enumerate(sa_seqs)
            println(io, ">SA_sample_$(i)")
            println(io, seq)
        end
    end

    # save stored sequences
    open(joinpath(_PATH_TO_DATA, "stored_sequences.fasta"), "w") do io
        for (i, seq) in enumerate(stored_seqs)
            name = i <= length(seq_names) ? seq_names[i] : "stored_$i"
            println(io, ">$name")
            println(io, seq)
        end
    end

    # save results dataframe
    CSV.write(joinpath(_PATH_TO_DATA, "protein_experiment_results.csv"), df)

    # ── Step 14: Figures ──────────────────────────────────────────────────────
    @info "Step 14: Generating figures …"

    # Figure 1: Amino acid frequency heatmap (stored vs generated)
    function aa_freq_matrix(seqs, L)
        freq = zeros(N_AA, L)
        for seq in seqs
            for pos in 1:min(L, length(seq))
                idx = get(AA_TO_IDX, seq[pos], 0)
                idx > 0 && (freq[idx, pos] += 1)
            end
        end
        # normalize columns
        for j in 1:L
            s = sum(freq[:, j])
            s > 0 && (freq[:, j] ./= s)
        end
        return freq
    end

    freq_stored = aa_freq_matrix(stored_seqs, L)
    freq_sa     = aa_freq_matrix(sa_seqs, L)

    p1 = heatmap(freq_stored, xlabel="Position", ylabel="Amino acid",
                 yticks=(1:N_AA, string.(AA_ALPHABET)),
                 title="Stored sequences ($K patterns)",
                 color=:viridis, clims=(0, 1), size=(800, 300))
    p2 = heatmap(freq_sa, xlabel="Position", ylabel="Amino acid",
                 yticks=(1:N_AA, string.(AA_ALPHABET)),
                 title="SA-generated sequences ($S samples, β=$β_retrieval)",
                 color=:viridis, clims=(0, 1), size=(800, 300))
    p_freq = plot(p1, p2, layout=(2, 1), size=(800, 600))
    savefig(p_freq, joinpath(figpath, "Fig_protein_aa_frequencies.pdf"))
    @info "  Saved amino acid frequency heatmap"

    # Figure 2: Sequence identity distribution
    sa_ids  = [nearest_sequence_identity(s, stored_seqs) for s in sa_seqs]
    bs_ids  = [nearest_sequence_identity(s, stored_seqs) for s in bs_seqs]
    gmm_ids = [nearest_sequence_identity(s, stored_seqs) for s in gmm_seqs]

    p_ids = histogram(sa_ids, bins=30, alpha=0.6, label="SA (β=$β_retrieval)",
                      xlabel="Sequence identity to nearest stored",
                      ylabel="Count", title="Novelty in sequence space — $pfam_id",
                      color=:coral, size=(600, 400))
    histogram!(gmm_ids, bins=30, alpha=0.5, label="GMM-PCA", color=:steelblue)
    vline!([1.0], label="Exact copy", ls=:dash, color=:black, lw=2)
    savefig(p_ids, joinpath(figpath, "Fig_protein_sequence_identity.pdf"))
    @info "  Saved sequence identity histogram"

    # ── Step 15: Phase transition summary ─────────────────────────────────────
    println("\n══════════════════════════════════════════════════════")
    println("PHASE TRANSITION VALIDATION")
    println("══════════════════════════════════════════════════════")
    println("  Dimension (PCA):     d = $d")
    println("  Empirical β*:        $(round(pt.β_star, digits=2))")
    println("  Theoretical β*=√d:   $(round(pt.β_star_theory, digits=2))")
    println("  Empirical SNR*:      $(round(pt.snr_star, digits=4))")
    println("  Theoretical SNR*:    $(round(pt.snr_star_theory, digits=4))")
    println("  Ratio (empirical/theoretical β*): $(round(pt.β_star/pt.β_star_theory, digits=2))")
    println("══════════════════════════════════════════════════════")

    @info "Done."
    return df, pt, mala_mean_ar
end

# ── Run ───────────────────────────────────────────────────────────────────────
df, pt, ar = run_protein_experiment()
