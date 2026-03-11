#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Phase 4: Protein HMM Validation & Extended Protein Analysis
# Addresses Review v2 W5 / Q3:
#   - Pfam HMM validation (E-value filter)
#   - Position-specific scoring (not just marginal AA frequencies)
#   - Pairwise coupling statistics (mutual information)
#   - MALA at generation regime for protein sequences
#
# This script:
#   1. Runs the full protein experiment (all baselines + SA at multiple betas)
#   2. Saves all generated sequences as FASTA
#   3. Downloads the Pfam RRM HMM profile
#   4. Runs hmmsearch on all generated sequences
#   5. Computes position-specific and pairwise metrics
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment …"
flush(stdout); flush(stderr)
include(joinpath(@__DIR__, "Include-Sequence.jl"))
using Printf
@info "Environment loaded."
flush(stdout); flush(stderr)

# ══════════════════════════════════════════════════════════════════════════════
# Part 1: Amino acid alphabet and helpers (from run_protein_experiment.jl)
# ══════════════════════════════════════════════════════════════════════════════
const AA_ALPHABET = collect("ACDEFGHIKLMNPQRSTVWY")
const AA_TO_IDX = Dict(aa => i for (i, aa) in enumerate(AA_ALPHABET))
const N_AA = length(AA_ALPHABET)

# ── Data loading functions ────────────────────────────────────────────────────
function download_pfam_seed(pfam_id::String; cache_dir::String=_PATH_TO_DATA)
    mkpath(cache_dir)
    cache_file = joinpath(cache_dir, "$(pfam_id)_seed.sto")
    if isfile(cache_file)
        @info "  Using cached alignment: $cache_file"
        return cache_file
    end
    url = "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/$(pfam_id)/?annotation=alignment:seed"
    gz_file = cache_file * ".gz"
    @info "  Downloading seed alignment from InterPro …"
    try
        Downloads.download(url, gz_file)
        run(`gunzip -f $gz_file`)
    catch e
        isfile(gz_file) && rm(gz_file)
        error("Download failed: $e")
    end
    return cache_file
end

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
            sequences[name] *= seq
        else
            sequences[name] = seq
            push!(seq_order, name)
        end
    end
    return [(name, sequences[name]) for name in seq_order]
end

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
    col_gap_frac = [count(is_gap, char_mat[:, j]) / K_raw for j in 1:L_raw]
    keep_cols = findall(f -> f <= max_gap_frac_col, col_gap_frac)
    char_mat = char_mat[:, keep_cols]
    L = length(keep_cols)
    seq_gap_frac = [count(is_gap, char_mat[i, :]) / L for i in 1:K_raw]
    keep_seqs = findall(f -> f <= max_gap_frac_seq, seq_gap_frac)
    char_mat = char_mat[keep_seqs, :]
    names = names[keep_seqs]
    @info "  Alignment: $(K_raw) seqs × $(L_raw) cols → $(length(keep_seqs)) seqs × $L cols"
    return char_mat, names
end

# ── Encoding/decoding ─────────────────────────────────────────────────────────
function onehot_encode(char_mat::Matrix{Char})
    K, L = size(char_mat)
    d_full = N_AA * L
    X = zeros(Float64, d_full, K)
    for k in 1:K, pos in 1:L
        idx = get(AA_TO_IDX, char_mat[k, pos], 0)
        idx > 0 && (X[(pos-1)*N_AA + idx, k] = 1.0)
    end
    return X
end

function decode_onehot(x::Vector{Float64}, L::Int)
    seq = Char[]
    for pos in 1:L
        block = x[(pos-1)*N_AA+1 : pos*N_AA]
        best_idx = argmax(block)
        push!(seq, maximum(block) < 1e-10 ? '-' : AA_ALPHABET[best_idx])
    end
    return String(seq)
end

# ── Protein metrics ───────────────────────────────────────────────────────────
function sequence_identity(seq1::String, seq2::String)
    L = min(length(seq1), length(seq2))
    matches = 0; compared = 0
    for i in 1:L
        (seq1[i] == '-' || seq2[i] == '-') && continue
        compared += 1
        seq1[i] == seq2[i] && (matches += 1)
    end
    return compared > 0 ? matches / compared : 0.0
end

function nearest_sequence_identity(gen_seq::String, stored_seqs::Vector{String})
    return maximum(sequence_identity(gen_seq, s) for s in stored_seqs)
end

function valid_residue_fraction(seq::String)
    non_gap = count(c -> c != '-', seq)
    non_gap == 0 && return 0.0
    return count(c -> c in AA_ALPHABET, seq) / non_gap
end

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
    p = aa_freqs(stored_seqs); q = aa_freqs(gen_seqs)
    eps = 1e-10
    p .+= eps; p ./= sum(p)
    q .+= eps; q ./= sum(q)
    return sum(p[i] * log(p[i] / q[i]) for i in 1:N_AA)
end

# ══════════════════════════════════════════════════════════════════════════════
# Part 2: NEW — Position-specific and pairwise metrics
# ══════════════════════════════════════════════════════════════════════════════

"""
    position_specific_kl(gen_seqs, stored_seqs, L) -> (mean_kl, per_position_kl)

KL divergence computed per-position, then averaged.
This captures position-specific conservation patterns, not just global composition.
"""
function position_specific_kl(gen_seqs::Vector{String}, stored_seqs::Vector{String}, L::Int)
    eps = 1e-10
    per_pos_kl = zeros(L)

    for pos in 1:L
        # count AA frequencies at this position
        p_counts = zeros(N_AA)  # stored (reference)
        q_counts = zeros(N_AA)  # generated
        for seq in stored_seqs
            pos <= length(seq) || continue
            idx = get(AA_TO_IDX, seq[pos], 0)
            idx > 0 && (p_counts[idx] += 1)
        end
        for seq in gen_seqs
            pos <= length(seq) || continue
            idx = get(AA_TO_IDX, seq[pos], 0)
            idx > 0 && (q_counts[idx] += 1)
        end

        # normalize
        p = p_counts .+ eps; p ./= sum(p)
        q = q_counts .+ eps; q ./= sum(q)

        per_pos_kl[pos] = sum(p[i] * log(p[i] / q[i]) for i in 1:N_AA)
    end

    return mean(per_pos_kl), per_pos_kl
end

"""
    pairwise_mutual_information(seqs, L) -> Matrix{Float64}

Compute mutual information between all position pairs.
MI(i,j) = sum_{a,b} p(a,b) log(p(a,b) / (p(a)*p(b)))
This captures pairwise coupling (co-evolution) patterns.
"""
function pairwise_mutual_information(seqs::Vector{String}, L::Int; n_positions::Int=20)
    # subsample positions for computational tractability
    positions = round.(Int, range(1, L, length=min(n_positions, L)))
    np = length(positions)
    MI = zeros(np, np)
    eps = 1e-10
    n_seqs = length(seqs)

    for (ii, pi) in enumerate(positions), (jj, pj) in enumerate(positions)
        ii >= jj && continue  # upper triangle only

        # joint and marginal counts
        joint = zeros(N_AA, N_AA)
        marg_i = zeros(N_AA)
        marg_j = zeros(N_AA)

        for seq in seqs
            (pi > length(seq) || pj > length(seq)) && continue
            ai = get(AA_TO_IDX, seq[pi], 0)
            aj = get(AA_TO_IDX, seq[pj], 0)
            (ai == 0 || aj == 0) && continue
            joint[ai, aj] += 1
            marg_i[ai] += 1
            marg_j[aj] += 1
        end

        total = sum(joint)
        total < 2 && continue

        # normalize
        p_joint = (joint .+ eps) ./ (total + N_AA^2 * eps)
        p_i = (marg_i .+ eps) ./ (sum(marg_i) + N_AA * eps)
        p_j = (marg_j .+ eps) ./ (sum(marg_j) + N_AA * eps)

        mi = 0.0
        for a in 1:N_AA, b in 1:N_AA
            mi += p_joint[a, b] * log(p_joint[a, b] / (p_i[a] * p_j[b]))
        end
        MI[ii, jj] = mi
        MI[jj, ii] = mi
    end

    return MI, positions
end

"""
    mi_correlation(MI_gen, MI_stored) -> Float64

Pearson correlation between upper-triangle MI matrices.
Measures how well generated sequences reproduce pairwise coupling patterns.
"""
function mi_correlation(MI_gen::Matrix{Float64}, MI_stored::Matrix{Float64})
    n = size(MI_gen, 1)
    gen_vals = Float64[]
    stored_vals = Float64[]
    for i in 1:n, j in (i+1):n
        push!(gen_vals, MI_gen[i, j])
        push!(stored_vals, MI_stored[i, j])
    end
    length(gen_vals) < 3 && return 0.0
    return cor(gen_vals, stored_vals)
end

# ══════════════════════════════════════════════════════════════════════════════
# Part 3: NEW — HMMER validation
# ══════════════════════════════════════════════════════════════════════════════

"""
    download_pfam_hmm(pfam_id; cache_dir) -> String

Download the Pfam HMM profile for hmmsearch validation.
"""
function download_pfam_hmm(pfam_id::String; cache_dir::String=_PATH_TO_DATA)
    mkpath(cache_dir)
    hmm_file = joinpath(cache_dir, "$(pfam_id).hmm")

    if isfile(hmm_file)
        @info "  Using cached HMM: $hmm_file"
        return hmm_file
    end

    # Try InterPro API for HMM
    url = "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/$(pfam_id)/?annotation=hmm"
    gz_file = hmm_file * ".gz"
    @info "  Downloading HMM profile from InterPro …"

    try
        Downloads.download(url, gz_file)
        run(`gunzip -f $gz_file`)
        @info "  Saved HMM to $hmm_file"
    catch e
        isfile(gz_file) && rm(gz_file)
        @warn "  InterPro download failed, trying Pfam legacy …"

        # Fallback: try Pfam legacy URL
        url2 = "https://pfam.xfam.org/family/$(pfam_id)/hmm"
        try
            Downloads.download(url2, hmm_file)
            @info "  Saved HMM to $hmm_file (legacy URL)"
        catch e2
            @error "  Could not download HMM profile. Please download manually:"
            @error "  Place the HMM at: $hmm_file"
            error("Cannot proceed without HMM profile.")
        end
    end

    return hmm_file
end

"""
    write_fasta(filepath, seqs, prefix)

Write sequences to FASTA format.
"""
function write_fasta(filepath::String, seqs::Vector{String}, prefix::String)
    open(filepath, "w") do io
        for (i, seq) in enumerate(seqs)
            # Remove gaps for hmmsearch (it expects unaligned sequences)
            clean_seq = replace(seq, "-" => "", "." => "")
            if length(clean_seq) > 0
                println(io, ">$(prefix)_$(i)")
                println(io, clean_seq)
            end
        end
    end
end

"""
    run_hmmsearch(hmm_file, fasta_file; evalue_threshold) -> (n_hits, n_total, evalues)

Run hmmsearch and return hit statistics.
"""
function run_hmmsearch(hmm_file::String, fasta_file::String;
                       evalue_threshold::Float64=0.01)

    # run hmmsearch with table output
    tblout = fasta_file * ".hmmsearch.tbl"
    try
        run(pipeline(`hmmsearch --tblout $tblout -E 1000 --noali $hmm_file $fasta_file`,
                     devnull))
    catch e
        @warn "hmmsearch failed: $e"
        return (0, 0, Float64[])
    end

    # parse tblout
    evalues = Float64[]
    hit_names = String[]
    for line in eachline(tblout)
        startswith(line, "#") && continue
        parts = split(line)
        length(parts) >= 5 || continue
        try
            ev = parse(Float64, parts[5])  # E-value is column 5 in tblout
            push!(evalues, ev)
            push!(hit_names, parts[1])
        catch
            continue
        end
    end

    # count sequences in FASTA
    n_total = count(line -> startswith(line, ">"), readlines(fasta_file))

    # hits passing threshold
    n_hits = count(ev -> ev < evalue_threshold, evalues)

    return (n_hits=n_hits, n_total=n_total, evalues=evalues, hit_names=hit_names)
end

# ══════════════════════════════════════════════════════════════════════════════
# Part 4: Entropy inflection + GMM helpers (from original)
# ══════════════════════════════════════════════════════════════════════════════
function find_entropy_inflection(X̂::Matrix{Float64};
                                  α::Float64=0.01, n_betas::Int=50,
                                  β_range::Tuple{Float64,Float64}=(0.1, 500.0))
    d, K = size(X̂)
    βs = 10 .^ range(log10(β_range[1]), log10(β_range[2]), length=n_betas)
    n_probes = min(K, 20)
    Hs = zeros(n_betas)
    for (bi, β) in enumerate(βs)
        H_sum = 0.0
        for k in 1:n_probes
            H_sum += attention_entropy(X̂[:, k], X̂, β)
        end
        Hs[bi] = H_sum / n_probes
    end
    log_βs = log.(βs)
    dH = diff(Hs) ./ diff(log_βs)
    d2H = diff(dH) ./ diff(log_βs[1:end-1])
    inflection_idx = argmin(d2H) + 1
    β_star = βs[inflection_idx]
    snr_star = sqrt(α * β_star / (2 * d))
    β_star_theory = sqrt(d)
    snr_star_theory = sqrt(α / (2 * sqrt(d)))
    @info "  Phase transition: β*=$(round(β_star, digits=2)) (empirical), √d=$(round(β_star_theory, digits=2)) (theory)"
    return (β_star=β_star, snr_star=snr_star, β_star_theory=β_star_theory,
            snr_star_theory=snr_star_theory, βs=βs, Hs=Hs)
end

struct _ProteinVAE
    enc_shared; enc_μ; enc_logσ²; dec
end
Flux.@layer _ProteinVAE

function _logsumexp_cols(A::Matrix{Float64})
    m = maximum(A, dims=1)
    return m .+ log.(sum(exp.(A .- m), dims=1))
end

function _fit_gmm(X::Matrix{Float64}, C::Int; n_iters=500, tol=1e-8, seed=42)
    r, N = size(X)
    rng  = MersenneTwister(seed)
    μ = zeros(r, C)
    μ[:, 1] = X[:, rand(rng, 1:N)]
    for k in 2:C
        dists = [minimum(sum((X[:, n] .- μ[:, j]).^2) for j in 1:(k-1)) for n in 1:N]
        μ[:, k] = X[:, StatsBase.sample(rng, 1:N, Weights(dists ./ sum(dists)))]
    end
    log_σ² = log.(repeat(max.(vec(var(X, dims=2)), 1e-6), 1, C))
    π_k = ones(C) / C
    log_r = zeros(C, N)
    log_lik_prev = -Inf
    for _ in 1:n_iters
        σ² = exp.(log_σ²)
        for k in 1:C
            diff = X .- μ[:, k]
            log_r[k, :] = log(π_k[k]) .-
                           0.5 .* vec(sum(diff.^2 ./ σ²[:, k], dims=1)) .-
                           0.5 .* sum(log.(2π .* σ²[:, k]))
        end
        lse = _logsumexp_cols(log_r)
        log_lik = sum(lse)
        log_r .-= lse
        r_mat = exp.(log_r)
        abs(log_lik - log_lik_prev) / (abs(log_lik_prev) + 1e-30) < tol && break
        log_lik_prev = log_lik
        N_k = vec(sum(r_mat, dims=2))
        π_k = max.(N_k / N, 1e-8); π_k ./= sum(π_k)
        for k in 1:C
            μ[:, k] = X * r_mat[k, :] / N_k[k]
            diff = X .- μ[:, k]
            log_σ²[:, k] = log.(max.(vec(sum(diff.^2 .* r_mat[k, :]', dims=2)) / N_k[k], 1e-6))
        end
    end
    return μ, log_σ², π_k
end

function _sample_gmm(μ, log_σ², π_k, n_samples; seed=9999)
    rng = MersenneTwister(seed)
    r = size(μ, 1)
    wts = Weights(π_k)
    [let k = StatsBase.sample(rng, 1:length(π_k), wts)
         μ[:, k] .+ sqrt.(exp.(log_σ²[:, k])) .* randn(rng, r)
     end for _ in 1:n_samples]
end

# ══════════════════════════════════════════════════════════════════════════════
# Part 5: MAIN EXPERIMENT
# ══════════════════════════════════════════════════════════════════════════════

const PFAM_ID = "PF00076"
const K_MAX = 100
const α_step = 0.01
const S = 150
const n_chains = 30
const T_per_chain = 5000
const T_burnin = 2000
const thin_interval = 100
const samples_per_chain = 5
const σ_init = 0.01

function main()
    figpath = _PATH_TO_FIG
    datapath = _PATH_TO_DATA
    mkpath(figpath); mkpath(datapath)

    # ── Step 1: Load alignment ─────────────────────────────────────────────────
    @info "Step 1: Loading Pfam $PFAM_ID alignment …"
    flush(stdout); flush(stderr)
    sto_file = download_pfam_seed(PFAM_ID)
    raw_seqs = parse_stockholm(sto_file)
    isempty(raw_seqs) && (raw_seqs = parse_fasta(sto_file))
    @info "  Parsed $(length(raw_seqs)) sequences"

    char_mat, seq_names = clean_alignment(raw_seqs)
    K_total, L = size(char_mat)
    if K_total > K_MAX
        Random.seed!(42)
        keep = StatsBase.sample(1:K_total, K_MAX, replace=false) |> sort
        char_mat = char_mat[keep, :]
        seq_names = seq_names[keep]
    end
    K = size(char_mat, 1)
    stored_seqs = [String(char_mat[k, :]) for k in 1:K]
    @info "  Using K=$K sequences, L=$L positions"
    flush(stdout); flush(stderr)

    # ── Step 2: Encode → PCA → normalize ───────────────────────────────────────
    @info "Step 2: Encoding …"
    X_onehot = onehot_encode(char_mat)
    pca_model = MultivariateStats.fit(PCA, X_onehot; pratio=0.95)
    d_pca = outdim(pca_model)
    Z = MultivariateStats.transform(pca_model, X_onehot)

    ϵ = 1e-12
    X̂ = copy(Z)
    for k in 1:K
        X̂[:, k] ./= (norm(X̂[:, k]) + ϵ)
    end
    d = size(X̂, 1)
    @info "  PCA: $(size(X_onehot,1)) → $d dimensions"

    decode_sample(ξ) = decode_onehot(vec(MultivariateStats.reconstruct(pca_model, ξ)), L)

    # ── Step 3: Phase transition ───────────────────────────────────────────────
    @info "Step 3: Phase transition analysis …"
    pt = find_entropy_inflection(X̂; α=α_step)
    β_ret = Float64(round(Int, 20 * pt.β_star))
    β_gen = Float64(round(Int, 2 * pt.β_star))
    @info "  β_ret=$β_ret, β_gen=$β_gen"
    flush(stdout); flush(stderr)

    # ── Step 4: Run all samplers ───────────────────────────────────────────────
    function run_multichain(X̂, β; use_mala=false)
        samples = Vector{Vector{Float64}}()
        accept_rates = Float64[]
        Random.seed!(42)
        pattern_indices = StatsBase.sample(1:size(X̂,2), n_chains, replace=(n_chains > size(X̂,2)))
        for (c, k) in enumerate(pattern_indices)
            Random.seed!(12345 + c)
            d_loc = size(X̂, 1)
            sₒ = X̂[:, k] .+ σ_init .* randn(d_loc)
            if use_mala
                (_, Ξ, ar) = mala_sample(X̂, sₒ, T_per_chain; β=β, α=α_step, seed=12345+c)
                push!(accept_rates, ar)
            else
                (_, Ξ) = sample(X̂, sₒ, T_per_chain; β=β, α=α_step, seed=12345+c)
            end
            chain_pool = Vector{Vector{Float64}}()
            for tᵢ in (T_burnin+1):thin_interval:T_per_chain
                push!(chain_pool, Ξ[tᵢ, :])
            end
            n_avail = length(chain_pool)
            idxs = round.(Int, range(1, n_avail, length=min(samples_per_chain, n_avail)))
            for idx in idxs
                push!(samples, chain_pool[idx])
            end
        end
        ar_mean = isempty(accept_rates) ? NaN : mean(accept_rates)
        return samples, ar_mean
    end

    @info "Step 4: Running samplers …"
    flush(stdout); flush(stderr)

    @info "  SA (retrieval, β=$β_ret) …"
    sa_ret_samps, _ = run_multichain(X̂, β_ret)
    @info "  SA (generation, β=$β_gen) …"
    sa_gen_samps, _ = run_multichain(X̂, β_gen)
    @info "  MALA (retrieval, β=$β_ret) …"
    mala_ret_samps, mala_ret_ar = run_multichain(X̂, β_ret; use_mala=true)
    @info "  MALA (generation, β=$β_gen) …"
    mala_gen_samps, mala_gen_ar = run_multichain(X̂, β_gen; use_mala=true)

    # Baselines
    @info "  Baselines …"
    bs_samps = [copy(X̂[:, rand(1:K)]) for _ in 1:S]
    Random.seed!(12345)

    σ_noise = sqrt(2 * α_step / β_ret)
    gp_samps = [X̂[:, rand(1:K)] .+ σ_noise .* randn(d) for _ in 1:S]

    dirichlet_dist = Dirichlet(K, 1.0)
    Random.seed!(12345)
    rc_samps = [X̂ * rand(dirichlet_dist) for _ in 1:S]

    # GMM-PCA
    R_gmm = min(50, d - 1)
    pca_gmm = MultivariateStats.fit(PCA, X̂; maxoutdim=R_gmm)
    Z_gmm = MultivariateStats.transform(pca_gmm, X̂)
    C_actual = min(10, K - 1)
    μ_gmm, log_σ²_gmm, π_gmm = _fit_gmm(Z_gmm, C_actual; seed=42)
    pca_s = _sample_gmm(μ_gmm, log_σ²_gmm, π_gmm, S; seed=7777)
    gmm_samps = [vec(MultivariateStats.reconstruct(pca_gmm, z)) for z in pca_s]

    # VAE
    @info "  Training VAE …"
    _h1 = max(d ÷ 2, 16); _h2 = max(d ÷ 4, 8)
    pvae = _ProteinVAE(
        Chain(Dense(d => _h1, relu), Dense(_h1 => _h2, relu)),
        Dense(_h2 => 8), Dense(_h2 => 8),
        Chain(Dense(8 => _h2, relu), Dense(_h2 => _h1, relu), Dense(_h1 => d))
    )
    X_t = Float32.(X̂)
    opt_vae = Flux.setup(Adam(1f-3), pvae)
    for epoch in 1:2000
        ε = randn(Float32, 8, K)
        _, grads = Flux.withgradient(pvae) do m
            h = m.enc_shared(X_t); μ = m.enc_μ(h); lσ² = m.enc_logσ²(h)
            z = μ .+ exp.(0.5f0 .* lσ²) .* ε
            o = m.dec(z); x̂ = o ./ (sqrt.(sum(o .^ 2; dims=1)) .+ 1f-8)
            mean(sum((X_t .- x̂) .^ 2; dims=1))
        end
        Flux.update!(opt_vae, pvae, grads[1])
    end
    for epoch in 1:2000
        kl_w = 0.0001f0 * Float32(epoch) / 2000f0
        ε = randn(Float32, 8, K)
        _, grads = Flux.withgradient(pvae) do m
            h = m.enc_shared(X_t); μ = m.enc_μ(h); lσ² = m.enc_logσ²(h)
            z = μ .+ exp.(0.5f0 .* lσ²) .* ε
            o = m.dec(z); x̂ = o ./ (sqrt.(sum(o .^ 2; dims=1)) .+ 1f-8)
            recon = mean(sum((X_t .- x̂) .^ 2; dims=1))
            kl = -0.5f0 * mean(sum(1f0 .+ lσ² .- μ .^ 2 .- exp.(lσ²); dims=1))
            recon + kl_w * kl
        end
        Flux.update!(opt_vae, pvae, grads[1])
    end
    Random.seed!(9999)
    raw_vae = let o = pvae.dec(randn(Float32, 8, S)); o ./ (sqrt.(sum(o .^ 2; dims=1)) .+ 1f-8) end
    vae_samps = [Float64.(raw_vae[:, i]) for i in 1:S]

    flush(stdout); flush(stderr)

    # ── Step 5: Decode all to amino acid sequences ─────────────────────────────
    @info "Step 5: Decoding to amino acid sequences …"
    all_methods = [
        ("Bootstrap",           bs_samps),
        ("Gaussian_perturb",    gp_samps),
        ("Convex_combination",  rc_samps),
        ("GMM_PCA",             gmm_samps),
        ("VAE_lat8",            vae_samps),
        ("MALA_ret_b$(Int(β_ret))",  mala_ret_samps),
        ("MALA_gen_b$(Int(β_gen))",  mala_gen_samps),
        ("SA_ret_b$(Int(β_ret))",    sa_ret_samps),
        ("SA_gen_b$(Int(β_gen))",    sa_gen_samps),
    ]

    decoded = Dict{String, Vector{String}}()
    for (name, samps) in all_methods
        decoded[name] = [decode_sample(ξ) for ξ in samps]
    end
    flush(stdout); flush(stderr)

    # ── Step 6: Write FASTA files for HMMER ────────────────────────────────────
    @info "Step 6: Writing FASTA files …"
    fasta_dir = joinpath(datapath, "hmm_validation")
    mkpath(fasta_dir)

    # Write stored sequences
    write_fasta(joinpath(fasta_dir, "stored.fasta"), stored_seqs, "stored")

    # Write each method's sequences
    fasta_files = Dict{String, String}()
    for (name, _) in all_methods
        fpath = joinpath(fasta_dir, "$(name).fasta")
        write_fasta(fpath, decoded[name], name)
        fasta_files[name] = fpath
    end
    flush(stdout); flush(stderr)

    # ── Step 7: Download HMM and run hmmsearch ─────────────────────────────────
    @info "Step 7: HMMER validation …"
    flush(stdout); flush(stderr)
    hmm_file = download_pfam_hmm(PFAM_ID)

    println("\n" * "="^80)
    println("HMMER VALIDATION RESULTS (Pfam $PFAM_ID RRM)")
    println("E-value threshold: 0.01")
    println("="^80)

    hmm_results = Dict{String, NamedTuple}()

    # First validate stored sequences
    stored_fasta = joinpath(fasta_dir, "stored.fasta")
    stored_hmm = run_hmmsearch(hmm_file, stored_fasta)
    @printf("%-30s | Hits: %3d / %3d (%.1f%%) | Median E-value: %.2e\n",
        "Stored (reference)", stored_hmm.n_hits, stored_hmm.n_total,
        100.0 * stored_hmm.n_hits / max(stored_hmm.n_total, 1),
        length(stored_hmm.evalues) > 0 ? median(stored_hmm.evalues) : NaN)

    for (name, _) in all_methods
        result = run_hmmsearch(hmm_file, fasta_files[name])
        hmm_results[name] = result
        n_t = max(result.n_total, 1)
        med_ev = length(result.evalues) > 0 ? median(result.evalues) : NaN
        @printf("%-30s | Hits: %3d / %3d (%.1f%%) | Median E-value: %.2e\n",
            name, result.n_hits, result.n_total,
            100.0 * result.n_hits / n_t, med_ev)
    end
    flush(stdout)

    # ── Step 8: Position-specific KL ───────────────────────────────────────────
    @info "Step 8: Position-specific KL divergence …"
    println("\n" * "="^80)
    println("POSITION-SPECIFIC KL DIVERGENCE (per-position, then averaged)")
    println("="^80)

    for (name, _) in all_methods
        mean_kl, _ = position_specific_kl(decoded[name], stored_seqs, L)
        @printf("%-30s | Mean per-position KL: %.4f\n", name, mean_kl)
    end
    flush(stdout)

    # ── Step 9: Pairwise coupling (mutual information) ─────────────────────────
    @info "Step 9: Pairwise mutual information analysis …"
    println("\n" * "="^80)
    println("PAIRWISE MUTUAL INFORMATION CORRELATION")
    println("(Pearson r between MI matrices: generated vs stored)")
    println("="^80)

    MI_stored, mi_positions = pairwise_mutual_information(stored_seqs, L)

    for (name, _) in all_methods
        MI_gen, _ = pairwise_mutual_information(decoded[name], L)
        r_mi = mi_correlation(MI_gen, MI_stored)
        @printf("%-30s | MI correlation: %.4f\n", name, r_mi)
    end
    flush(stdout)

    # ── Step 10: Global amino acid KL ──────────────────────────────────────────
    println("\n" * "="^80)
    println("GLOBAL AMINO ACID COMPOSITION KL + SEQUENCE IDENTITY")
    println("="^80)

    for (name, _) in all_methods
        kl = aa_composition_kl(decoded[name], stored_seqs)
        seq_ids = [nearest_sequence_identity(s, stored_seqs) for s in decoded[name]]
        valid = mean(valid_residue_fraction(s) for s in decoded[name])
        @printf("%-30s | KL=%.4f | SeqID=%.3f±%.3f | ValidAA=%.3f\n",
            name, kl, mean(seq_ids), std(seq_ids), valid)
    end
    flush(stdout)

    # ── Step 11: MALA comparison at generation regime ──────────────────────────
    println("\n" * "="^80)
    println("MALA vs SA AT GENERATION REGIME (beta_gen=$(Int(β_gen)))")
    println("="^80)
    println("MALA acceptance rate (ret, beta=$(Int(β_ret))): $(round(mala_ret_ar, digits=4))")
    println("MALA acceptance rate (gen, beta=$(Int(β_gen))): $(round(mala_gen_ar, digits=4))")

    for (name, samps) in [("SA_gen", sa_gen_samps), ("MALA_gen", mala_gen_samps)]
        seqs = [decode_sample(ξ) for ξ in samps]
        kl = aa_composition_kl(seqs, stored_seqs)
        seq_ids = [nearest_sequence_identity(s, stored_seqs) for s in seqs]
        nov = mean(sample_novelty(ξ, X̂) for ξ in samps)
        div = sample_diversity(samps)
        @printf("%-15s | KL=%.4f | SeqID=%.3f | Nov=%.4f | Div=%.4f\n",
            name, kl, mean(seq_ids), nov, div)
    end
    flush(stdout)

    # ── Step 12: Combined summary table ────────────────────────────────────────
    println("\n" * "="^80)
    println("COMBINED SUMMARY TABLE")
    println("="^80)
    println("Method                         | AA_KL  | Pos_KL | MI_corr | HMM_pass | SeqID  | Nov    | Div")
    println("-"^110)

    for (name, samps) in all_methods
        seqs = decoded[name]
        kl_aa = aa_composition_kl(seqs, stored_seqs)
        mean_pos_kl, _ = position_specific_kl(seqs, stored_seqs, L)
        MI_gen, _ = pairwise_mutual_information(seqs, L)
        r_mi = mi_correlation(MI_gen, MI_stored)
        hmm = hmm_results[name]
        hmm_pct = 100.0 * hmm.n_hits / max(hmm.n_total, 1)
        seq_ids = [nearest_sequence_identity(s, stored_seqs) for s in seqs]
        nov = mean(sample_novelty(ξ, X̂) for ξ in samps)
        div = sample_diversity(samps)

        @printf("%-30s | %.4f | %.4f | %.4f  | %5.1f%%  | %.3f  | %.4f | %.4f\n",
            name, kl_aa, mean_pos_kl, r_mi, hmm_pct, mean(seq_ids), nov, div)
    end
    flush(stdout)

    println("\n" * "="^80)
    println("EXPERIMENT COMPLETE")
    println("="^80)
    println("Data saved to: $fasta_dir")
    println("Key files: stored.fasta, SA_*.fasta, VAE_*.fasta")
end

main()
