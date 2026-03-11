#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Pfam Families Experiment (Paper Section 3a)
#
# Run stochastic attention on multiple Pfam families with varying size,
# conservation, and sequence length. Compare against simple baselines.
# Protein-specific learned baselines (DeepSequence, EvoDiff, etc.) are
# handled separately in ../baselines/.
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment …"
const _EXPERIMENT_ROOT = @__DIR__
include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))

# override paths for this experiment (Include.jl sets PFAM_DATA_DIR and PFAM_FIG_DIR
# for the top-level code/ directory; we want experiment-local paths)
const PFAM_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "data")
const PFAM_FIG_DIR  = joinpath(_EXPERIMENT_ROOT, "figs")
@info "Environment loaded."

# ══════════════════════════════════════════════════════════════════════════════
# Family definitions
# ══════════════════════════════════════════════════════════════════════════════

const FAMILIES = [
    (id="PF00076", name="RRM",     desc="RNA Recognition Motif"),
    (id="PF00018", name="SH3",     desc="SH3 domain"),
    (id="PF00397", name="WW",      desc="WW domain"),
    (id="PF00014", name="Kunitz",  desc="Kunitz/BPTI domain"),
    (id="PF00096", name="zf-C2H2", desc="Zinc finger C2H2"),
    (id="PF00595", name="PDZ",     desc="PDZ domain"),
    (id="PF00069", name="Pkinase", desc="Protein kinase domain"),
]

# ══════════════════════════════════════════════════════════════════════════════
# Experiment parameters
# ══════════════════════════════════════════════════════════════════════════════

const α_step        = 0.01       # Langevin step size
const n_chains      = 30         # number of independent chains
const T_per_chain   = 5000       # iterations per chain
const T_burnin      = 2000       # burn-in iterations to discard
const thin_interval = 100        # thinning interval
const samples_per_chain = 5      # samples to keep per chain
const S             = n_chains * samples_per_chain  # total samples (150)
const σ_init        = 0.01       # initialization noise
const R_pca         = 0.95       # fraction of variance retained in PCA
const K_MAX         = 500        # max sequences per family (subsample if larger)

# ══════════════════════════════════════════════════════════════════════════════
# Multi-chain SA sampler (shared logic)
# ══════════════════════════════════════════════════════════════════════════════

"""
    run_sa_chains(X̂, β; n_chains, T, T_burn, thin, spc, σ_init, α) -> Vector{Vector{Float64}}

Run multi-chain stochastic attention sampling on memory matrix X̂.
Returns a flat vector of PCA-space samples.
"""
function run_sa_chains(X̂::Matrix{Float64}, β::Float64;
                       n_chains::Int=n_chains, T::Int=T_per_chain,
                       T_burn::Int=T_burnin, thin::Int=thin_interval,
                       spc::Int=samples_per_chain, σ_init::Float64=σ_init,
                       α::Float64=α_step)
    d, K = size(X̂)
    all_samples = Vector{Vector{Float64}}()

    Random.seed!(42)
    pattern_indices = StatsBase.sample(1:K, n_chains, replace=(n_chains > K))

    for (c, k) in enumerate(pattern_indices)
        seed_c = 12345 + c
        Random.seed!(seed_c)
        ξ₀ = X̂[:, k] .+ σ_init .* randn(d)
        (_, Ξ) = sample(X̂, ξ₀, T; β=β, α=α, seed=seed_c)

        # collect post-burn-in, thinned samples
        chain_pool = Vector{Vector{Float64}}()
        for tᵢ in (T_burn+1):thin:T
            push!(chain_pool, Ξ[tᵢ, :])
        end

        # evenly space spc samples from the pool
        n_avail = length(chain_pool)
        n_avail == 0 && continue
        idxs = round.(Int, range(1, n_avail, length=min(spc, n_avail)))
        for idx in idxs
            push!(all_samples, chain_pool[idx])
        end
    end

    return all_samples
end

"""
    run_mala_chains(X̂, β; ...) -> (Vector{Vector{Float64}}, Float64)

Run multi-chain MALA. Returns (samples, mean_accept_rate).
"""
function run_mala_chains(X̂::Matrix{Float64}, β::Float64;
                         n_chains::Int=n_chains, T::Int=T_per_chain,
                         T_burn::Int=T_burnin, thin::Int=thin_interval,
                         spc::Int=samples_per_chain, σ_init::Float64=σ_init,
                         α::Float64=α_step)
    d, K = size(X̂)
    all_samples = Vector{Vector{Float64}}()
    accept_rates = Float64[]

    Random.seed!(42)
    pattern_indices = StatsBase.sample(1:K, n_chains, replace=(n_chains > K))

    for (c, k) in enumerate(pattern_indices)
        seed_c = 12345 + c
        Random.seed!(seed_c)
        ξ₀ = X̂[:, k] .+ σ_init .* randn(d)
        (_, Ξ, ar) = mala_sample(X̂, ξ₀, T; β=β, α=α, seed=seed_c)
        push!(accept_rates, ar)

        chain_pool = Vector{Vector{Float64}}()
        for tᵢ in (T_burn+1):thin:T
            push!(chain_pool, Ξ[tᵢ, :])
        end

        n_avail = length(chain_pool)
        n_avail == 0 && continue
        idxs = round.(Int, range(1, n_avail, length=min(spc, n_avail)))
        for idx in idxs
            push!(all_samples, chain_pool[idx])
        end
    end

    return all_samples, mean(accept_rates)
end

# ══════════════════════════════════════════════════════════════════════════════
# Simple baselines
# ══════════════════════════════════════════════════════════════════════════════

"""
    run_simple_baselines(X̂, β, S) -> Dict{String, Vector{Vector{Float64}}}

Run bootstrap, Gaussian perturbation, and random convex combination baselines.
"""
function run_simple_baselines(X̂::Matrix{Float64}, β::Float64, S::Int)
    d, K = size(X̂)
    σ_noise = sqrt(2 * α_step / β)

    # Bootstrap (replay)
    Random.seed!(12345)
    bs = [copy(X̂[:, rand(1:K)]) for _ in 1:S]

    # Gaussian perturbation
    Random.seed!(12345)
    gp = [X̂[:, rand(1:K)] .+ σ_noise .* randn(d) for _ in 1:S]

    # Random convex combination
    dirichlet = Dirichlet(K, 1.0)
    Random.seed!(12345)
    rc = [X̂ * rand(dirichlet) for _ in 1:S]

    return Dict("Bootstrap" => bs, "Gaussian perturbation" => gp, "Convex combination" => rc)
end

# ══════════════════════════════════════════════════════════════════════════════
# Evaluation
# ══════════════════════════════════════════════════════════════════════════════

"""
    evaluate_method(samps, seqs, X̂, β, stored_seqs, L) -> NamedTuple

Compute all metrics for a set of generated samples/sequences.
"""
function evaluate_method(samps::Vector{Vector{Float64}}, seqs::Vector{String},
                         X̂::Matrix{Float64}, β::Float64,
                         stored_seqs::Vector{String}, L::Int)

    novelty_vals  = [sample_novelty(ξ, X̂) for ξ in samps]
    energy_vals   = [hopfield_energy(ξ, X̂, β) for ξ in samps]
    seq_id_vals   = [nearest_sequence_identity(s, stored_seqs) for s in seqs]

    return (
        Novelty   = mean(novelty_vals),
        Diversity = sample_diversity(samps),
        Energy    = mean(energy_vals),
        SeqID     = mean(seq_id_vals),
        ValidAA   = mean(valid_residue_fraction(s) for s in seqs),
        KL_AA     = aa_composition_kl(seqs, stored_seqs),
    )
end

# ══════════════════════════════════════════════════════════════════════════════
# Main: run experiment for each family
# ══════════════════════════════════════════════════════════════════════════════

function run_pfam_experiment()

    mkpath(PFAM_DATA_DIR)
    mkpath(PFAM_FIG_DIR)

    # master results table
    all_rows = Vector{NamedTuple}()

    # family-level summary
    family_summary = Vector{NamedTuple}()

    for fam in FAMILIES
        @info "════════════════════════════════════════════════════════════"
        @info "  Family: $(fam.name) ($(fam.id)) — $(fam.desc)"
        @info "════════════════════════════════════════════════════════════"

        # -- Step 1: Download and parse --
        fam_data_dir = joinpath(PFAM_DATA_DIR, fam.id)
        mkpath(fam_data_dir)

        @info "Step 1: Loading alignment …"
        sto_file = download_pfam_seed(fam.id; cache_dir=fam_data_dir)
        raw_seqs = parse_stockholm(sto_file)
        if isempty(raw_seqs)
            @info "  No Stockholm sequences, trying FASTA …"
            raw_seqs = parse_fasta(sto_file)
        end
        @info "  Parsed $(length(raw_seqs)) sequences"

        if isempty(raw_seqs)
            @warn "  Skipping $(fam.name): no sequences parsed"
            continue
        end

        # -- Step 2: Clean alignment --
        @info "Step 2: Cleaning alignment …"
        char_mat, seq_names = clean_alignment(raw_seqs)
        K_total, L = size(char_mat)

        if K_total > K_MAX
            Random.seed!(42)
            keep = StatsBase.sample(1:K_total, K_MAX, replace=false) |> sort
            char_mat = char_mat[keep, :]
            seq_names = seq_names[keep]
            @info "  Subsampled to $K_MAX sequences"
        end
        K = size(char_mat, 1)
        stored_seqs = [String(char_mat[k, :]) for k in 1:K]

        # -- Step 3: Build memory matrix --
        @info "Step 3: Building memory matrix …"
        X̂, pca_model, L, d_full = build_memory_matrix(char_mat; pratio=R_pca)
        d = size(X̂, 1)

        # -- Step 4: MSA statistics (for β* prediction) --
        @info "Step 4: Computing MSA statistics …"
        col_ent = msa_column_entropy(char_mat)
        H_col_mean = mean(col_ent)
        K_eff = effective_num_sequences(char_mat; threshold=0.8)
        # covariance spectrum
        X_onehot = onehot_encode(char_mat)
        C_mat = cov(X_onehot'; dims=1)  # covariance of one-hot vectors
        eigvals_C = sort(eigvals(C_mat), rev=true)
        λ1_ratio = eigvals_C[1] / sum(eigvals_C)
        @info "  H̄_col = $(round(H_col_mean, digits=3)), K_eff = $(round(K_eff, digits=1)), λ₁/tr(C) = $(round(λ1_ratio, digits=4))"

        # -- Step 5: Phase transition --
        @info "Step 5: Phase transition analysis …"
        pt = find_entropy_inflection(X̂; α=α_step)
        β_ret = round(Int, max(20 * pt.β_star, 50))   # retrieval regime
        β_gen = round(Int, max(2 * pt.β_star, 5))     # generation regime
        @info "  β* = $(round(pt.β_star, digits=2)), β_ret = $β_ret, β_gen = $β_gen"

        # save entropy curve
        p_ent = plot(pt.βs, pt.Hs ./ log(K),
            xlabel="β", ylabel="H(β) / log K",
            title="$(fam.name) (d=$d, K=$K)",
            xscale=:log10, label="H(β)/log K", lw=2, color=:coral,
            legend=:topright, size=(500, 350))
        vline!([pt.β_star], label="β* = $(round(pt.β_star, digits=1))", ls=:dash, color=:blue)
        savefig(p_ent, joinpath(PFAM_FIG_DIR, "entropy_$(fam.id).pdf"))

        # save family summary
        push!(family_summary, (
            Family=fam.name, PfamID=fam.id, K=K, L=L, d_PCA=d, d_full=d_full,
            H_col_mean=H_col_mean, K_eff=K_eff, λ1_ratio=λ1_ratio,
            β_star=pt.β_star, β_star_theory=pt.β_star_theory,
            β_retrieval=β_ret, β_generation=β_gen,
        ))

        # -- Step 6: Run SA --
        @info "Step 6: Running SA (retrieval β=$β_ret) …"
        sa_ret_samps = run_sa_chains(X̂, Float64(β_ret))
        @info "  SA retrieval: $(length(sa_ret_samps)) samples"

        @info "Step 6b: Running SA (generation β=$β_gen) …"
        sa_gen_samps = run_sa_chains(X̂, Float64(β_gen))
        @info "  SA generation: $(length(sa_gen_samps)) samples"

        # -- Step 7: Run MALA --
        @info "Step 7: Running MALA (β=$β_ret) …"
        mala_samps, mala_ar = run_mala_chains(X̂, Float64(β_ret))
        @info "  MALA: $(length(mala_samps)) samples, accept rate = $(round(mala_ar, digits=4))"

        # -- Step 8: Simple baselines --
        @info "Step 8: Running simple baselines …"
        baselines = run_simple_baselines(X̂, Float64(β_ret), S)

        # -- Step 9: Decode and evaluate --
        @info "Step 9: Decoding and evaluating …"

        # collect all methods
        methods = Dict{String, Vector{Vector{Float64}}}(
            "SA (retrieval)"  => sa_ret_samps,
            "SA (generation)" => sa_gen_samps,
            "MALA"            => mala_samps,
        )
        merge!(methods, baselines)

        β_eval = Float64(β_ret)
        for (method_name, samps) in methods
            seqs = [decode_sample(ξ, pca_model, L) for ξ in samps]
            metrics = evaluate_method(samps, seqs, X̂, β_eval, stored_seqs, L)
            push!(all_rows, (
                Family   = fam.name,
                PfamID   = fam.id,
                K        = K,
                L        = L,
                d_PCA    = d,
                Method   = method_name,
                β        = method_name == "SA (generation)" ? β_gen : β_ret,
                MALA_AR  = method_name == "MALA" ? mala_ar : missing,
                metrics...
            ))
        end

        # -- Step 10: Save generated sequences --
        @info "Step 10: Saving sequences …"
        for (label, samps) in [("sa_retrieval", sa_ret_samps), ("sa_generation", sa_gen_samps)]
            seqs = [decode_sample(ξ, pca_model, L) for ξ in samps]
            open(joinpath(fam_data_dir, "$(label)_sequences.fasta"), "w") do io
                for (i, seq) in enumerate(seqs)
                    println(io, ">$(label)_$(i)")
                    println(io, seq)
                end
            end
        end

        # save stored sequences
        open(joinpath(fam_data_dir, "stored_sequences.fasta"), "w") do io
            for (i, seq) in enumerate(stored_seqs)
                name = i <= length(seq_names) ? seq_names[i] : "stored_$i"
                println(io, ">$name")
                println(io, seq)
            end
        end

        # -- Step 11: Per-family AA frequency figure --
        @info "Step 11: Generating figures …"
        sa_ret_seqs = [decode_sample(ξ, pca_model, L) for ξ in sa_ret_samps]
        freq_stored = aa_freq_matrix(stored_seqs, L)
        freq_sa     = aa_freq_matrix(sa_ret_seqs, L)

        p1 = heatmap(freq_stored, xlabel="Position", ylabel="AA",
                     yticks=(1:N_AA, string.(AA_ALPHABET)),
                     title="Stored (K=$K)", color=:viridis, clims=(0,1), size=(700,250))
        p2 = heatmap(freq_sa, xlabel="Position", ylabel="AA",
                     yticks=(1:N_AA, string.(AA_ALPHABET)),
                     title="SA retrieval", color=:viridis, clims=(0,1), size=(700,250))
        p_freq = plot(p1, p2, layout=(2,1), size=(700,500),
                      plot_title="$(fam.name) — AA Frequencies")
        savefig(p_freq, joinpath(PFAM_FIG_DIR, "aa_freq_$(fam.id).pdf"))

        @info "  Done with $(fam.name)."
        println()
    end

    # ══════════════════════════════════════════════════════════════════════════
    # Aggregate results
    # ══════════════════════════════════════════════════════════════════════════

    df_results = DataFrame(all_rows)
    df_families = DataFrame(family_summary)

    # save
    CSV.write(joinpath(PFAM_DATA_DIR, "pfam_results.csv"), df_results)
    CSV.write(joinpath(PFAM_DATA_DIR, "pfam_family_summary.csv"), df_families)

    # print summary tables
    println("\n══════════════════════════════════════════════════════")
    println("FAMILY SUMMARY")
    println("══════════════════════════════════════════════════════")
    pretty_table(df_families[:, [:Family, :K, :L, :d_PCA, :H_col_mean, :K_eff, :β_star]];
                 formatters=PrettyTables.ft_printf("%.2f", [5,6,7]))

    println("\n══════════════════════════════════════════════════════")
    println("RESULTS BY FAMILY × METHOD")
    println("══════════════════════════════════════════════════════")
    pretty_table(df_results[:, [:Family, :Method, :Novelty, :Diversity, :SeqID, :KL_AA, :ValidAA]];
                 formatters=PrettyTables.ft_printf("%.4f", 3:7))

    # -- Cross-family comparison figure --
    @info "Generating cross-family comparison figure …"
    sa_ret_df = filter(r -> r.Method == "SA (retrieval)", df_results)
    sa_gen_df = filter(r -> r.Method == "SA (generation)", df_results)

    p_kl = groupedbar(sa_ret_df.Family,
        hcat(sa_ret_df.KL_AA, sa_gen_df.KL_AA),
        label=["SA retrieval" "SA generation"],
        ylabel="KL divergence (AA composition)",
        title="AA Composition Fidelity Across Families",
        color=[:coral :steelblue], alpha=0.8, size=(700, 400))
    savefig(p_kl, joinpath(PFAM_FIG_DIR, "cross_family_kl.pdf"))

    @info "All done. Results saved to $(PFAM_DATA_DIR)"
    return df_results, df_families
end

# ── Run ──────────────────────────────────────────────────────────────────────
df_results, df_families = run_pfam_experiment()
