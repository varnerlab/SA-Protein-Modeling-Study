#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Scaling Study (Paper Section 3c)
#
# Fix the WW domain family (PF00397, K=420 seeds), subsample at
# K = {20, 50, 100, 200, 400}, and track generation quality as K shrinks.
# Shows where SA holds up vs. where baselines degrade.
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment …"
const _EXPERIMENT_ROOT = @__DIR__
include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))

const SCALE_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "data")
const SCALE_FIG_DIR  = joinpath(_EXPERIMENT_ROOT, "figs")
@info "Environment loaded."

# ══════════════════════════════════════════════════════════════════════════════
# Configuration
# ══════════════════════════════════════════════════════════════════════════════

const PFAM_ID   = "PF00397"   # WW domain (largest family, ~420 seeds)
const K_VALUES  = [20, 50, 100, 200, 400]
const N_REPEATS = 5            # independent subsamples per K for error bars

# SA parameters
const α_step        = 0.01
const n_chains      = 30
const T_per_chain   = 5000
const T_burnin      = 2000
const thin_interval = 100
const samples_per_chain = 5
const S             = n_chains * samples_per_chain  # 150
const σ_init        = 0.01
const R_pca         = 0.95

# ══════════════════════════════════════════════════════════════════════════════
# Multi-chain SA sampler (same as pfam-families)
# ══════════════════════════════════════════════════════════════════════════════

function run_sa_chains(X̂::Matrix{Float64}, β::Float64;
                       n_ch::Int=n_chains, T::Int=T_per_chain,
                       T_burn::Int=T_burnin, thin::Int=thin_interval,
                       spc::Int=samples_per_chain, σ0::Float64=σ_init,
                       α::Float64=α_step, base_seed::Int=12345)
    d, K = size(X̂)
    all_samples = Vector{Vector{Float64}}()

    pattern_indices = StatsBase.sample(1:K, n_ch, replace=(n_ch > K))
    for (c, k) in enumerate(pattern_indices)
        seed_c = base_seed + c
        ξ₀ = X̂[:, k] .+ σ0 .* randn(d)
        (_, Ξ) = sample(X̂, ξ₀, T; β=β, α=α, seed=seed_c)

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
    return all_samples
end

# ══════════════════════════════════════════════════════════════════════════════
# Simple baselines
# ══════════════════════════════════════════════════════════════════════════════

function run_baselines(X̂::Matrix{Float64}, β::Float64, S::Int; seed::Int=12345)
    d, K = size(X̂)
    σ_noise = sqrt(2 * α_step / β)

    Random.seed!(seed)
    bs = [copy(X̂[:, rand(1:K)]) for _ in 1:S]

    Random.seed!(seed)
    gp = [X̂[:, rand(1:K)] .+ σ_noise .* randn(d) for _ in 1:S]

    dirichlet = Dirichlet(K, 1.0)
    Random.seed!(seed)
    rc = [X̂ * rand(dirichlet) for _ in 1:S]

    return Dict("Bootstrap" => bs, "Gaussian pert." => gp, "Convex comb." => rc)
end

# ══════════════════════════════════════════════════════════════════════════════
# Main experiment
# ══════════════════════════════════════════════════════════════════════════════

function run_scaling_study()

    mkpath(SCALE_DATA_DIR)
    mkpath(SCALE_FIG_DIR)

    # -- Step 1: Load full family --
    @info "Step 1: Loading full WW family …"
    fam_cache = joinpath(SCALE_DATA_DIR, PFAM_ID)
    mkpath(fam_cache)
    sto_file = download_pfam_seed(PFAM_ID; cache_dir=fam_cache)
    raw_seqs = parse_stockholm(sto_file)
    isempty(raw_seqs) && (raw_seqs = parse_fasta(sto_file))
    @info "  Parsed $(length(raw_seqs)) sequences"

    char_mat_full, seq_names_full = clean_alignment(raw_seqs)
    K_full, L = size(char_mat_full)
    stored_seqs_full = [String(char_mat_full[k, :]) for k in 1:K_full]
    @info "  Full family: K=$K_full, L=$L"

    # -- Step 2: Run at each K --
    all_rows = Vector{NamedTuple}()

    for K_sub in K_VALUES
        @info "══════════════════════════════════════════════════════"
        @info "  K = $K_sub"
        @info "══════════════════════════════════════════════════════"

        if K_sub > K_full
            @warn "  K_sub=$K_sub > K_full=$K_full, skipping"
            continue
        end

        for rep in 1:N_REPEATS
            @info "  Repeat $rep/$N_REPEATS …"
            rep_seed = 1000 * K_sub + rep

            # subsample
            Random.seed!(rep_seed)
            keep = StatsBase.sample(1:K_full, K_sub, replace=false) |> sort
            char_mat = char_mat_full[keep, :]
            stored_seqs = [String(char_mat[k, :]) for k in 1:K_sub]

            # build memory matrix
            X̂, pca_model, L_out, d_full = build_memory_matrix(char_mat; pratio=R_pca)
            d = size(X̂, 1)

            # phase transition → select β
            pt = find_entropy_inflection(X̂; α=α_step)
            β_ret = Float64(max(round(Int, 20 * pt.β_star), 50))
            β_gen = Float64(max(round(Int, 2 * pt.β_star), 5))

            # -- SA (generation + retrieval) --
            Random.seed!(rep_seed)
            sa_gen_samps = run_sa_chains(X̂, β_gen; base_seed=rep_seed)
            Random.seed!(rep_seed + 10000)
            sa_ret_samps = run_sa_chains(X̂, β_ret; base_seed=rep_seed + 10000)

            # -- Baselines --
            baselines = run_baselines(X̂, β_ret, S; seed=rep_seed + 20000)

            # -- Evaluate all methods --
            methods = Dict(
                "SA (generation)" => (sa_gen_samps, β_gen),
                "SA (retrieval)"  => (sa_ret_samps, β_ret),
            )
            for (name, samps) in baselines
                methods[name] = (samps, β_ret)
            end

            for (method_name, (samps, β_eval)) in methods
                seqs = [decode_sample(ξ, pca_model, L_out) for ξ in samps]

                novelty_vals = [sample_novelty(ξ, X̂) for ξ in samps]
                energy_vals  = [hopfield_energy(ξ, X̂, β_eval) for ξ in samps]
                seq_id_vals  = [nearest_sequence_identity(s, stored_seqs) for s in seqs]

                push!(all_rows, (
                    K        = K_sub,
                    Repeat   = rep,
                    d_PCA    = d,
                    Method   = method_name,
                    β        = β_eval,
                    β_star   = pt.β_star,
                    Novelty  = mean(novelty_vals),
                    Diversity = sample_diversity(samps),
                    Energy   = mean(energy_vals),
                    SeqID    = mean(seq_id_vals),
                    ValidAA  = mean(valid_residue_fraction(s) for s in seqs),
                    KL_AA    = aa_composition_kl(seqs, stored_seqs),
                ))
            end
        end
    end

    df = DataFrame(all_rows)
    CSV.write(joinpath(SCALE_DATA_DIR, "scaling_results.csv"), df)
    @info "Results saved to $(joinpath(SCALE_DATA_DIR, "scaling_results.csv"))"

    # ══════════════════════════════════════════════════════════════════════════
    # Figures
    # ══════════════════════════════════════════════════════════════════════════

    @info "Generating figures …"

    # aggregate: mean ± std over repeats
    gdf = groupby(df, [:K, :Method])
    agg = combine(gdf,
        :KL_AA    => mean => :KL_mean,     :KL_AA    => std => :KL_std,
        :Novelty  => mean => :Nov_mean,    :Novelty  => std => :Nov_std,
        :SeqID    => mean => :SID_mean,    :SeqID    => std => :SID_std,
        :Diversity => mean => :Div_mean,   :Diversity => std => :Div_std,
        :β_star   => mean => :β_star_mean,
    )

    method_colors = Dict(
        "SA (generation)" => :steelblue,
        "SA (retrieval)"  => :coral,
        "Bootstrap"       => :gray,
        "Gaussian pert."  => :orange,
        "Convex comb."    => :purple,
    )
    method_order = ["SA (generation)", "SA (retrieval)", "Bootstrap", "Gaussian pert.", "Convex comb."]

    # -- KL divergence vs K --
    p_kl = plot(xlabel="Family size K", ylabel="KL divergence (AA composition)",
                title="Generation Quality vs. Family Size (WW domain)",
                xscale=:log10, size=(700, 450), legend=:topright)
    for m in method_order
        sub = filter(r -> r.Method == m, agg)
        isempty(sub) && continue
        plot!(p_kl, sub.K, sub.KL_mean, ribbon=sub.KL_std,
              label=m, lw=2, color=method_colors[m], alpha=0.8,
              markershape=:circle, markersize=5)
    end
    savefig(p_kl, joinpath(SCALE_FIG_DIR, "scaling_kl_vs_k.pdf"))

    # -- Novelty vs K --
    p_nov = plot(xlabel="Family size K", ylabel="Novelty",
                 title="Sample Novelty vs. Family Size",
                 xscale=:log10, size=(700, 450), legend=:topright)
    for m in method_order
        sub = filter(r -> r.Method == m, agg)
        isempty(sub) && continue
        plot!(p_nov, sub.K, sub.Nov_mean, ribbon=sub.Nov_std,
              label=m, lw=2, color=method_colors[m], alpha=0.8,
              markershape=:circle, markersize=5)
    end
    savefig(p_nov, joinpath(SCALE_FIG_DIR, "scaling_novelty_vs_k.pdf"))

    # -- Sequence identity vs K --
    p_sid = plot(xlabel="Family size K", ylabel="Nearest sequence identity",
                 title="Sequence Identity vs. Family Size",
                 xscale=:log10, size=(700, 450), legend=:bottomright)
    for m in method_order
        sub = filter(r -> r.Method == m, agg)
        isempty(sub) && continue
        plot!(p_sid, sub.K, sub.SID_mean, ribbon=sub.SID_std,
              label=m, lw=2, color=method_colors[m], alpha=0.8,
              markershape=:circle, markersize=5)
    end
    savefig(p_sid, joinpath(SCALE_FIG_DIR, "scaling_seqid_vs_k.pdf"))

    # -- β* vs K --
    sa_gen_agg = filter(r -> r.Method == "SA (generation)", agg)
    p_beta = plot(sa_gen_agg.K, sa_gen_agg.β_star_mean,
                  xlabel="Family size K", ylabel="β*",
                  title="Critical Temperature vs. Family Size",
                  xscale=:log10, lw=2, color=:steelblue,
                  markershape=:circle, markersize=6, label="β* (empirical)",
                  size=(600, 400), legend=:topleft)
    hline!([sqrt(mean(sa_gen_agg.β_star_mean .^ 2))], ls=:dash, color=:gray, label="mean β*")
    savefig(p_beta, joinpath(SCALE_FIG_DIR, "scaling_beta_star_vs_k.pdf"))

    @info "All figures saved to $SCALE_FIG_DIR"

    # -- Print summary table --
    println("\n══════════════════════════════════════════════════════")
    println("SCALING STUDY SUMMARY (WW domain, PF00397)")
    println("══════════════════════════════════════════════════════")
    for m in method_order
        sub = filter(r -> r.Method == m, agg)
        isempty(sub) && continue
        println("\n  $m:")
        for row in eachrow(sub)
            println("    K=$(row.K):  KL=$(round(row.KL_mean, digits=4)) ± $(round(row.KL_std, digits=4)),  " *
                    "Novelty=$(round(row.Nov_mean, digits=3)),  SeqID=$(round(row.SID_mean, digits=3))")
        end
    end

    return df, agg
end

# ── Run ──────────────────────────────────────────────────────────────────────
df, agg = run_scaling_study()
