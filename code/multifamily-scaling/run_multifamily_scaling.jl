#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Multi-Family Scaling Study at K=20
#
# Subsample SH3 (K=55), Kunitz (K=99), and zf-C2H2 (K=151) to K=20
# with 5 independent replicates per family. Complements the WW-only scaling
# study to show that SA robustness at small K generalizes across families.
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment …"
const _EXPERIMENT_ROOT = @__DIR__
include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))

const MFSCALE_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "data")
const MFSCALE_FIG_DIR  = joinpath(_EXPERIMENT_ROOT, "figs")
const PFAM_DATA_DIR    = joinpath(_EXPERIMENT_ROOT, "..", "pfam-families", "data")

# ══════════════════════════════════════════════════════════════════════════════
# Configuration
# ══════════════════════════════════════════════════════════════════════════════

const FAMILIES = [
    (id="PF00018", name="SH3",     K_full=55),
    (id="PF00014", name="Kunitz",  K_full=99),
    (id="PF00096", name="zf-C2H2", K_full=151),
]

const K_TARGET  = 20
const N_REPEATS = 5

# SA parameters (identical to main experiment)
const α_mix        = 0.01
const n_chains      = 30
const T_per_chain   = 5000
const T_burnin      = 2000
const thin_interval = 100
const samples_per_chain = 5
const S             = n_chains * samples_per_chain  # 150
const σ_init        = 0.01
const R_pca         = 0.95

# ══════════════════════════════════════════════════════════════════════════════
# Multi-chain SA sampler
# ══════════════════════════════════════════════════════════════════════════════

function run_sa_chains(X̂::Matrix{Float64}, β::Float64;
                       n_ch::Int=n_chains, T::Int=T_per_chain,
                       T_burn::Int=T_burnin, thin::Int=thin_interval,
                       spc::Int=samples_per_chain, σ0::Float64=σ_init,
                       α::Float64=α_mix, base_seed::Int=12345)
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
# Main experiment
# ══════════════════════════════════════════════════════════════════════════════

function run_multifamily_scaling()
    mkpath(MFSCALE_DATA_DIR)
    mkpath(MFSCALE_FIG_DIR)

    all_rows = Vector{NamedTuple}()

    for fam in FAMILIES
        @info "══════════════════════════════════════════════════════"
        @info "  $(fam.name) ($(fam.id)), K_full=$(fam.K_full) → K=$(K_TARGET)"
        @info "══════════════════════════════════════════════════════"

        # Load full family alignment
        sto_file = joinpath(PFAM_DATA_DIR, fam.id, "$(fam.id)_seed.sto")
        if !isfile(sto_file)
            @warn "  Missing alignment: $sto_file — skipping"
            continue
        end
        raw_seqs = parse_stockholm(sto_file)
        isempty(raw_seqs) && (raw_seqs = parse_fasta(sto_file))
        char_mat_full, seq_names_full = clean_alignment(raw_seqs)
        K_full, L = size(char_mat_full)
        @info "  Full family: K=$K_full, L=$L"

        if K_TARGET >= K_full
            @warn "  K_TARGET=$K_TARGET >= K_full=$K_full, skipping"
            continue
        end

        for rep in 1:N_REPEATS
            @info "  Repeat $rep/$N_REPEATS …"
            rep_seed = 1000 * K_TARGET + 100 * findfirst(f -> f.id == fam.id, FAMILIES) + rep

            # subsample
            Random.seed!(rep_seed)
            keep = StatsBase.sample(1:K_full, K_TARGET, replace=false) |> sort
            char_mat = char_mat_full[keep, :]
            stored_seqs = [String(char_mat[k, :]) for k in 1:K_TARGET]

            # build memory matrix
            X̂, pca_model, L_out, d_full = build_memory_matrix(char_mat; pratio=R_pca)
            d = size(X̂, 1)

            # phase transition
            pt = find_entropy_inflection(X̂; α=α_mix)
            β_gen = Float64(max(round(Int, 2 * pt.β_star), 5))

            # SA generation
            Random.seed!(rep_seed)
            sa_samps = run_sa_chains(X̂, β_gen; base_seed=rep_seed)
            sa_seqs = [decode_sample(ξ, pca_model, L_out) for ξ in sa_samps]

            # evaluate
            novelty_vals = [sample_novelty(ξ, X̂) for ξ in sa_samps]
            seq_id_vals  = [nearest_sequence_identity(s, stored_seqs) for s in sa_seqs]

            push!(all_rows, (
                Family   = fam.name,
                PfamID   = fam.id,
                K_full   = K_full,
                K_sub    = K_TARGET,
                Repeat   = rep,
                d_PCA    = d,
                β_star   = pt.β_star,
                β_gen    = β_gen,
                Novelty  = mean(novelty_vals),
                Diversity = sample_diversity(sa_samps),
                SeqID    = mean(seq_id_vals),
                ValidAA  = mean(valid_residue_fraction(s) for s in sa_seqs),
                KL_AA    = aa_composition_kl(sa_seqs, stored_seqs),
            ))

            @info "    KL=$(round(all_rows[end].KL_AA, digits=4)), Nov=$(round(all_rows[end].Novelty, digits=3)), SeqID=$(round(all_rows[end].SeqID, digits=3))"
        end
    end

    df = DataFrame(all_rows)
    CSV.write(joinpath(MFSCALE_DATA_DIR, "multifamily_scaling_k20.csv"), df)
    @info "Results saved."

    # ══════════════════════════════════════════════════════════════════════════
    # Summary table
    # ══════════════════════════════════════════════════════════════════════════

    println("\n══════════════════════════════════════════════════════")
    println("MULTI-FAMILY SCALING STUDY AT K=20")
    println("══════════════════════════════════════════════════════")
    gdf = groupby(df, :Family)
    for (key, sub) in pairs(gdf)
        fam = key.Family
        kl_m = round(mean(sub.KL_AA), digits=4)
        kl_s = round(std(sub.KL_AA), digits=4)
        nov_m = round(mean(sub.Novelty), digits=3)
        nov_s = round(std(sub.Novelty), digits=3)
        sid_m = round(mean(sub.SeqID), digits=3)
        sid_s = round(std(sub.SeqID), digits=3)
        println("  $fam:  KL=$kl_m ± $kl_s,  Nov=$nov_m ± $nov_s,  SeqID=$sid_m ± $sid_s")
    end

    # Also include WW K=20 from original scaling study for comparison
    ww_csv = joinpath(_EXPERIMENT_ROOT, "..", "scaling-study", "data", "scaling_results.csv")
    if isfile(ww_csv)
        ww_df = CSV.read(ww_csv, DataFrame)
        ww_k20 = filter(r -> r.K == 20 && r.Method == "SA (generation)", ww_df)
        if nrow(ww_k20) > 0
            kl_m = round(mean(ww_k20.KL_AA), digits=4)
            kl_s = round(std(ww_k20.KL_AA), digits=4)
            nov_m = round(mean(ww_k20.Novelty), digits=3)
            nov_s = round(std(ww_k20.Novelty), digits=3)
            sid_m = round(mean(ww_k20.SeqID), digits=3)
            sid_s = round(std(ww_k20.SeqID), digits=3)
            println("  WW (from original):  KL=$kl_m ± $kl_s,  Nov=$nov_m ± $nov_s,  SeqID=$sid_m ± $sid_s")
        end
    end

    return df
end

# ── Run ──────────────────────────────────────────────────────────────────────
df = run_multifamily_scaling()
