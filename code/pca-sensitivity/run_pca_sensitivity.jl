#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# PCA variance threshold sensitivity study (Reviewer Q2)
#
# Runs SA generation at pratio ∈ {0.80, 0.85, 0.90, 0.95, 0.99} for all 8
# families, comparing KL divergence, novelty, nearest sequence identity,
# PCA dimension d, and critical temperature β*.
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment …"
const _EXPERIMENT_ROOT = @__DIR__

using Pkg
Pkg.activate(joinpath(_EXPERIMENT_ROOT, ".."))
include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))

using Printf

# ── SA chain runner (from run_pfam_families.jl) ───────────────────────────────

function run_sa_chains(X̂::Matrix{Float64}, β::Float64;
    n_chains::Int=30, T::Int=5000, T_burn::Int=2000,
    thin::Int=100, spc::Int=5, σ_init::Float64=0.01, α::Float64=0.01)

    d, K = size(X̂)
    all_samples = Vector{Vector{Float64}}()

    Random.seed!(42)
    pattern_indices = StatsBase.sample(1:K, n_chains, replace=(n_chains > K))

    for (c, k) in enumerate(pattern_indices)
        seed_c = 12345 + c
        Random.seed!(seed_c)
        ξ₀ = X̂[:, k] .+ σ_init .* randn(d)
        (_, Ξ) = sample(X̂, ξ₀, T; β=β, α=α, seed=seed_c)

        chain_pool = Vector{Vector{Float64}}()
        for tᵢ in (T_burn+1):thin:T
            push!(chain_pool, Ξ[tᵢ, :])
        end

        n_avail = length(chain_pool)
        idxs = round.(Int, range(1, n_avail, length=min(spc, n_avail)))
        for idx in idxs
            push!(all_samples, chain_pool[idx])
        end
    end
    return all_samples
end

# ── Family configuration ──────────────────────────────────────────────────────

families = [
    ("PF00076", "RRM"),
    ("PF00018", "SH3"),
    ("PF00397", "WW"),
    ("PF00014", "Kunitz"),
    ("PF00096", "zf-C2H2"),
    ("PF00595", "PDZ"),
    ("PF00069", "Pkinase"),
    ("PF00711", "Defensin_beta"),
]

pratios = [0.80, 0.85, 0.90, 0.95, 0.99]
α_step = 0.01

data_dir = joinpath(_EXPERIMENT_ROOT, "..", "pfam-families", "data")
output_dir = joinpath(_EXPERIMENT_ROOT, "results")
mkpath(output_dir)

# ── Main loop ─────────────────────────────────────────────────────────────────

results = []

for (pfam_id, family_name) in families
    @info "═══ $family_name ($pfam_id) ═══"

    # load alignment
    sto_file = joinpath(data_dir, pfam_id, "$(pfam_id)_seed.sto")
    raw_seqs = parse_stockholm(sto_file)
    char_mat, retained_names = clean_alignment(raw_seqs)
    K, L = size(char_mat)
    stored_seqs = [String(char_mat[i, :]) for i in 1:K]

    for pratio in pratios
        @info "  pratio=$pratio"

        # build memory matrix
        X̂, pca_model, _, d_full = build_memory_matrix(char_mat; pratio=pratio)
        d = size(X̂, 1)

        # find β*
        pt = find_entropy_inflection(X̂; α=α_step)
        β_star = pt.β_star
        β_gen = Float64(round(Int, max(2 * β_star, 5)))

        @info "    d=$d, β*=$(round(β_star, digits=2)), β_gen=$β_gen"

        # run SA generation
        samps = run_sa_chains(X̂, β_gen; α=α_step)
        seqs = [decode_sample(ξ, pca_model, L) for ξ in samps]

        # compute metrics
        kl = aa_composition_kl(seqs, stored_seqs)
        novelty_vals = [1.0 - maximum(X̂' * (s ./ (norm(s) + 1e-12))) for s in samps]
        mean_novelty = mean(novelty_vals)
        seqid_vals = [nearest_sequence_identity(s, stored_seqs) for s in seqs]
        mean_seqid = mean(seqid_vals)
        valid_aa = mean(valid_residue_fraction(s) for s in seqs)

        @info @sprintf("    KL=%.4f, Novelty=%.3f, SeqID=%.3f, ValidAA=%.3f",
            kl, mean_novelty, mean_seqid, valid_aa)

        push!(results, (
            family=family_name, pfam_id=pfam_id, pratio=pratio,
            K=K, L=L, d=d, d_full=d_full,
            beta_star=β_star, beta_gen=β_gen,
            KL=kl, novelty=mean_novelty, seqid=mean_seqid, valid_aa=valid_aa,
        ))
    end
    println()
end

# ── Write CSV ─────────────────────────────────────────────────────────────────

df = DataFrame(results)
csv_path = joinpath(output_dir, "pca_sensitivity_results.csv")
CSV.write(csv_path, df)
@info "Results saved to $csv_path"

# ── Print summary table ──────────────────────────────────────────────────────

println("\n" * "="^90)
println("PCA VARIANCE THRESHOLD SENSITIVITY")
println("="^90)
@printf("%-15s %6s %5s %6s %8s %8s %8s %8s\n",
    "Family", "pratio", "d", "β*", "β_gen", "KL", "Novelty", "SeqID")
println("-"^90)

for row in results
    @printf("%-15s %6.2f %5d %6.1f %8.0f %8.4f %8.3f %8.3f\n",
        row.family, row.pratio, row.d, row.beta_star, row.beta_gen,
        row.KL, row.novelty, row.seqid)
end

# ── Per-family summary: how much do metrics change? ───────────────────────────

println("\n\nSensitivity summary (range across pratio values):")
println("-"^70)
@printf("%-15s  %12s  %12s  %12s\n", "Family", "KL range", "Nov range", "SeqID range")
println("-"^70)

for (pfam_id, family_name) in families
    fam_rows = filter(r -> r.family == family_name, results)
    kls = [r.KL for r in fam_rows]
    novs = [r.novelty for r in fam_rows]
    sids = [r.seqid for r in fam_rows]
    @printf("%-15s  %5.3f-%-5.3f  %5.3f-%-5.3f  %5.3f-%-5.3f\n",
        family_name,
        minimum(kls), maximum(kls),
        minimum(novs), maximum(novs),
        minimum(sids), maximum(sids))
end

@info "\nDone."
