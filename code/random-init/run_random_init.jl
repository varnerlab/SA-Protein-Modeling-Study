#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Random sphere initialization experiment (Reviewer Q4)
#
# Compares SA generation with two initialization strategies:
#   (A) Standard: initialize near stored patterns (σ_init = 0.01)
#   (B) Random:   initialize from random points on the unit sphere
#
# Same chain parameters, same β_gen, same number of samples.
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment …"
const _EXPERIMENT_ROOT = @__DIR__

using Pkg
Pkg.activate(joinpath(_EXPERIMENT_ROOT, ".."))
include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))

using Printf, LinearAlgebra, MultivariateStats

"""
    project_to_pca(seq::String, pca_model, L::Int) -> Vector{Float64}

One-hot encode a single sequence, project into PCA space, and unit-normalize.
"""
function project_to_pca(seq::String, pca_model, L::Int)
    AA_ALPHABET = collect("ACDEFGHIKLMNPQRSTVWY")
    n_aa = length(AA_ALPHABET)
    aa_to_idx = Dict(aa => i for (i, aa) in enumerate(AA_ALPHABET))
    oh = zeros(Float64, n_aa * L)
    for (pos, ch) in enumerate(seq)
        pos > L && break
        idx = get(aa_to_idx, ch, 0)
        idx > 0 && (oh[(pos-1)*n_aa + idx] = 1.0)
    end
    v = MultivariateStats.transform(pca_model, oh)
    nv = norm(v)
    nv > 1e-12 ? v ./ nv : v
end

# ── SA chain runners ─────────────────────────────────────────────────────────

function run_sa_chains_stored(X̂::Matrix{Float64}, β::Float64;
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

function run_sa_chains_random(X̂::Matrix{Float64}, β::Float64;
    n_chains::Int=30, T::Int=5000, T_burn::Int=2000,
    thin::Int=100, spc::Int=5, α::Float64=0.01)

    d, K = size(X̂)
    all_samples = Vector{Vector{Float64}}()

    for c in 1:n_chains
        seed_c = 12345 + c
        Random.seed!(seed_c)
        # Random point on unit sphere
        z = randn(d)
        ξ₀ = z ./ norm(z)
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

    # build memory matrix (standard 95%)
    X̂, pca_model, _, d_full = build_memory_matrix(char_mat; pratio=0.95)
    d = size(X̂, 1)

    # find β*
    pt = find_entropy_inflection(X̂; α=α_step)
    β_star = pt.β_star
    β_gen = Float64(round(Int, max(2 * β_star, 5)))

    @info "  d=$d, K=$K, β*=$(round(β_star, digits=2)), β_gen=$β_gen"

    for (init_name, runner) in [("stored", run_sa_chains_stored), ("random", run_sa_chains_random)]
        @info "  init=$init_name"

        samps = runner(X̂, β_gen; α=α_step)
        seqs = [decode_sample(ξ, pca_model, L) for ξ in samps]

        # compute metrics
        kl = aa_composition_kl(seqs, stored_seqs)
        # novelty: re-encode decoded sequences to PCA space (round-trip consistent with tab:results)
        gen_pca = [project_to_pca(s, pca_model, L) for s in seqs]
        novelty_vals = [1.0 - nearest_cosine_similarity(g, X̂) for g in gen_pca]
        mean_novelty = mean(novelty_vals)
        seqid_vals = [nearest_sequence_identity(s, stored_seqs) for s in seqs]
        mean_seqid = mean(seqid_vals)
        valid_aa = mean(valid_residue_fraction(s) for s in seqs)

        # measure how close final samples are to stored patterns (cosine sim, also round-trip)
        max_cos_vals = [nearest_cosine_similarity(g, X̂) for g in gen_pca]
        mean_max_cos = mean(max_cos_vals)

        @info @sprintf("    KL=%.4f, Novelty=%.3f, SeqID=%.3f, ValidAA=%.3f, MaxCos=%.3f",
            kl, mean_novelty, mean_seqid, valid_aa, mean_max_cos)

        push!(results, (
            family=family_name, pfam_id=pfam_id,
            init=init_name, K=K, L=L, d=d,
            beta_star=β_star, beta_gen=β_gen,
            KL=kl, novelty=mean_novelty, seqid=mean_seqid,
            valid_aa=valid_aa, max_cos=mean_max_cos,
        ))
    end
    println()
end

# ── Write CSV ─────────────────────────────────────────────────────────────────

df = DataFrame(results)
csv_path = joinpath(output_dir, "random_init_results.csv")
CSV.write(csv_path, df)
@info "Results saved to $csv_path"

# ── Print comparison table ────────────────────────────────────────────────────

println("\n" * "="^95)
println("INITIALIZATION COMPARISON: STORED vs RANDOM SPHERE")
println("="^95)
@printf("%-15s %7s %8s %8s %8s %8s %8s\n",
    "Family", "Init", "KL", "Novelty", "SeqID", "ValidAA", "MaxCos")
println("-"^95)

for (pfam_id, family_name) in families
    fam_rows = filter(r -> r.family == family_name, results)
    for r in fam_rows
        @printf("%-15s %7s %8.4f %8.3f %8.3f %8.3f %8.3f\n",
            r.family, r.init, r.KL, r.novelty, r.seqid, r.valid_aa, r.max_cos)
    end
end

# ── Difference summary ────────────────────────────────────────────────────────

println("\n\nDifference (random - stored):")
println("-"^70)
@printf("%-15s %10s %10s %10s\n", "Family", "ΔKL", "ΔNovelty", "ΔSeqID")
println("-"^70)

for (pfam_id, family_name) in families
    fam_rows = filter(r -> r.family == family_name, results)
    stored = filter(r -> r.init == "stored", fam_rows)[1]
    random = filter(r -> r.init == "random", fam_rows)[1]
    @printf("%-15s %+10.4f %+10.3f %+10.3f\n",
        family_name,
        random.KL - stored.KL,
        random.novelty - stored.novelty,
        random.seqid - stored.seqid)
end

@info "\nDone."
