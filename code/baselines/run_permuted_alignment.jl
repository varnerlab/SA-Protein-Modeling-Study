#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Column-Permuted Alignment Control
#
# For each family:
#   1. Create a permuted alignment by independently shuffling each column
#      (preserves per-position AA frequencies, destroys inter-position correlations)
#   2. Run the full SA pipeline on the permuted alignment
#   3. Compare KL, novelty, SeqID, and MI preservation between SA-on-real vs
#      SA-on-permuted
#
# This tests whether SA captures more than marginal per-position frequencies.
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment ..."
const _EXPERIMENT_ROOT = @__DIR__
include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))

const BASELINE_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "data")
const PFAM_DATA_DIR     = joinpath(_EXPERIMENT_ROOT, "..", "pfam-families", "data")
mkpath(BASELINE_DATA_DIR)

const FAMILIES = [
    (id="PF00076", name="RRM"),
    (id="PF00018", name="SH3"),
    (id="PF00397", name="WW"),
    (id="PF00014", name="Kunitz"),
    (id="PF00096", name="zf-C2H2"),
    (id="PF00595", name="PDZ"),
    (id="PF00069", name="Pkinase"),
    (id="PF00711", name="Defensin_beta"),
]

# SA parameters (identical to main experiment)
const alpha_step      = 0.01
const n_chains        = 30
const T_per_chain     = 5000
const T_burnin        = 2000
const thin_interval   = 100
const samples_per_chain = 5
const S_total         = n_chains * samples_per_chain  # 150
const sigma_init      = 0.01
const R_pca           = 0.95
const n_bootstrap     = 1000

# ══════════════════════════════════════════════════════════════════════════════
# Multi-chain SA sampler (same as main experiment)
# ══════════════════════════════════════════════════════════════════════════════

function run_sa_chains(X_hat::Matrix{Float64}, beta::Float64;
                       n_ch::Int=n_chains, T::Int=T_per_chain,
                       T_burn::Int=T_burnin, thin::Int=thin_interval,
                       spc::Int=samples_per_chain, sig0::Float64=sigma_init,
                       alpha::Float64=alpha_step, base_seed::Int=12345)
    d, K = size(X_hat)
    all_samples = Vector{Vector{Float64}}()

    pattern_indices = StatsBase.sample(1:K, n_ch, replace=(n_ch > K))
    for (c, k) in enumerate(pattern_indices)
        seed_c = base_seed + c
        xi0 = X_hat[:, k] .+ sig0 .* randn(d)
        (_, Xi) = sample(X_hat, xi0, T; β=beta, α=alpha, seed=seed_c)

        chain_pool = Vector{Vector{Float64}}()
        for ti in (T_burn+1):thin:T
            push!(chain_pool, Xi[ti, :])
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
# Pairwise MI computation
# ══════════════════════════════════════════════════════════════════════════════

"""
    compute_mi_vector(seqs, L) -> Vector{Float64}

Compute the upper-triangle pairwise MI vector for an alignment.
"""
function compute_mi_vector(seqs::Vector{String}, L::Int)
    n_pairs = L * (L - 1) / 2
    mi_vec = Float64[]

    for i in 1:(L-1)
        for j in (i+1):L
            # joint and marginal counts
            joint = zeros(N_AA, N_AA)
            marg_i = zeros(N_AA)
            marg_j = zeros(N_AA)
            total = 0

            for seq in seqs
                length(seq) < max(i, j) && continue
                ai = get(AA_TO_IDX, seq[i], 0)
                aj = get(AA_TO_IDX, seq[j], 0)
                (ai == 0 || aj == 0) && continue
                joint[ai, aj] += 1
                marg_i[ai] += 1
                marg_j[aj] += 1
                total += 1
            end

            total < 2 && (push!(mi_vec, 0.0); continue)

            # MI = sum p(x,y) log(p(x,y) / (p(x)p(y)))
            mi = 0.0
            for a in 1:N_AA, b in 1:N_AA
                pxy = joint[a, b] / total
                px = marg_i[a] / total
                py = marg_j[b] / total
                if pxy > 0 && px > 0 && py > 0
                    mi += pxy * log(pxy / (px * py))
                end
            end
            push!(mi_vec, mi)
        end
    end
    return mi_vec
end

"""
    bootstrap_kl_se(gen_seqs, stored_seqs; n_boot) -> (kl, se)
"""
function bootstrap_kl_se(gen_seqs, stored_seqs; n_boot=n_bootstrap)
    kl_point = aa_composition_kl(gen_seqs, stored_seqs)
    n = length(gen_seqs)
    boot_kls = [aa_composition_kl(gen_seqs[rand(1:n, n)], stored_seqs) for _ in 1:n_boot]
    return (kl_point, std(boot_kls))
end

"""
    project_to_pca(seq, pca_model, L) -> Vector{Float64}
"""
function project_to_pca(seq::String, pca_model, L::Int)
    oh = zeros(Float64, N_AA * L)
    for (pos, ch) in enumerate(seq)
        pos > L && break
        idx = get(AA_TO_IDX, ch, 0)
        idx > 0 && (oh[(pos-1)*N_AA + idx] = 1.0)
    end
    v = MultivariateStats.transform(pca_model, oh)
    nv = norm(v)
    nv > 1e-12 ? v ./ nv : v
end

# ══════════════════════════════════════════════════════════════════════════════
# Main experiment
# ══════════════════════════════════════════════════════════════════════════════

@info "Running column-permuted alignment control ..."

rows = Vector{NamedTuple}()

for fam in FAMILIES
    @info "  $(fam.name) ($(fam.id)) ..."
    fam_dir = joinpath(PFAM_DATA_DIR, fam.id)

    # load stored sequences
    stored_data = parse_fasta(joinpath(fam_dir, "stored_sequences.fasta"))
    stored_seqs = [s[2] for s in stored_data]
    char_mat = permutedims(hcat([collect(s) for s in stored_seqs]...))
    K, L = size(char_mat)

    # ── Real alignment SA (load existing results) ──
    X_hat_real, pca_model_real, L_real, d_full_real = build_memory_matrix(char_mat; pratio=R_pca)

    sa_data = parse_fasta(joinpath(fam_dir, "sa_generation_sequences.fasta"))
    sa_seqs_real = [s[2] for s in sa_data]
    sa_pca_real  = [project_to_pca(s, pca_model_real, L) for s in sa_seqs_real]

    kl_real, kl_real_se = bootstrap_kl_se(sa_seqs_real, stored_seqs)

    real_novelties = Float64[]
    real_seqids    = Float64[]
    for c in 1:n_chains
        si = (c - 1) * samples_per_chain + 1
        ei = c * samples_per_chain
        ei > length(sa_seqs_real) && break
        push!(real_novelties, mean(1.0 - nearest_cosine_similarity(sa_pca_real[j], X_hat_real) for j in si:ei))
        push!(real_seqids, mean(nearest_sequence_identity(sa_seqs_real[j], stored_seqs) for j in si:ei))
    end

    # MI for real SA
    mi_stored = compute_mi_vector(stored_seqs, L)
    mi_real   = compute_mi_vector(sa_seqs_real, L)
    mi_corr_real = length(mi_stored) > 1 ? cor(mi_stored, mi_real) : NaN

    # ── Permuted alignment: shuffle each column independently ──
    Random.seed!(54321)
    char_mat_perm = copy(char_mat)
    for j in 1:L
        char_mat_perm[:, j] = char_mat[randperm(K), j]
    end
    stored_seqs_perm = [String(char_mat_perm[k, :]) for k in 1:K]

    # build memory matrix on permuted alignment
    X_hat_perm, pca_model_perm, L_perm, d_full_perm = build_memory_matrix(char_mat_perm; pratio=R_pca)
    d_perm = size(X_hat_perm, 1)

    # find beta-star for permuted alignment
    pt_perm = find_entropy_inflection(X_hat_perm; α=alpha_step)
    beta_gen_perm = Float64(max(round(Int, 2 * pt_perm.β_star), 5))
    @info "    Permuted: d=$(d_perm), beta*=$(round(pt_perm.β_star, digits=2)), beta_gen=$(beta_gen_perm)"

    # run SA on permuted alignment
    Random.seed!(12345)
    sa_samps_perm = run_sa_chains(X_hat_perm, beta_gen_perm; base_seed=12345)
    sa_seqs_perm  = [decode_sample(xi, pca_model_perm, L_perm) for xi in sa_samps_perm]

    # evaluate against original stored sequences (not permuted)
    kl_perm, kl_perm_se = bootstrap_kl_se(sa_seqs_perm, stored_seqs)

    # project permuted SA sequences into the *original* PCA space for novelty comparison
    sa_pca_perm = [project_to_pca(s, pca_model_real, L) for s in sa_seqs_perm]

    perm_novelties = Float64[]
    perm_seqids    = Float64[]
    for c in 1:n_chains
        si = (c - 1) * samples_per_chain + 1
        ei = c * samples_per_chain
        ei > length(sa_seqs_perm) && break
        push!(perm_novelties, mean(1.0 - nearest_cosine_similarity(sa_pca_perm[j], X_hat_real) for j in si:ei))
        push!(perm_seqids, mean(nearest_sequence_identity(sa_seqs_perm[j], stored_seqs) for j in si:ei))
    end

    # MI for permuted SA
    mi_perm = compute_mi_vector(sa_seqs_perm, L)
    mi_corr_perm = length(mi_stored) > 1 ? cor(mi_stored, mi_perm) : NaN

    n_ch = min(length(real_novelties), length(perm_novelties))

    push!(rows, (
        Family       = fam.name,
        Method       = "SA (real)",
        KL_val       = kl_real,
        KL_se        = kl_real_se,
        Nov_val      = mean(real_novelties),
        Nov_se       = std(real_novelties) / sqrt(n_ch),
        SID_val      = mean(real_seqids),
        SID_se       = std(real_seqids) / sqrt(n_ch),
        MI_corr      = mi_corr_real,
    ))

    push!(rows, (
        Family       = fam.name,
        Method       = "SA (permuted)",
        KL_val       = kl_perm,
        KL_se        = kl_perm_se,
        Nov_val      = mean(perm_novelties),
        Nov_se       = std(perm_novelties) / sqrt(n_ch),
        SID_val      = mean(perm_seqids),
        SID_se       = std(perm_seqids) / sqrt(n_ch),
        MI_corr      = mi_corr_perm,
    ))

    @info "    Real:     KL=$(round(kl_real, digits=4)), MI_r=$(round(mi_corr_real, digits=3))"
    @info "    Permuted: KL=$(round(kl_perm, digits=4)), MI_r=$(round(mi_corr_perm, digits=3))"
end

# ══════════════════════════════════════════════════════════════════════════════
# Save results
# ══════════════════════════════════════════════════════════════════════════════

df = DataFrame(rows)
CSV.write(joinpath(BASELINE_DATA_DIR, "permuted_alignment_results.csv"), df)

println("\n======================================================")
println("PERMUTED ALIGNMENT CONTROL RESULTS")
println("======================================================")
println()

fmt(v, se) = "$(round(v, digits=3)) +/- $(round(se, digits=3))"

for fam in FAMILIES
    sub = filter(r -> r.Family == fam.name, df)
    println("$(fam.name):")
    for r in eachrow(sub)
        println("  $(rpad(r.Method, 16)) KL=$(fmt(r.KL_val, r.KL_se))  Nov=$(fmt(r.Nov_val, r.Nov_se))  " *
                "SID=$(fmt(r.SID_val, r.SID_se))  MI_r=$(round(r.MI_corr, digits=3))")
    end
end

@info "Done. Results saved to permuted_alignment_results.csv"
