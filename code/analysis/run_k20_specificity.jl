#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# K=20 Inter-Family Specificity Analysis
#
# Shows that SA-generated sequences at K=20 are family-specific, not collapsed
# to a generic mode. Uses the multifamily scaling study results.
#
# Method:
#   1. For each of the 4 families at K=20, subsample and generate sequences
#   2. Compute inter-family cosine similarity in a shared embedding space
#   3. Compare per-family amino acid frequency profiles
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment ..."
const _EXPERIMENT_ROOT = @__DIR__
include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))

const ANALYSIS_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "data")
const PFAM_DATA_DIR     = joinpath(_EXPERIMENT_ROOT, "..", "pfam-families", "data")
const SCALE_DATA_DIR    = joinpath(_EXPERIMENT_ROOT, "..", "scaling-study", "data")
mkpath(ANALYSIS_DATA_DIR)

# SA parameters (identical to multifamily scaling)
const alpha_step        = 0.01
const n_chains          = 30
const T_per_chain       = 5000
const T_burnin          = 2000
const thin_interval     = 100
const samples_per_chain = 5
const S_total           = n_chains * samples_per_chain
const sigma_init        = 0.01
const R_pca             = 0.95
const K_TARGET          = 20

const FAMILIES = [
    (id="PF00018", name="SH3",     K_full=55),
    (id="PF00014", name="Kunitz",  K_full=99),
    (id="PF00096", name="zf-C2H2", K_full=151),
]

# ══════════════════════════════════════════════════════════════════════════════
# Multi-chain SA sampler
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
# Generate K=20 sequences for each family
# ══════════════════════════════════════════════════════════════════════════════

@info "Generating K=20 sequences for each family ..."

family_seqs = Dict{String, Vector{String}}()
family_freq = Dict{String, Vector{Float64}}()  # global AA frequency profile

for fam in FAMILIES
    @info "  $(fam.name) ..."
    fam_dir = joinpath(PFAM_DATA_DIR, fam.id)

    # load and subsample
    sto_file = joinpath(fam_dir, "$(fam.id)_seed.sto")
    raw_seqs = parse_stockholm(sto_file)
    isempty(raw_seqs) && (raw_seqs = parse_fasta(sto_file))
    char_mat_full, _ = clean_alignment(raw_seqs)
    K_full = size(char_mat_full, 1)

    Random.seed!(20001)
    keep = StatsBase.sample(1:K_full, K_TARGET, replace=false) |> sort
    char_mat = char_mat_full[keep, :]
    stored_seqs = [String(char_mat[k, :]) for k in 1:K_TARGET]

    # build memory matrix and run SA
    X_hat, pca_model, L, d_full = build_memory_matrix(char_mat; pratio=R_pca)
    d = size(X_hat, 1)
    pt = find_entropy_inflection(X_hat; α=alpha_step)
    beta_gen = Float64(max(round(Int, 2 * pt.β_star), 5))

    Random.seed!(20001)
    sa_samps = run_sa_chains(X_hat, beta_gen; base_seed=20001)
    sa_seqs = [decode_sample(xi, pca_model, L) for xi in sa_samps]

    family_seqs[fam.name] = sa_seqs

    # global AA frequency profile
    counts = zeros(N_AA)
    for seq in sa_seqs, c in seq
        idx = get(AA_TO_IDX, c, 0)
        idx > 0 && (counts[idx] += 1)
    end
    total = sum(counts)
    family_freq[fam.name] = total > 0 ? counts ./ total : counts

    @info "    Generated $(length(sa_seqs)) sequences (d=$d, beta_gen=$beta_gen)"
end

# Also include WW at K=20 from the scaling study
@info "  WW (from scaling study) ..."
ww_dir = joinpath(SCALE_DATA_DIR, "PF00397")
sto_file_ww = joinpath(ww_dir, "PF00397_seed.sto")
raw_seqs_ww = parse_stockholm(sto_file_ww)
isempty(raw_seqs_ww) && (raw_seqs_ww = parse_fasta(sto_file_ww))
char_mat_ww_full, _ = clean_alignment(raw_seqs_ww)

Random.seed!(20001)
keep_ww = StatsBase.sample(1:size(char_mat_ww_full, 1), K_TARGET, replace=false) |> sort
char_mat_ww = char_mat_ww_full[keep_ww, :]
stored_seqs_ww = [String(char_mat_ww[k, :]) for k in 1:K_TARGET]

X_hat_ww, pca_model_ww, L_ww, _ = build_memory_matrix(char_mat_ww; pratio=R_pca)
pt_ww = find_entropy_inflection(X_hat_ww; α=alpha_step)
beta_gen_ww = Float64(max(round(Int, 2 * pt_ww.β_star), 5))
Random.seed!(20001)
sa_samps_ww = run_sa_chains(X_hat_ww, beta_gen_ww; base_seed=20001)
sa_seqs_ww = [decode_sample(xi, pca_model_ww, L_ww) for xi in sa_samps_ww]
family_seqs["WW"] = sa_seqs_ww

counts_ww = zeros(N_AA)
for seq in sa_seqs_ww, c in seq
    idx = get(AA_TO_IDX, c, 0)
    idx > 0 && (counts_ww[idx] += 1)
end
total_ww = sum(counts_ww)
family_freq["WW"] = total_ww > 0 ? counts_ww ./ total_ww : counts_ww

# ══════════════════════════════════════════════════════════════════════════════
# Inter-family similarity analysis
# ══════════════════════════════════════════════════════════════════════════════

@info "Computing inter-family similarity ..."

all_fam_names = ["SH3", "Kunitz", "zf-C2H2", "WW"]
n_fam = length(all_fam_names)

# 1. AA frequency cosine similarity between families
freq_sim_matrix = zeros(n_fam, n_fam)
for i in 1:n_fam
    for j in 1:n_fam
        fi = family_freq[all_fam_names[i]]
        fj = family_freq[all_fam_names[j]]
        freq_sim_matrix[i, j] = dot(fi, fj) / (norm(fi) * norm(fj))
    end
end

# 2. Jensen-Shannon divergence between AA frequency profiles
function js_divergence(p::Vector{Float64}, q::Vector{Float64})
    eps = 1e-10
    p_safe = p .+ eps; p_safe ./= sum(p_safe)
    q_safe = q .+ eps; q_safe ./= sum(q_safe)
    m = 0.5 .* (p_safe .+ q_safe)
    kl_pm = sum(p_safe[k] * log(p_safe[k] / m[k]) for k in eachindex(p_safe))
    kl_qm = sum(q_safe[k] * log(q_safe[k] / m[k]) for k in eachindex(q_safe))
    return 0.5 * (kl_pm + kl_qm)
end

jsd_matrix = zeros(n_fam, n_fam)
for i in 1:n_fam
    for j in 1:n_fam
        jsd_matrix[i, j] = js_divergence(family_freq[all_fam_names[i]], family_freq[all_fam_names[j]])
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# Save results
# ══════════════════════════════════════════════════════════════════════════════

# frequency similarity matrix
sim_rows = Vector{NamedTuple}()
for i in 1:n_fam
    for j in 1:n_fam
        push!(sim_rows, (
            Family1 = all_fam_names[i],
            Family2 = all_fam_names[j],
            FreqCosSim = freq_sim_matrix[i, j],
            JSD = jsd_matrix[i, j],
        ))
    end
end
sim_df = DataFrame(sim_rows)
CSV.write(joinpath(ANALYSIS_DATA_DIR, "k20_interfamily_similarity.csv"), sim_df)

# per-family AA frequency profiles
freq_rows = Vector{NamedTuple}()
for fam_name in all_fam_names
    for (aa_idx, aa) in enumerate(AA_ALPHABET)
        push!(freq_rows, (
            Family = fam_name,
            AA     = string(aa),
            Freq   = family_freq[fam_name][aa_idx],
        ))
    end
end
freq_df = DataFrame(freq_rows)
CSV.write(joinpath(ANALYSIS_DATA_DIR, "k20_aa_freq_profiles.csv"), freq_df)

# ══════════════════════════════════════════════════════════════════════════════
# Display results
# ══════════════════════════════════════════════════════════════════════════════

println("\n======================================================")
println("K=20 INTER-FAMILY SPECIFICITY")
println("======================================================")

println("\nAA Frequency Cosine Similarity Matrix:")
print("            ")
for name in all_fam_names
    print(rpad(name, 10))
end
println()
for i in 1:n_fam
    print(rpad(all_fam_names[i], 12))
    for j in 1:n_fam
        print(rpad(round(freq_sim_matrix[i, j], digits=3), 10))
    end
    println()
end

println("\nJensen-Shannon Divergence Matrix:")
print("            ")
for name in all_fam_names
    print(rpad(name, 10))
end
println()
for i in 1:n_fam
    print(rpad(all_fam_names[i], 12))
    for j in 1:n_fam
        print(rpad(round(jsd_matrix[i, j], digits=4), 10))
    end
    println()
end

# off-diagonal statistics
off_diag_sim = [freq_sim_matrix[i, j] for i in 1:n_fam for j in 1:n_fam if i != j]
off_diag_jsd = [jsd_matrix[i, j] for i in 1:n_fam for j in 1:n_fam if i != j]
println("\nOff-diagonal cosine similarity: mean=$(round(mean(off_diag_sim), digits=3)), " *
        "range=[$(round(minimum(off_diag_sim), digits=3)), $(round(maximum(off_diag_sim), digits=3))]")
println("Off-diagonal JSD: mean=$(round(mean(off_diag_jsd), digits=4)), " *
        "range=[$(round(minimum(off_diag_jsd), digits=4)), $(round(maximum(off_diag_jsd), digits=4))]")

@info "Done. Results saved to k20_*.csv"
