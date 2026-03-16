#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Consensus-with-Noise Baseline
#
# For each family:
#   1. Compute position-wise consensus (argmax per column) from stored alignment
#   2. One-hot encode the consensus into 20L space
#   3. Project to PCA space using the same PCA model as SA
#   4. Generate 150 sequences by adding calibrated Gaussian noise (sigma matched
#      to the empirical spread of stored patterns in PCA space) and decoding
#   5. Compute KL, novelty, SeqID using the same functions as SA
#
# This tests whether SA produces meaningfully different sequences than simply
# perturbing the family consensus.
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

const S_total    = 150
const n_chains   = 30
const spc        = 5   # samples_per_chain
const R_pca      = 0.95
const n_bootstrap = 1000

# ══════════════════════════════════════════════════════════════════════════════
# Consensus computation
# ══════════════════════════════════════════════════════════════════════════════

"""
    compute_consensus(char_mat) -> String

Compute the position-wise consensus sequence (most frequent AA per column).
"""
function compute_consensus(char_mat::Matrix{Char})
    K, L = size(char_mat)
    consensus = Char[]
    for j in 1:L
        counts = zeros(Int, N_AA)
        for i in 1:K
            idx = get(AA_TO_IDX, char_mat[i, j], 0)
            idx > 0 && (counts[idx] += 1)
        end
        best = argmax(counts)
        push!(consensus, AA_ALPHABET[best])
    end
    return String(consensus)
end

"""
    project_to_pca(seq, pca_model, L) -> Vector{Float64}

One-hot encode a single sequence, project to PCA space, and unit-normalize.
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

"""
    bootstrap_kl_se(gen_seqs, stored_seqs; n_boot) -> (kl, se)
"""
function bootstrap_kl_se(gen_seqs, stored_seqs; n_boot=n_bootstrap, rng=make_experiment_rng(314))
    kl_point = aa_composition_kl(gen_seqs, stored_seqs)
    n = length(gen_seqs)
    boot_kls = [aa_composition_kl(gen_seqs[rand(rng, 1:n, n)], stored_seqs) for _ in 1:n_boot]
    return (kl_point, std(boot_kls))
end

# ══════════════════════════════════════════════════════════════════════════════
# Main experiment
# ══════════════════════════════════════════════════════════════════════════════

@info "Running consensus-with-noise baseline ..."

# Validate that pfam-families results exist
let required = String[]
    for fam in FAMILIES
        push!(required, joinpath(PFAM_DATA_DIR, fam.id, "stored_sequences.fasta"))
        push!(required, joinpath(PFAM_DATA_DIR, fam.id, "sa_generation_sequences.fasta"))
    end
    assert_inputs_exist(required; context="run_consensus_baseline — run pfam-families experiment first")
end

rows = Vector{NamedTuple}()

for fam in FAMILIES
    @info "  $(fam.name) ($(fam.id)) ..."
    fam_dir = joinpath(PFAM_DATA_DIR, fam.id)

    # load stored sequences
    stored_data = parse_fasta(joinpath(fam_dir, "stored_sequences.fasta"))
    stored_seqs = [s[2] for s in stored_data]

    # build memory matrix (same as SA)
    char_mat = permutedims(hcat([collect(s) for s in stored_seqs]...))
    X_hat, pca_model, L, d_full = build_memory_matrix(char_mat; pratio=R_pca)
    d, K = size(X_hat)

    # compute consensus
    consensus_seq = compute_consensus(char_mat)
    @info "    Consensus: $(consensus_seq[1:min(30, length(consensus_seq))])..."

    # project consensus to PCA space
    consensus_pca = project_to_pca(consensus_seq, pca_model, L)

    # calibrate noise: match the empirical spread of stored patterns
    # sigma = mean distance of stored patterns from the centroid in PCA space
    centroid = vec(mean(X_hat, dims=2))
    dists = [norm(X_hat[:, k] .- centroid) for k in 1:K]
    sigma_empirical = mean(dists)
    @info "    Empirical spread sigma = $(round(sigma_empirical, digits=4))"

    # generate 150 sequences by adding calibrated noise to consensus
    rng = make_experiment_rng(42)
    consensus_noisy_pca = [consensus_pca .+ sigma_empirical .* randn(rng, d) for _ in 1:S_total]

    # decode to AA sequences
    consensus_noisy_seqs = [decode_sample(xi, pca_model, L) for xi in consensus_noisy_pca]

    # also load SA generation sequences for comparison
    sa_data = parse_fasta(joinpath(fam_dir, "sa_generation_sequences.fasta"))
    sa_seqs = [s[2] for s in sa_data]
    sa_pca  = [project_to_pca(s, pca_model, L) for s in sa_seqs]

    # ── Evaluate consensus-with-noise ──
    kl_cn, kl_cn_se = bootstrap_kl_se(consensus_noisy_seqs, stored_seqs)

    cn_novelties = Float64[]
    cn_seqids    = Float64[]
    for c in 1:n_chains
        si = (c - 1) * spc + 1
        ei = c * spc
        ei > S_total && break
        push!(cn_novelties, mean(1.0 - nearest_cosine_similarity(consensus_noisy_pca[j], X_hat) for j in si:ei))
        push!(cn_seqids, mean(nearest_sequence_identity(consensus_noisy_seqs[j], stored_seqs) for j in si:ei))
    end

    # ── Evaluate SA generation (for comparison) ──
    kl_sa, kl_sa_se = bootstrap_kl_se(sa_seqs, stored_seqs)

    sa_novelties = Float64[]
    sa_seqids    = Float64[]
    for c in 1:n_chains
        si = (c - 1) * spc + 1
        ei = c * spc
        ei > length(sa_seqs) && break
        push!(sa_novelties, mean(1.0 - nearest_cosine_similarity(sa_pca[j], X_hat) for j in si:ei))
        push!(sa_seqids, mean(nearest_sequence_identity(sa_seqs[j], stored_seqs) for j in si:ei))
    end

    n_ch = min(length(cn_novelties), length(sa_novelties))

    push!(rows, (
        Family       = fam.name,
        Method       = "Consensus+noise",
        KL_val       = kl_cn,
        KL_se        = kl_cn_se,
        Nov_val      = mean(cn_novelties),
        Nov_se       = std(cn_novelties) / sqrt(n_ch),
        SID_val      = mean(cn_seqids),
        SID_se       = std(cn_seqids) / sqrt(n_ch),
    ))

    push!(rows, (
        Family       = fam.name,
        Method       = "SA generation",
        KL_val       = kl_sa,
        KL_se        = kl_sa_se,
        Nov_val      = mean(sa_novelties),
        Nov_se       = std(sa_novelties) / sqrt(n_ch),
        SID_val      = mean(sa_seqids),
        SID_se       = std(sa_seqids) / sqrt(n_ch),
    ))

    @info "    CN:  KL=$(round(kl_cn, digits=4)), Nov=$(round(mean(cn_novelties), digits=3)), SID=$(round(mean(cn_seqids), digits=3))"
    @info "    SA:  KL=$(round(kl_sa, digits=4)), Nov=$(round(mean(sa_novelties), digits=3)), SID=$(round(mean(sa_seqids), digits=3))"
end

# ══════════════════════════════════════════════════════════════════════════════
# Save and display results
# ══════════════════════════════════════════════════════════════════════════════

df = DataFrame(rows)
CSV.write(joinpath(BASELINE_DATA_DIR, "consensus_baseline_results.csv"), df)

println("\n======================================================")
println("CONSENSUS-WITH-NOISE BASELINE RESULTS")
println("======================================================")
println()

fmt(v, se) = "$(round(v, digits=3)) +/- $(round(se, digits=3))"

println("% LaTeX table rows")
println("\\begin{tabular}{@{}llccc@{}}")
println("\\toprule")
println("Family & Method & KL & Novelty & SeqID \\\\")
println("\\midrule")
for fam in FAMILIES
    sub = filter(r -> r.Family == fam.name, df)
    for r in eachrow(sub)
        println("$(r.Family) & $(r.Method) & \$$(fmt(r.KL_val, r.KL_se))\$ & " *
                "\$$(fmt(r.Nov_val, r.Nov_se))\$ & \$$(fmt(r.SID_val, r.SID_se))\$ \\\\")
    end
    println("\\midrule")
end
println("\\bottomrule")
println("\\end{tabular}")

@info "Done. Results saved to consensus_baseline_results.csv"
