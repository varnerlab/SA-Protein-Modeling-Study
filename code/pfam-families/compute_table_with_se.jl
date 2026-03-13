#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Compute Table 2 values with proper ±SE for each metric
#
# KL:      bootstrap SE — resample 150 sequences 1000× with replacement,
#          compute KL each time, SE = std of bootstrap distribution
# Novelty: per-chain mean of 1 - max_k cos(ξ, m_k) in PCA space
# SeqID:   per-chain mean of nearest sequence identity in AA space
#
# To get PCA-space novelty we rebuild the memory matrix from stored sequences,
# then project generated sequences into PCA space.
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment …"
const _EXPERIMENT_ROOT = @__DIR__
include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))

const PFAM_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "data")

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

const n_chains = 30
const samples_per_chain = 5
const S_total = n_chains * samples_per_chain  # 150
const n_bootstrap = 1000
const R_pca = 0.95
const α_step = 0.01

# load family summary for β values
const fam_summary = CSV.read(joinpath(PFAM_DATA_DIR, "pfam_family_summary.csv"), DataFrame)

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

"""
    bootstrap_kl_se(gen_seqs, stored_seqs; n_boot=1000) -> (kl_point, kl_se)

Compute KL on full sample (point estimate) and bootstrap SE.
"""
function bootstrap_kl_se(gen_seqs, stored_seqs; n_boot=n_bootstrap)
    kl_point = aa_composition_kl(gen_seqs, stored_seqs)

    n = length(gen_seqs)
    boot_kls = Float64[]
    for b in 1:n_boot
        idx = rand(1:n, n)
        boot_seqs = gen_seqs[idx]
        push!(boot_kls, aa_composition_kl(boot_seqs, stored_seqs))
    end

    return (kl_point, std(boot_kls))
end

@info "Computing table values with SE …"

rows = Vector{NamedTuple}()

for fam in FAMILIES
    @info "  $(fam.name) ($(fam.id)) …"
    fam_dir = joinpath(PFAM_DATA_DIR, fam.id)

    # load stored sequences
    stored_data = parse_fasta(joinpath(fam_dir, "stored_sequences.fasta"))
    stored_seqs = [s[2] for s in stored_data]

    # rebuild memory matrix for PCA-space novelty
    char_mat = permutedims(hcat([collect(s) for s in stored_seqs]...))
    X̂, pca_model, L_out, d_full = build_memory_matrix(char_mat; pratio=R_pca)
    L = size(char_mat, 2)
    d, K = size(X̂)

    for (label, method_name) in [("sa_generation", "SA gen"),
                                  ("sa_retrieval", "SA ret")]
        fasta_file = joinpath(fam_dir, "$(label)_sequences.fasta")
        isfile(fasta_file) || continue
        gen_data = parse_fasta(fasta_file)
        gen_seqs = [s[2] for s in gen_data]

        # --- KL: bootstrap SE on full 150 samples ---
        kl_point, kl_se = bootstrap_kl_se(gen_seqs, stored_seqs)

        # --- Novelty (PCA cosine): per-chain SE ---
        # project all generated seqs to PCA space
        gen_pca = [project_to_pca(s, pca_model, L) for s in gen_seqs]

        n_total = length(gen_seqs)
        spc = min(samples_per_chain, n_total ÷ n_chains)
        chain_novelties = Float64[]
        chain_seqids = Float64[]

        for c in 1:n_chains
            si = (c - 1) * spc + 1
            ei = c * spc
            ei > n_total && break

            # novelty: 1 - max cosine similarity in PCA space
            chain_nov = mean(1.0 - nearest_cosine_similarity(gen_pca[j], X̂) for j in si:ei)
            push!(chain_novelties, chain_nov)

            # seqid: nearest sequence identity in AA space
            chain_sid = mean(nearest_sequence_identity(gen_seqs[j], stored_seqs) for j in si:ei)
            push!(chain_seqids, chain_sid)
        end

        n_ch = length(chain_novelties)
        push!(rows, (
            Family   = fam.name,
            Method   = method_name,
            KL_val   = kl_point,
            KL_se    = kl_se,
            Nov_val  = mean(chain_novelties),
            Nov_se   = std(chain_novelties) / sqrt(n_ch),
            SID_val  = mean(chain_seqids),
            SID_se   = std(chain_seqids) / sqrt(n_ch),
        ))
    end

    # --- Helper: evaluate a set of PCA-space samples ---
    function evaluate_pca_samples(pca_samps, method_label)
        # decode to AA sequences
        gen_seqs_bl = [decode_sample(ξ, pca_model, L) for ξ in pca_samps]

        # KL: bootstrap SE on decoded sequences
        kl_pt, kl_se_val = bootstrap_kl_se(gen_seqs_bl, stored_seqs)

        # per-chain novelty (PCA cosine) and seqid
        spc = samples_per_chain
        ch_nov = Float64[]
        ch_sid = Float64[]
        for c in 1:n_chains
            si = (c - 1) * spc + 1
            ei = c * spc
            ei > length(pca_samps) && break
            push!(ch_nov, mean(1.0 - nearest_cosine_similarity(pca_samps[j], X̂) for j in si:ei))
            push!(ch_sid, mean(nearest_sequence_identity(gen_seqs_bl[j], stored_seqs) for j in si:ei))
        end
        n_ch = length(ch_nov)
        push!(rows, (
            Family   = fam.name,
            Method   = method_label,
            KL_val   = kl_pt,
            KL_se    = kl_se_val,
            Nov_val  = mean(ch_nov),
            Nov_se   = std(ch_nov) / sqrt(n_ch),
            SID_val  = mean(ch_sid),
            SID_se   = std(ch_sid) / sqrt(n_ch),
        ))
    end

    # --- Bootstrap baseline: resample stored patterns ---
    Random.seed!(12345)
    bs_pca = [copy(X̂[:, rand(1:K)]) for _ in 1:S_total]
    evaluate_pca_samples(bs_pca, "BS")

    # --- Gaussian perturbation: stored pattern + matched noise ---
    fam_row = filter(r -> r.PfamID == fam.id, fam_summary)
    β_ret = Float64(fam_row.β_retrieval[1])
    σ_noise = sqrt(2 * α_step / β_ret)
    Random.seed!(12345)
    gp_pca = [X̂[:, rand(1:K)] .+ σ_noise .* randn(d) for _ in 1:S_total]
    evaluate_pca_samples(gp_pca, "GP")

    # --- Convex combination: Dirichlet-weighted average ---
    dirichlet = Dirichlet(K, 1.0)
    Random.seed!(12345)
    cc_pca = [X̂ * rand(dirichlet) for _ in 1:S_total]
    evaluate_pca_samples(cc_pca, "CC")
end

df = DataFrame(rows)
CSV.write(joinpath(PFAM_DATA_DIR, "table2_with_se.csv"), df)

# ── Print LaTeX table rows ────────────────────────────────────────────────────
println("\n% ── Table 2 with ±SE ──")
fmt(v, se) = string(round(v, digits=3), " \\pm ", round(se, digits=3))

for fam in FAMILIES
    sub = filter(r -> r.Family == fam.name, df)
    println("\\midrule")
    println("\\multirow{5}{*}{$(fam.name)}")
    for r in eachrow(sub)
        println(" & $(r.Method) & \$$(fmt(r.KL_val, r.KL_se))\$ & \$$(fmt(r.Nov_val, r.Nov_se))\$ & \$$(fmt(r.SID_val, r.SID_se))\$ \\\\")
    end
end

@info "Done. Saved table2_with_se.csv"
