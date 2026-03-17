#!/usr/bin/env julia
# Evaluate Potts (plmDCA) baseline sequences using the same metrics as other baselines

@info "Loading environment …"
const _EXPERIMENT_ROOT = @__DIR__
include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))

const BASELINE_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "data")
const PFAM_DATA_DIR     = joinpath(_EXPERIMENT_ROOT, "..", "pfam-families", "data")

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

const n_chains_eval = 30
const spc_eval = 5

function evaluate_sequences(gen_seqs::Vector{String}, stored_seqs::Vector{String},
                             X̂::Matrix{Float64}, pca_model, L::Int)
    gen_pca = Vector{Vector{Float64}}()
    for seq in gen_seqs
        x_onehot = zeros(Float64, N_AA * L)
        for pos in 1:min(L, length(seq))
            idx = get(AA_TO_IDX, seq[pos], 0)
            if idx > 0
                x_onehot[(pos-1)*N_AA + idx] = 1.0
            end
        end
        ξ_pca = vec(MultivariateStats.transform(pca_model, x_onehot))
        push!(gen_pca, ξ_pca)
    end

    kl_full = aa_composition_kl(gen_seqs, stored_seqs)
    nov_full = mean(sample_novelty(ξ, X̂) for ξ in gen_pca)
    sid_full = mean(nearest_sequence_identity(s, stored_seqs) for s in gen_seqs)

    n_total = length(gen_seqs)
    spc = min(spc_eval, n_total ÷ n_chains_eval)

    chain_kls = Float64[]
    chain_novs = Float64[]
    chain_sids = Float64[]
    for c in 1:n_chains_eval
        i1 = (c-1)*spc + 1
        i2 = c*spc
        i2 > n_total && break
        ch_seqs = gen_seqs[i1:i2]
        ch_pca  = gen_pca[i1:i2]
        push!(chain_kls, aa_composition_kl(ch_seqs, stored_seqs))
        push!(chain_novs, mean(sample_novelty(ξ, X̂) for ξ in ch_pca))
        push!(chain_sids, mean(nearest_sequence_identity(s, stored_seqs) for s in ch_seqs))
    end

    n_ch = length(chain_kls)
    return (
        KL_mean  = kl_full,
        KL_se    = n_ch > 1 ? std(chain_kls)  / sqrt(n_ch) : 0.0,
        Nov_mean = nov_full,
        Nov_se   = n_ch > 1 ? std(chain_novs) / sqrt(n_ch) : 0.0,
        SID_mean = sid_full,
        SID_se   = n_ch > 1 ? std(chain_sids) / sqrt(n_ch) : 0.0,
        ValidAA  = mean(valid_residue_fraction(s) for s in gen_seqs),
    )
end

function load_and_clean_seqs(fasta_path::String, L::Int)
    data = parse_fasta(fasta_path)
    seqs = String[]
    for (_, raw_seq) in data
        cleaned = filter(c -> c in AA_ALPHABET, uppercase(raw_seq))
        if length(cleaned) >= L
            push!(seqs, cleaned[1:L])
        elseif length(cleaned) > 0
            push!(seqs, cleaned * repeat("A", L - length(cleaned)))
        end
    end
    return seqs
end

function main()
    println("\n" * "="^70)
    println("Potts (plmDCA) Baseline Evaluation")
    println("="^70)

    for fam in FAMILIES
        println("\n--- $(fam.name) ($(fam.id)) ---")

        # Load stored sequences
        stored_fasta = joinpath(PFAM_DATA_DIR, fam.id, "stored_sequences.fasta")
        stored_data = parse_fasta(stored_fasta)
        stored_seqs = [s[2] for s in stored_data]

        # Build memory matrix for PCA-space novelty
        sto_file = joinpath(PFAM_DATA_DIR, fam.id, "$(fam.id)_seed.sto")
        raw_seqs = parse_stockholm(sto_file)
        char_mat, _ = clean_alignment(raw_seqs)
        K, L = size(char_mat)
        X̂, pca_model, L, d_full = build_memory_matrix(char_mat; pratio=0.95)
        println("  Family: K=$K, L=$L, d=$d_full")

        # Evaluate Potts sequences
        potts_path = joinpath(BASELINE_DATA_DIR, fam.id, "potts_sequences.fasta")
        if !isfile(potts_path)
            @warn "  Potts sequences not found at $potts_path"
            continue
        end
        potts_seqs = load_and_clean_seqs(potts_path, L)
        println("  Loaded $(length(potts_seqs)) Potts sequences")

        pm = evaluate_sequences(potts_seqs, stored_seqs, X̂, pca_model, L)
        println("  Potts plmDCA:")
        println("    KL  = $(round(pm.KL_mean, digits=4)) ± $(round(pm.KL_se, digits=4))")
        println("    Nov = $(round(pm.Nov_mean, digits=4)) ± $(round(pm.Nov_se, digits=4))")
        println("    SID = $(round(pm.SID_mean, digits=4)) ± $(round(pm.SID_se, digits=4))")
        println("    ValidAA = $(round(pm.ValidAA, digits=4))")

        # Also evaluate SA generation for direct comparison
        sa_path = joinpath(PFAM_DATA_DIR, fam.id, "sa_generation_sequences.fasta")
        sa_data = parse_fasta(sa_path)
        sa_seqs = [s[2] for s in sa_data]
        sm = evaluate_sequences(sa_seqs, stored_seqs, X̂, pca_model, L)
        println("  SA generation:")
        println("    KL  = $(round(sm.KL_mean, digits=4)) ± $(round(sm.KL_se, digits=4))")
        println("    Nov = $(round(sm.Nov_mean, digits=4)) ± $(round(sm.Nov_se, digits=4))")
        println("    SID = $(round(sm.SID_mean, digits=4)) ± $(round(sm.SID_se, digits=4))")
        println("    ValidAA = $(round(sm.ValidAA, digits=4))")

        # Also show HMM for reference
        hmm_path = joinpath(BASELINE_DATA_DIR, fam.id, "hmm_sequences.fasta")
        if isfile(hmm_path)
            hmm_seqs = load_and_clean_seqs(hmm_path, L)
            hm = evaluate_sequences(hmm_seqs, stored_seqs, X̂, pca_model, L)
            println("  HMM emit:")
            println("    KL  = $(round(hm.KL_mean, digits=4)) ± $(round(hm.KL_se, digits=4))")
            println("    Nov = $(round(hm.Nov_mean, digits=4)) ± $(round(hm.Nov_se, digits=4))")
            println("    SID = $(round(hm.SID_mean, digits=4)) ± $(round(hm.SID_se, digits=4))")
            println("    ValidAA = $(round(hm.ValidAA, digits=4))")
        end
    end
end

main()
