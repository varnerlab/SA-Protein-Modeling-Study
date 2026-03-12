#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# HMM Baseline: build profile HMMs from Pfam seed alignments and emit sequences
#
# Uses HMMER3 (hmmbuild + hmmemit). Requires hmmer to be installed.
# This is the classical bioinformatics baseline for protein sequence generation.
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment …"
const _EXPERIMENT_ROOT = @__DIR__
include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))

const BASELINE_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "data")
const PFAM_DATA_DIR     = joinpath(_EXPERIMENT_ROOT, "..", "pfam-families", "data")

# Family definitions (same as pfam-families experiment)
const FAMILIES = [
    (id="PF00076", name="RRM"),
    (id="PF00018", name="SH3"),
    (id="PF00397", name="WW"),
    (id="PF00014", name="Kunitz"),
    (id="PF00096", name="zf-C2H2"),
    (id="PF00595", name="PDZ"),
    (id="PF00069", name="Pkinase"),
]

const S = 150  # number of sequences to generate (match SA output)

# ══════════════════════════════════════════════════════════════════════════════
# Check HMMER3 installation
# ══════════════════════════════════════════════════════════════════════════════

function check_hmmer()
    try
        run(pipeline(`hmmbuild -h`, devnull))
        run(pipeline(`hmmemit -h`, devnull))
        @info "HMMER3 found."
        return true
    catch
        @error "HMMER3 not found. Install with: brew install hmmer"
        return false
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# Main: build HMM and emit sequences for each family
# ══════════════════════════════════════════════════════════════════════════════

function run_hmm_baseline()
    check_hmmer() || return

    rows = Vector{NamedTuple}()

    for fam in FAMILIES
        @info "════════════════════════════════════════════════════════════"
        @info "  $(fam.name) ($(fam.id))"
        @info "════════════════════════════════════════════════════════════"

        fam_data_dir = joinpath(BASELINE_DATA_DIR, fam.id)
        mkpath(fam_data_dir)

        # locate the Stockholm alignment from the pfam-families experiment
        sto_file = joinpath(PFAM_DATA_DIR, fam.id, "$(fam.id)_seed.sto")
        if !isfile(sto_file)
            @warn "  Missing alignment: $sto_file — skipping"
            continue
        end

        # load stored sequences for evaluation
        stored_fasta = joinpath(PFAM_DATA_DIR, fam.id, "stored_sequences.fasta")
        if !isfile(stored_fasta)
            @warn "  Missing stored sequences: $stored_fasta — skipping"
            continue
        end
        stored_data = parse_fasta(stored_fasta)
        stored_seqs = [s[2] for s in stored_data]

        # build memory matrix for PCA-space novelty
        raw_seqs = parse_stockholm(sto_file)
        char_mat, _ = clean_alignment(raw_seqs)
        K, L = size(char_mat)
        X̂, pca_model, L, d_full = build_memory_matrix(char_mat; pratio=0.95)

        # ── Step 1: Build HMM ────────────────────────────────────────────
        hmm_file = joinpath(fam_data_dir, "$(fam.id).hmm")
        @info "  Building HMM …"
        run(pipeline(`hmmbuild --amino $hmm_file $sto_file`, devnull))
        @info "  HMM built: $hmm_file"

        # ── Step 2: Emit sequences ───────────────────────────────────────
        hmm_fasta = joinpath(fam_data_dir, "hmm_sequences.fasta")
        @info "  Emitting $S sequences …"
        run(`hmmemit -N $S --seed 42 -o $hmm_fasta $hmm_file`)
        @info "  Saved: $hmm_fasta"

        # ── Step 3: Parse and clean emitted sequences ────────────────────
        hmm_data = parse_fasta(hmm_fasta)
        hmm_seqs_raw = [s[2] for s in hmm_data]

        # hmmemit can produce variable-length sequences with insertions
        # Align to the HMM to get fixed-length sequences, or just trim/pad
        # For fair comparison, we strip non-standard characters and truncate/pad to L
        hmm_seqs = String[]
        for seq in hmm_seqs_raw
            # keep only standard amino acids
            cleaned = filter(c -> c in AA_ALPHABET, uppercase(seq))
            # truncate or pad to alignment length L
            if length(cleaned) >= L
                push!(hmm_seqs, cleaned[1:L])
            else
                # pad with most common residue (or gap-equivalent)
                push!(hmm_seqs, cleaned * repeat("A", L - length(cleaned)))
            end
        end
        @info "  Processed $(length(hmm_seqs)) HMM sequences (aligned to L=$L)"

        # ── Step 4: Encode into PCA space for novelty computation ────────
        # One-hot encode the HMM sequences, project through PCA
        hmm_pca_samples = Vector{Vector{Float64}}()
        for seq in hmm_seqs
            # build a 1-row char matrix for encoding
            char_vec = collect(seq)
            x_onehot = zeros(Float64, N_AA * L)
            for pos in 1:min(L, length(char_vec))
                idx = get(AA_TO_IDX, char_vec[pos], 0)
                if idx > 0
                    x_onehot[(pos-1)*N_AA + idx] = 1.0
                end
            end
            # project to PCA space
            ξ_pca = vec(MultivariateStats.transform(pca_model, x_onehot))
            push!(hmm_pca_samples, ξ_pca)
        end

        # ── Step 5: Evaluate ─────────────────────────────────────────────
        @info "  Evaluating …"

        # KL divergence
        kl = aa_composition_kl(hmm_seqs, stored_seqs)

        # Novelty (PCA-space cosine)
        novelties = [sample_novelty(ξ, X̂) for ξ in hmm_pca_samples]
        nov_mean = mean(novelties)

        # Sequence identity
        seqids = [nearest_sequence_identity(s, stored_seqs) for s in hmm_seqs]
        sid_mean = mean(seqids)

        # Valid residue fraction
        valid_frac = mean(valid_residue_fraction(s) for s in hmm_seqs)

        # Per-chain SE (split 150 sequences into 30 groups of 5)
        n_chains = 30
        spc = 5
        chain_kls = Float64[]
        chain_novs = Float64[]
        chain_sids = Float64[]
        for c in 1:n_chains
            i1 = (c-1)*spc + 1
            i2 = c*spc
            i2 > length(hmm_seqs) && break
            chain_seqs = hmm_seqs[i1:i2]
            chain_pca  = hmm_pca_samples[i1:i2]
            push!(chain_kls, aa_composition_kl(chain_seqs, stored_seqs))
            push!(chain_novs, mean(sample_novelty(ξ, X̂) for ξ in chain_pca))
            push!(chain_sids, mean(nearest_sequence_identity(s, stored_seqs) for s in chain_seqs))
        end

        n_ch = length(chain_kls)
        kl_se  = n_ch > 1 ? std(chain_kls)  / sqrt(n_ch) : 0.0
        nov_se = n_ch > 1 ? std(chain_novs) / sqrt(n_ch) : 0.0
        sid_se = n_ch > 1 ? std(chain_sids) / sqrt(n_ch) : 0.0

        @info "  Results: KL=$(round(kl, digits=4)), Nov=$(round(nov_mean, digits=4)), SeqID=$(round(sid_mean, digits=4)), Valid=$(round(valid_frac, digits=4))"

        push!(rows, (
            Family   = fam.name,
            PfamID   = fam.id,
            Method   = "HMM emit",
            KL_mean  = kl,
            KL_se    = kl_se,
            Nov_mean = nov_mean,
            Nov_se   = nov_se,
            SID_mean = sid_mean,
            SID_se   = sid_se,
            ValidAA  = valid_frac,
        ))
    end

    # ══════════════════════════════════════════════════════════════════════════
    # Save results
    # ══════════════════════════════════════════════════════════════════════════

    df = DataFrame(rows)
    results_dir = joinpath(_EXPERIMENT_ROOT, "results")
    mkpath(results_dir)
    CSV.write(joinpath(results_dir, "hmm_baseline_results.csv"), df)
    @info "Results saved to results/hmm_baseline_results.csv"

    pretty_table(df[:, [:Family, :Method, :KL_mean, :KL_se, :Nov_mean, :Nov_se, :SID_mean, :SID_se]])

    return df
end

# ── Run ──────────────────────────────────────────────────────────────────────
df_hmm = run_hmm_baseline()
