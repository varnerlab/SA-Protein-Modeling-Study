#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Wall-Clock Timing Benchmark
#
# Measures generation time for 150 sequences across all 8 families for SA
# and HMM emit. Reports seconds per 150 sequences. EvoDiff/MSAT/Potts times
# are pulled from logged runs or estimated from per-sequence timing.
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment …"
const _EXPERIMENT_ROOT = @__DIR__
include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))

const TIMING_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "data")
const PFAM_DATA_DIR   = joinpath(_EXPERIMENT_ROOT, "..", "pfam-families", "data")

# ══════════════════════════════════════════════════════════════════════════════
# Configuration
# ══════════════════════════════════════════════════════════════════════════════

const FAMILIES = [
    (id="PF00076", name="RRM",            L=71),
    (id="PF00018", name="SH3",            L=48),
    (id="PF00397", name="WW",             L=31),
    (id="PF00014", name="Kunitz",         L=53),
    (id="PF00096", name="zf-C2H2",        L=23),
    (id="PF00595", name="PDZ",            L=83),
    (id="PF00069", name="Pkinase",        L=262),
    (id="PF00711", name="Defensin_beta",  L=36),
]

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
# Multi-chain SA sampler (timing version: no extra evaluation)
# ══════════════════════════════════════════════════════════════════════════════

function run_sa_chains_timed(X̂::Matrix{Float64}, β::Float64;
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
# Main benchmark
# ══════════════════════════════════════════════════════════════════════════════

function run_timing_benchmark()
    mkpath(TIMING_DATA_DIR)
    rows = Vector{NamedTuple}()

    for fam in FAMILIES
        @info "════════════════════════════════════════════════════════════"
        @info "  $(fam.name) ($(fam.id))"
        @info "════════════════════════════════════════════════════════════"

        # Load alignment
        sto_file = joinpath(PFAM_DATA_DIR, fam.id, "$(fam.id)_seed.sto")
        if !isfile(sto_file)
            @warn "  Missing alignment: $sto_file — skipping"
            continue
        end
        raw_seqs = parse_stockholm(sto_file)
        isempty(raw_seqs) && (raw_seqs = parse_fasta(sto_file))
        char_mat, _ = clean_alignment(raw_seqs)
        K, L = size(char_mat)

        # ── SA timing (includes PCA + phase transition + sampling + decoding) ──
        t_sa = @elapsed begin
            X̂, pca_model, L_out, d_full = build_memory_matrix(char_mat; pratio=R_pca)
            pt = find_entropy_inflection(X̂; α=α_mix)
            β_gen = Float64(max(round(Int, 2 * pt.β_star), 5))

            Random.seed!(42)
            sa_samps = run_sa_chains_timed(X̂, β_gen; base_seed=42)
            sa_seqs = [decode_sample(ξ, pca_model, L_out) for ξ in sa_samps]
        end
        @info "  SA: $(round(t_sa, digits=1))s for $(length(sa_samps)) sequences"

        push!(rows, (
            Family = fam.name,
            PfamID = fam.id,
            K = K,
            L = L,
            Method = "SA",
            Time_s = round(t_sa, digits=2),
            N_seqs = length(sa_samps),
            Hardware = "CPU (laptop)",
        ))

        # ── HMM emit timing ──
        hmm_available = true
        try
            run(pipeline(`hmmbuild -h`, devnull))
        catch
            hmm_available = false
            @warn "  HMMER3 not found, skipping HMM timing"
        end

        if hmm_available
            hmm_dir = joinpath(TIMING_DATA_DIR, fam.id)
            mkpath(hmm_dir)
            hmm_file = joinpath(hmm_dir, "$(fam.id).hmm")
            hmm_fasta = joinpath(hmm_dir, "hmm_sequences.fasta")

            t_hmm = @elapsed begin
                run(pipeline(`hmmbuild --amino $hmm_file $sto_file`, devnull))
                run(pipeline(`hmmemit -N $S --seed 42 -o $hmm_fasta $hmm_file`, devnull))
            end
            @info "  HMM: $(round(t_hmm, digits=2))s"

            push!(rows, (
                Family = fam.name,
                PfamID = fam.id,
                K = K,
                L = L,
                Method = "HMM emit",
                Time_s = round(t_hmm, digits=2),
                N_seqs = S,
                Hardware = "CPU (laptop)",
            ))
        end
    end

    # ══════════════════════════════════════════════════════════════════════════
    # Save and display
    # ══════════════════════════════════════════════════════════════════════════

    df = DataFrame(rows)
    CSV.write(joinpath(TIMING_DATA_DIR, "timing_benchmark.csv"), df)
    @info "Results saved."

    println("\n══════════════════════════════════════════════════════")
    println("WALL-CLOCK TIMING (150 sequences, CPU)")
    println("══════════════════════════════════════════════════════")
    pretty_table(df[:, [:Family, :Method, :L, :K, :Time_s]],
                 header=["Family", "Method", "L", "K", "Time (s)"])

    return df
end

# ── Run ──────────────────────────────────────────────────────────────────────
df_timing = run_timing_benchmark()
