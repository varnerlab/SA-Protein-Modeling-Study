#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Sampling Diagnostics (Reviewer #6)
#
# For each family, runs ULA and MALA chains with full trajectory recording.
# Computes:
#   - Autocorrelation function of the Hopfield energy E(ξ_t)
#   - Integrated autocorrelation time τ_int
#   - Effective sample size (ESS)
#   - MALA acceptance rate
#   - Validates burn-in and thinning choices
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment …"
const _EXPERIMENT_ROOT = @__DIR__

using Pkg
Pkg.activate(joinpath(_EXPERIMENT_ROOT, ".."))
include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))

using Printf, LinearAlgebra

# ── Hopfield energy ──────────────────────────────────────────────────────────

function hopfield_energy(ξ::Vector{Float64}, X̂::Matrix{Float64}, β::Float64)
    scores = X̂' * ξ  # K-vector
    lse = log(sum(exp.(β .* scores)))
    return 0.5 * dot(ξ, ξ) - lse / β
end

# ── Autocorrelation and ESS ──────────────────────────────────────────────────

function autocorrelation(x::Vector{Float64}; max_lag::Int=1000)
    n = length(x)
    x̄ = mean(x)
    var_x = var(x; corrected=false)
    var_x ≈ 0 && return ones(min(max_lag+1, n))

    max_lag = min(max_lag, n - 1)
    acf = zeros(max_lag + 1)
    for lag in 0:max_lag
        s = 0.0
        for i in 1:(n - lag)
            s += (x[i] - x̄) * (x[i + lag] - x̄)
        end
        acf[lag + 1] = s / (n * var_x)
    end
    return acf
end

function integrated_autocorrelation_time(acf::Vector{Float64})
    # Sokal's windowed estimator: sum until acf drops below threshold
    # Use the initial positive sequence estimator (Geyer 1992)
    τ = 1.0  # lag 0 contributes 1
    for lag in 1:(length(acf) - 1)
        if acf[lag + 1] < 0.05  # stop when autocorrelation is negligible
            break
        end
        τ += 2.0 * acf[lag + 1]
    end
    return τ
end

function ess_from_acf(n::Int, τ_int::Float64)
    return n / τ_int
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
T_chain = 5000
T_burn = 2000
thin = 100
n_diag_chains = 10  # fewer chains, but analyze each one fully
σ_init = 0.01

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

    # build memory matrix
    X̂, pca_model, _, d_full = build_memory_matrix(char_mat; pratio=0.95)
    d = size(X̂, 1)

    # find β*
    pt = find_entropy_inflection(X̂; α=α_step)
    β_star = pt.β_star
    β_gen = Float64(round(Int, max(2 * β_star, 5)))

    @info "  d=$d, K=$K, β*=$(round(β_star, digits=2)), β_gen=$β_gen"

    # run diagnostic chains
    Random.seed!(42)
    pattern_indices = StatsBase.sample(1:K, n_diag_chains, replace=(n_diag_chains > K))

    chain_τ_ula = Float64[]
    chain_ess_ula = Float64[]
    chain_τ_mala = Float64[]
    chain_ess_mala = Float64[]
    chain_ar_mala = Float64[]
    burn_converged = Int[]  # iteration where energy stabilizes

    for (c, k) in enumerate(pattern_indices)
        seed_c = 12345 + c
        Random.seed!(seed_c)
        ξ₀ = X̂[:, k] .+ σ_init .* randn(d)

        # ULA chain
        (_, Ξ_ula) = sample(X̂, ξ₀, T_chain; β=β_gen, α=α_step, seed=seed_c)

        # MALA chain
        Random.seed!(seed_c)
        ξ₀_mala = X̂[:, k] .+ σ_init .* randn(d)
        (_, Ξ_mala, ar) = mala_sample(X̂, ξ₀_mala, T_chain; β=β_gen, α=α_step, seed=seed_c)
        push!(chain_ar_mala, ar)

        # compute energy traces
        E_ula = [hopfield_energy(Ξ_ula[t, :], X̂, β_gen) for t in 1:(T_chain+1)]
        E_mala = [hopfield_energy(Ξ_mala[t, :], X̂, β_gen) for t in 1:(T_chain+1)]

        # autocorrelation on post-burn-in energy (ULA)
        E_post_burn_ula = E_ula[(T_burn+1):end]
        acf_ula = autocorrelation(E_post_burn_ula; max_lag=min(1000, length(E_post_burn_ula)-1))
        τ_ula = integrated_autocorrelation_time(acf_ula)
        push!(chain_τ_ula, τ_ula)
        push!(chain_ess_ula, ess_from_acf(length(E_post_burn_ula), τ_ula))

        # autocorrelation on post-burn-in energy (MALA)
        E_post_burn_mala = E_mala[(T_burn+1):end]
        acf_mala = autocorrelation(E_post_burn_mala; max_lag=min(1000, length(E_post_burn_mala)-1))
        τ_mala = integrated_autocorrelation_time(acf_mala)
        push!(chain_τ_mala, τ_mala)
        push!(chain_ess_mala, ess_from_acf(length(E_post_burn_mala), τ_mala))

        # estimate burn-in convergence: find iteration where running mean
        # of energy is within 1% of final mean
        E_final_mean = mean(E_ula[(T_burn+1):end])
        window = 200
        converged_at = T_chain  # default
        for t in window:T_chain
            running_mean = mean(E_ula[(t-window+1):t])
            if abs(running_mean - E_final_mean) / abs(E_final_mean) < 0.01
                converged_at = t
                break
            end
        end
        push!(burn_converged, converged_at)
    end

    mean_τ_ula = mean(chain_τ_ula)
    mean_ess_ula = mean(chain_ess_ula)
    mean_τ_mala = mean(chain_τ_mala)
    mean_ess_mala = mean(chain_ess_mala)
    mean_ar = mean(chain_ar_mala)
    mean_converge = mean(burn_converged)
    n_post_burn = T_chain - T_burn  # 3000
    n_thinned = length((T_burn+1):thin:T_chain)  # 30
    ess_per_thinned = mean_ess_ula / n_thinned  # ratio: how many independent samples per thinned sample

    @info @sprintf("  ULA:  τ_int=%.1f, ESS=%.0f (of %d post-burn-in), ESS/thin_sample=%.1f",
        mean_τ_ula, mean_ess_ula, n_post_burn, ess_per_thinned)
    @info @sprintf("  MALA: τ_int=%.1f, ESS=%.0f, accept=%.4f",
        mean_τ_mala, mean_ess_mala, mean_ar)
    @info @sprintf("  Burn-in convergence: mean iter %.0f (budget: %d)", mean_converge, T_burn)

    push!(results, (
        family=family_name, pfam_id=pfam_id, d=d, K=K,
        beta_gen=β_gen,
        τ_int_ula=mean_τ_ula, ess_ula=mean_ess_ula,
        τ_int_mala=mean_τ_mala, ess_mala=mean_ess_mala,
        mala_accept=mean_ar,
        n_post_burn=n_post_burn,
        n_thinned=n_thinned,
        ess_per_thinned=ess_per_thinned,
        burn_converge_iter=mean_converge,
        τ_int_ula_std=std(chain_τ_ula),
        ess_ula_std=std(chain_ess_ula),
    ))
end

# ── Write CSV ─────────────────────────────────────────────────────────────────

df = DataFrame(results)
csv_path = joinpath(output_dir, "sampling_diagnostics.csv")
CSV.write(csv_path, df)
@info "Results saved to $csv_path"

# ── Print summary table ──────────────────────────────────────────────────────

println("\n" * "="^110)
println("SAMPLING DIAGNOSTICS")
println("="^110)
@printf("%-15s %5s %5s %7s %8s %8s %8s %8s %8s %10s\n",
    "Family", "d", "K", "β_gen", "τ_ULA", "ESS_ULA", "τ_MALA", "ESS_MALA", "MALA_AR", "Burn_conv")
println("-"^110)

for r in results
    @printf("%-15s %5d %5d %7.0f %8.1f %8.0f %8.1f %8.0f %8.4f %10.0f\n",
        r.family, r.d, r.K, r.beta_gen,
        r.τ_int_ula, r.ess_ula,
        r.τ_int_mala, r.ess_mala,
        r.mala_accept, r.burn_converge_iter)
end

# ── Thinning justification ────────────────────────────────────────────────────

println("\n\nThinning justification (thin=$thin):")
println("-"^80)
@printf("%-15s %10s %10s %12s %15s\n",
    "Family", "τ_int", "thin", "thin/τ_int", "ESS_per_thin_samp")
println("-"^80)
for r in results
    @printf("%-15s %10.1f %10d %12.1f %15.1f\n",
        r.family, r.τ_int_ula, thin, thin / r.τ_int_ula, r.ess_per_thinned)
end

println("\n\nBurn-in justification (burn=$T_burn):")
println("-"^80)
@printf("%-15s %12s %10s %15s\n",
    "Family", "Conv_iter", "Burn-in", "Margin(×)")
println("-"^80)
for r in results
    margin = T_burn / r.burn_converge_iter
    @printf("%-15s %12.0f %10d %15.1f\n",
        r.family, r.burn_converge_iter, T_burn, margin)
end

@info "\nDone."
