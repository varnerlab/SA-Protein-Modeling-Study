#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Analytical β* Investigation
#
# Part 1: Derive τ*(K) — the dimensionless critical temperature for K i.i.d.
#          Gaussian similarities. The entropy inflection of the softmax over
#          K standard normals defines τ*(K) such that β* = τ*(K)/σ where
#          σ = √Var(m_k^T ξ) is the similarity std dev.
#
# Part 2: Compute the similarity distribution moments for each protein family
#          (variance, skewness, kurtosis) and test whether β* = τ*(K)/σ_family.
#
# Part 3: Check if higher moments (kurtosis) explain residuals.
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment …"
const _EXPERIMENT_ROOT = @__DIR__

using Pkg
Pkg.activate(joinpath(_EXPERIMENT_ROOT, ".."))

include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))

using Printf
using LinearAlgebra
using MultivariateStats

const PFAM_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "..", "pfam-families", "data")
const OUTPUT_DIR = joinpath(_EXPERIMENT_ROOT, "data")
const FIGS_DIR = joinpath(_EXPERIMENT_ROOT, "figs")

const FAMILIES = [
    (id="PF00076", name="RRM"),
    (id="PF00018", name="SH3"),
    (id="PF00397", name="WW"),
    (id="PF00014", name="Kunitz"),
    (id="PF00096", name="zf-C2H2"),
    (id="PF00595", name="PDZ"),
    (id="PF00069", name="Pkinase"),
]

# ══════════════════════════════════════════════════════════════════════════════
# PART 1: Dimensionless critical temperature τ*(K) for Gaussian similarities
# ══════════════════════════════════════════════════════════════════════════════

"""
Compute attention entropy for softmax over scores z at inverse temperature τ.
H(τ) = log Z - τ <z>_p where Z = Σ exp(τ z_k), <z>_p = Σ p_k z_k
"""
function softmax_entropy(z::Vector{Float64}, τ::Float64)
    τz = τ .* z
    # Numerically stable log-sum-exp
    mx = maximum(τz)
    log_Z = mx + log(sum(exp.(τz .- mx)))
    p = exp.(τz .- mx) ./ sum(exp.(τz .- mx))
    H = log_Z - τ * dot(p, z)
    return H
end

"""
Compute softmax-weighted variance of scores z at inverse temperature τ.
"""
function softmax_variance(z::Vector{Float64}, τ::Float64)
    τz = τ .* z
    mx = maximum(τz)
    p = exp.(τz .- mx) ./ sum(exp.(τz .- mx))
    mean_z = dot(p, z)
    var_z = dot(p, z.^2) - mean_z^2
    return var_z
end

"""
Find τ* for a single realization of K i.i.d. scores by entropy inflection.
"""
function find_tau_star(z::Vector{Float64}; n_tau=200, τ_range=(0.01, 50.0))
    τs = exp.(range(log(τ_range[1]), log(τ_range[2]), length=n_tau))
    Hs = [softmax_entropy(z, τ) for τ in τs]

    # Second derivative in log-τ space (same as the main code)
    log_τs = log.(τs)
    dH = diff(Hs) ./ diff(log_τs)
    d2H = diff(dH) ./ diff(log_τs[1:end-1])

    # Find minimum of d²H (most negative = steepest entropy drop)
    idx = argmin(d2H)
    return τs[idx+1]
end

"""
Compute τ*(K) by averaging over many Gaussian realizations.
"""
function compute_tau_star_gaussian(K::Int; n_realizations=500, rng=Random.MersenneTwister(42))
    τ_stars = Float64[]
    for _ in 1:n_realizations
        z = randn(rng, K)  # K i.i.d. N(0,1)
        τ_star = find_tau_star(z)
        push!(τ_stars, τ_star)
    end
    return mean(τ_stars), std(τ_stars) / sqrt(n_realizations)
end

function part1_gaussian_inflection()
    @info "\n══════════════════════════════════════════════════════════"
    @info "  PART 1: τ*(K) for Gaussian similarities"
    @info "══════════════════════════════════════════════════════════\n"

    K_values = [10, 20, 30, 37, 44, 50, 55, 68, 99, 100, 151, 200, 300, 400, 420, 500, 1000]
    tau_results = Vector{NamedTuple}()

    for K in K_values
        τ_mean, τ_se = compute_tau_star_gaussian(K)
        @info @sprintf("  K=%4d:  τ* = %.4f ± %.4f,  τ*/√(2 log K) = %.4f",
            K, τ_mean, τ_se, τ_mean / sqrt(2*log(K)))
        push!(tau_results, (K=K, τ_star=τ_mean, τ_se=τ_se,
                            sqrt_2logK=sqrt(2*log(K)),
                            ratio=τ_mean/sqrt(2*log(K))))
    end

    # Fit: τ*(K) = a + b √(2 log K)  or  τ*(K) = c √(2 log K)
    X_reg = hcat(ones(length(tau_results)), [r.sqrt_2logK for r in tau_results])
    y_reg = [r.τ_star for r in tau_results]
    coeffs = X_reg \ y_reg
    y_pred = X_reg * coeffs
    SS_res = sum((y_reg .- y_pred).^2)
    SS_tot = sum((y_reg .- mean(y_reg)).^2)
    R2 = 1 - SS_res / SS_tot

    @info @sprintf("\n  Fit: τ*(K) = %.4f + %.4f √(2 log K),  R² = %.4f", coeffs[1], coeffs[2], R2)

    # Also fit proportional: τ* = c √(2 log K)
    X_prop = reshape([r.sqrt_2logK for r in tau_results], :, 1)
    c_prop = (X_prop' * X_prop) \ (X_prop' * y_reg)
    y_prop = X_prop * c_prop
    SS_res_p = sum((y_reg .- y_prop).^2)
    R2_p = 1 - SS_res_p / SS_tot

    @info @sprintf("  Fit: τ*(K) = %.4f √(2 log K),  R² = %.4f", c_prop[1], R2_p)

    return tau_results, coeffs
end

# ══════════════════════════════════════════════════════════════════════════════
# PART 2: Similarity distribution moments for protein families
# ══════════════════════════════════════════════════════════════════════════════

"""
Compute similarity distribution moments by sampling random queries on S^{d-1}.
For each random query ξ ∈ S^{d-1}, compute {m_k^T ξ}_{k=1}^K and collect stats.
"""
function similarity_moments(X̂::Matrix{Float64}; n_probes=1000, rng=Random.MersenneTwister(123))
    d, K = size(X̂)

    all_vars = Float64[]
    all_kurt = Float64[]
    all_skew = Float64[]

    for _ in 1:n_probes
        ξ = randn(rng, d)
        ξ ./= norm(ξ)  # uniform on S^{d-1}

        e = X̂' * ξ  # K similarities
        μ_e = mean(e)
        σ²_e = var(e; corrected=false)
        push!(all_vars, σ²_e)

        # Standardized moments
        if σ²_e > 1e-12
            σ_e = sqrt(σ²_e)
            z = (e .- μ_e) ./ σ_e
            push!(all_skew, mean(z.^3))
            push!(all_kurt, mean(z.^4) - 3.0)  # excess kurtosis
        end
    end

    return (
        σ² = mean(all_vars),         # average variance of similarities
        σ²_std = std(all_vars),      # std of variance across probes
        σ = sqrt(mean(all_vars)),
        skewness = mean(all_skew),
        kurtosis = mean(all_kurt),   # excess kurtosis (0 for Gaussian)
        d_effective = 1.0 / mean(all_vars),  # 1/σ² ≈ d for random
    )
end

"""
Also compute τ*(K) directly from the actual similarity distribution
(not assuming Gaussian). Sample probes, compute entropy curves, find inflection.
"""
function empirical_tau_star(X̂::Matrix{Float64}; n_probes=200, rng=Random.MersenneTwister(456))
    d, K = size(X̂)
    τ_stars = Float64[]

    for _ in 1:n_probes
        ξ = randn(rng, d)
        ξ ./= norm(ξ)
        e = X̂' * ξ  # K similarities

        # Normalize to unit variance for τ comparison
        σ_e = std(e; corrected=false)
        z = e ./ (σ_e + 1e-12)  # standardized: Var(z) ≈ 1

        τ_star = find_tau_star(z)
        push!(τ_stars, τ_star)
    end

    return mean(τ_stars), std(τ_stars) / sqrt(length(τ_stars))
end

function part2_protein_moments()
    @info "\n══════════════════════════════════════════════════════════"
    @info "  PART 2: Similarity moments for protein families"
    @info "══════════════════════════════════════════════════════════\n"

    results = Vector{NamedTuple}()

    for fam in FAMILIES
        sto_path = joinpath(PFAM_DATA_DIR, fam.id, "$(fam.id)_seed.sto")
        if !isfile(sto_path); continue; end

        raw_seqs = parse_stockholm(sto_path)
        char_mat, _ = clean_alignment(raw_seqs)
        X̂, pca_model, L, d_full = build_memory_matrix(char_mat; pratio=0.95)
        d, K = size(X̂)

        # Empirical β*
        inflection = find_entropy_inflection(X̂)
        β_star = inflection.β_star

        # Similarity moments
        mom = similarity_moments(X̂)

        # Empirical τ* from actual (non-Gaussian) similarity distribution
        τ_emp, τ_emp_se = empirical_tau_star(X̂)

        @info @sprintf("  %s: d=%d, K=%d, β*=%.3f", fam.name, d, K, β_star)
        @info @sprintf("    σ²=%.6f (1/d=%.6f, ratio=%.3f), σ=%.4f",
            mom.σ², 1.0/d, mom.σ² * d, mom.σ)
        @info @sprintf("    skewness=%.4f, excess kurtosis=%.4f", mom.skewness, mom.kurtosis)
        @info @sprintf("    d_eff(1/σ²)=%.1f vs d=%d", mom.d_effective, d)
        @info @sprintf("    τ*_empirical=%.4f ± %.4f", τ_emp, τ_emp_se)
        @info @sprintf("    β*_predicted(τ*/σ)=%.3f vs β*_actual=%.3f",
            τ_emp / mom.σ, β_star)

        push!(results, (
            Family = fam.name, K = K, d = d, sqrt_d = sqrt(d),
            β_star = β_star,
            σ² = mom.σ², σ = mom.σ,
            one_over_d = 1.0/d,
            σ²_times_d = mom.σ² * d,  # should be ~1 for random
            skewness = mom.skewness,
            kurtosis = mom.kurtosis,
            d_effective = mom.d_effective,
            τ_star_emp = τ_emp,
            τ_star_se = τ_emp_se,
            β_predicted = τ_emp / mom.σ,
        ))
    end

    # Also do WW scaling subsamples
    @info "\n  WW scaling subsamples:"
    ww_sto = joinpath(PFAM_DATA_DIR, "PF00397", "PF00397_seed.sto")
    if isfile(ww_sto)
        raw_seqs = parse_stockholm(ww_sto)
        full_char_mat, _ = clean_alignment(raw_seqs)
        K_full = size(full_char_mat, 1)

        for K_sub in [20, 50, 100, 200, 400]
            rng = Random.MersenneTwister(1000 * K_sub + 1)  # rep=1
            K_actual = min(K_sub, K_full)
            indices = sort(Random.randperm(rng, K_full)[1:K_actual])
            char_sub = full_char_mat[indices, :]

            X̂, _, _, _ = build_memory_matrix(char_sub; pratio=0.95)
            d, K = size(X̂)
            inflection = find_entropy_inflection(X̂)
            β_star = inflection.β_star
            mom = similarity_moments(X̂)
            τ_emp, τ_emp_se = empirical_tau_star(X̂)

            @info @sprintf("  WW K=%d: d=%d, σ²=%.6f, σ²·d=%.3f, kurt=%.3f, τ*=%.4f, β*_pred=%.3f vs %.3f",
                K_sub, d, mom.σ², mom.σ²*d, mom.kurtosis, τ_emp, τ_emp/mom.σ, β_star)

            push!(results, (
                Family = "WW_K$(K_sub)", K = K, d = d, sqrt_d = sqrt(d),
                β_star = β_star,
                σ² = mom.σ², σ = mom.σ,
                one_over_d = 1.0/d,
                σ²_times_d = mom.σ² * d,
                skewness = mom.skewness,
                kurtosis = mom.kurtosis,
                d_effective = mom.d_effective,
                τ_star_emp = τ_emp,
                τ_star_se = τ_emp_se,
                β_predicted = τ_emp / mom.σ,
            ))
        end
    end

    return results
end

# ══════════════════════════════════════════════════════════════════════════════
# PART 3: Combined analysis and figures
# ══════════════════════════════════════════════════════════════════════════════

function part3_combined(tau_gauss, tau_coeffs, protein_results)
    @info "\n══════════════════════════════════════════════════════════"
    @info "  PART 3: Combined analysis"
    @info "══════════════════════════════════════════════════════════\n"

    # ── Test: β* = τ*(K)/σ_family ────────────────────────────────────────
    # Use the Gaussian τ*(K) function to predict β*
    @info @sprintf("  %-12s %4s %4s %7s %7s %7s %7s %7s %7s",
        "Family", "K", "d", "σ", "τ*_G", "τ*_emp", "β*_G/σ", "β*_e/σ", "β*_act")
    @info "  " * "-"^80

    β_predicted_gauss = Float64[]
    β_predicted_emp = Float64[]
    β_actual = Float64[]

    for r in protein_results
        # Get Gaussian τ*(K) from our fit
        τ_gauss = tau_coeffs[1] + tau_coeffs[2] * sqrt(2*log(r.K))
        β_gauss = τ_gauss / r.σ
        β_emp = r.β_predicted  # τ*_emp / σ

        push!(β_predicted_gauss, β_gauss)
        push!(β_predicted_emp, β_emp)
        push!(β_actual, r.β_star)

        @info @sprintf("  %-12s %4d %4d %7.4f %7.4f %7.4f %7.3f %7.3f %7.3f",
            r.Family, r.K, r.d, r.σ, τ_gauss, r.τ_star_emp,
            β_gauss, β_emp, r.β_star)
    end

    # Correlations
    r_gauss = cor(β_predicted_gauss, β_actual)
    r_emp = cor(β_predicted_emp, β_actual)

    SS_tot = sum((β_actual .- mean(β_actual)).^2)
    R2_gauss = 1 - sum((β_actual .- β_predicted_gauss).^2) / SS_tot
    R2_emp = 1 - sum((β_actual .- β_predicted_emp).^2) / SS_tot

    @info @sprintf("\n  β* = τ*_Gaussian(K) / σ_family:  r = %.3f,  R² = %.3f", r_gauss, R2_gauss)
    @info @sprintf("  β* = τ*_empirical / σ_family:    r = %.3f,  R² = %.3f", r_emp, R2_emp)

    # ── Check kurtosis as correction ─────────────────────────────────────
    kurtosis_vals = [r.kurtosis for r in protein_results]
    residuals_gauss = β_actual .- β_predicted_gauss
    r_kurt = cor(kurtosis_vals, residuals_gauss)
    @info @sprintf("  Kurtosis-residual correlation: r = %.3f", r_kurt)

    # ── Regression: β* = a/σ + b (one-parameter from similarity variance)
    σ_inv = [1.0/r.σ for r in protein_results]
    X_σ = hcat(ones(length(σ_inv)), σ_inv)
    c_σ = X_σ \ β_actual
    y_σ = X_σ * c_σ
    R2_σ = 1 - sum((β_actual .- y_σ).^2) / SS_tot
    @info @sprintf("  β* = %.3f + %.4f / σ:  R² = %.3f", c_σ[1], c_σ[2], R2_σ)

    # ── Figures ──────────────────────────────────────────────────────────
    generate_combined_figures(tau_gauss, tau_coeffs, protein_results,
                              β_predicted_gauss, β_predicted_emp, β_actual)

    # ── Save data ────────────────────────────────────────────────────────
    df = DataFrame(protein_results)
    df.β_pred_gauss = β_predicted_gauss
    df.β_pred_emp = β_predicted_emp
    CSV.write(joinpath(OUTPUT_DIR, "analytical_beta_star.csv"), df)
    @info "Saved: $(joinpath(OUTPUT_DIR, "analytical_beta_star.csv"))"
end

function generate_combined_figures(tau_gauss, tau_coeffs, protein_results,
                                    β_pred_gauss, β_pred_emp, β_actual)
    mkpath(FIGS_DIR)

    # ── Panel A: τ*(K) vs √(2 log K) ────────────────────────────────────
    K_vals = [r.K for r in tau_gauss]
    τ_vals = [r.τ_star for r in tau_gauss]
    x_vals = [r.sqrt_2logK for r in tau_gauss]

    p1 = scatter(x_vals, τ_vals,
                 xlabel="√(2 log K)", ylabel="τ*(K)",
                 title="(A) Gaussian τ*(K)", titlefontsize=9,
                 legend=:topleft, legendfontsize=7,
                 markersize=5, markercolor=RGB(0.122, 0.467, 0.706),
                 label="Numerical", grid=true)
    x_line = range(minimum(x_vals)*0.9, maximum(x_vals)*1.1, length=100)
    y_line = tau_coeffs[1] .+ tau_coeffs[2] .* x_line
    plot!(p1, x_line, y_line, color=:black, linewidth=1.5,
          label=@sprintf("%.3f + %.3f√(2logK)", tau_coeffs[1], tau_coeffs[2]))

    # ── Panel B: σ² vs 1/d (is similarity variance ≈ 1/d?) ──────────────
    d_vals = [r.d for r in protein_results]
    σ²_vals = [r.σ² for r in protein_results]
    inv_d = [1.0/r.d for r in protein_results]

    p2 = scatter(inv_d, σ²_vals,
                 xlabel="1/d", ylabel="σ² = Var(m_k^T ξ)",
                 title="(B) Similarity Variance", titlefontsize=9,
                 legend=:topleft, legendfontsize=7,
                 markersize=5, markercolor=RGB(0.173, 0.627, 0.173),
                 label="Protein families", grid=true)
    lo = minimum(inv_d)*0.8; hi = maximum(inv_d)*1.2
    plot!(p2, [lo, hi], [lo, hi], color=:gray, linestyle=:dash,
          linewidth=1, label="σ² = 1/d (random)")

    # ── Panel C: β*_predicted (Gaussian model) vs β*_actual ──────────────
    p3 = scatter(β_pred_gauss, β_actual,
                 xlabel="β*_predicted = τ*_G(K) / σ", ylabel="β* [empirical]",
                 title="(C) Gaussian Prediction", titlefontsize=9,
                 legend=:topleft, legendfontsize=7,
                 markersize=5, markercolor=RGB(0.839, 0.153, 0.157),
                 label="Data", grid=true)
    lo3 = min(minimum(β_pred_gauss), minimum(β_actual)) * 0.8
    hi3 = max(maximum(β_pred_gauss), maximum(β_actual)) * 1.1
    plot!(p3, [lo3, hi3], [lo3, hi3], color=:gray, linestyle=:dash,
          linewidth=1, label="y = x")
    R2_g = 1 - sum((β_actual .- β_pred_gauss).^2) / sum((β_actual .- mean(β_actual)).^2)
    annotate!(p3, [(lo3 + 0.1*(hi3-lo3), hi3 - 0.1*(hi3-lo3),
              text(@sprintf("R² = %.3f", R2_g), 9, :left))])

    # ── Panel D: β*_predicted (empirical τ*) vs β*_actual ────────────────
    p4 = scatter(β_pred_emp, β_actual,
                 xlabel="β*_predicted = τ*_emp / σ", ylabel="β* [empirical]",
                 title="(D) Empirical τ* Prediction", titlefontsize=9,
                 legend=:topleft, legendfontsize=7,
                 markersize=5, markercolor=RGB(1.0, 0.498, 0.055),
                 label="Data", grid=true)
    lo4 = min(minimum(β_pred_emp), minimum(β_actual)) * 0.8
    hi4 = max(maximum(β_pred_emp), maximum(β_actual)) * 1.1
    plot!(p4, [lo4, hi4], [lo4, hi4], color=:gray, linestyle=:dash,
          linewidth=1, label="y = x")
    R2_e = 1 - sum((β_actual .- β_pred_emp).^2) / sum((β_actual .- mean(β_actual)).^2)
    annotate!(p4, [(lo4 + 0.1*(hi4-lo4), hi4 - 0.1*(hi4-lo4),
              text(@sprintf("R² = %.3f", R2_e), 9, :left))])

    # ── Panel E: Excess kurtosis across families ─────────────────────────
    fam_names = [r.Family for r in protein_results]
    kurt_vals = [r.kurtosis for r in protein_results]
    p5 = bar(kurt_vals,
             xticks=(1:length(fam_names), fam_names), xrotation=45,
             ylabel="Excess Kurtosis",
             title="(E) Similarity Kurtosis", titlefontsize=9,
             legend=false, color=RGB(0.58, 0.404, 0.741),
             grid=:y)
    hline!(p5, [0.0], color=:gray, linestyle=:dash, linewidth=1, label=false)

    # ── Panel F: σ²·d ratio (deviation from random) ─────────────────────
    σ²d = [r.σ²_times_d for r in protein_results]
    p6 = bar(σ²d,
             xticks=(1:length(fam_names), fam_names), xrotation=45,
             ylabel="σ² × d",
             title="(F) Variance Ratio σ²d (=1 if random)", titlefontsize=9,
             legend=false, color=RGB(0.122, 0.467, 0.706),
             grid=:y)
    hline!(p6, [1.0], color=:gray, linestyle=:dash, linewidth=1, label=false)

    composite = plot(p1, p2, p3, p4, p5, p6,
                     layout=(2, 3), size=(1500, 900),
                     left_margin=8Plots.mm, bottom_margin=12Plots.mm)

    savefig(composite, joinpath(FIGS_DIR, "analytical_beta_star.pdf"))
    @info "Saved: $(joinpath(FIGS_DIR, "analytical_beta_star.pdf"))"

    paper_figs = joinpath(_EXPERIMENT_ROOT, "..", "..", "arxiv", "figs")
    if isdir(paper_figs)
        cp(joinpath(FIGS_DIR, "analytical_beta_star.pdf"),
           joinpath(paper_figs, "analytical_beta_star.pdf"), force=true)
        @info "Copied to arxiv/figs/"
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# PART 4: Stored-pattern probes — matching what find_entropy_inflection does
#
# The empirical β* uses STORED PATTERNS m_k as probes (not random ξ ∈ S^{d-1}).
# When ξ = m_k, the similarity vector has a special structure:
#   s_k = m_k^T m_k = 1       (self-similarity, since ||m_k|| = 1)
#   s_j = m_j^T m_k            (cross-similarity, j ≠ k)
#
# The self-similarity s_k = 1 creates a gap relative to the cross-similarities
# (which are ~O(1/√d)). The entropy inflection occurs when the self-term starts
# to dominate the softmax, requiring less β than for random probes.
# ══════════════════════════════════════════════════════════════════════════════

"""
Compute similarity distribution moments using STORED PATTERNS as probes,
matching the behavior of find_entropy_inflection.
Returns both the full distribution stats and the cross-only (j≠k) stats.
"""
function stored_pattern_moments(X̂::Matrix{Float64})
    d, K = size(X̂)

    # Gram matrix G[j,k] = m_j^T m_k
    G = X̂' * X̂  # K × K

    # Collect cross-similarity statistics (j ≠ k) across all probe patterns
    n_probes = min(K, 20)  # match find_entropy_inflection
    cross_sims = Float64[]
    all_gaps = Float64[]       # gap = 1 - max_{j≠k} s_j
    all_σ²_cross = Float64[]

    for k in 1:n_probes
        sims = G[:, k]         # similarities for probe m_k
        self_sim = sims[k]     # should be ≈ 1.0

        # Cross-similarities
        cross = [sims[j] for j in 1:K if j != k]
        append!(cross_sims, cross)

        push!(all_gaps, self_sim - maximum(cross))
        push!(all_σ²_cross, var(cross; corrected=false))
    end

    σ²_cross = var(cross_sims; corrected=false)
    μ_cross = mean(cross_sims)
    σ_cross = sqrt(σ²_cross)

    # Full distribution (including self-similarity)
    all_sims = Float64[]
    for k in 1:n_probes
        append!(all_sims, G[:, k])
    end
    σ²_full = var(all_sims; corrected=false)
    σ_full = sqrt(σ²_full)
    μ_full = mean(all_sims)

    # Effective gap statistics
    gap_mean = mean(all_gaps)
    gap_std = std(all_gaps)

    # Excess kurtosis of cross-similarities
    z_cross = (cross_sims .- μ_cross) ./ σ_cross
    kurt_cross = mean(z_cross.^4) - 3.0

    return (
        σ²_cross = σ²_cross, σ_cross = σ_cross, μ_cross = μ_cross,
        σ²_full = σ²_full, σ_full = σ_full, μ_full = μ_full,
        gap_mean = gap_mean, gap_std = gap_std,
        kurt_cross = kurt_cross,
        σ²_cross_per_probe = mean(all_σ²_cross),
    )
end

"""
Find τ* using STORED PATTERNS as probes (matching find_entropy_inflection exactly).
For each stored probe m_k, compute similarities {m_j^T m_k}_{j=1..K},
standardize them, and find the entropy inflection in τ-space.
"""
function stored_tau_star(X̂::Matrix{Float64})
    d, K = size(X̂)
    G = X̂' * X̂
    n_probes = min(K, 20)

    τ_stars = Float64[]
    for k in 1:n_probes
        sims = G[:, k]
        σ_s = std(sims; corrected=false)
        z = sims ./ (σ_s + 1e-12)  # standardize
        τ_star = find_tau_star(z; τ_range=(0.01, 100.0))
        push!(τ_stars, τ_star)
    end

    return mean(τ_stars), std(τ_stars) / sqrt(length(τ_stars))
end

"""
Direct β* prediction from stored-pattern probe:
Instead of standardizing, directly find the β inflection of the raw similarity scores
(this should match find_entropy_inflection closely).
"""
function stored_beta_star_direct(X̂::Matrix{Float64})
    d, K = size(X̂)
    G = X̂' * X̂
    n_probes = min(K, 20)

    β_stars = Float64[]
    for k in 1:n_probes
        sims = G[:, k]  # raw similarities (NOT standardized)
        # Find inflection in β-space directly
        β_star = find_tau_star(sims; n_tau=200, τ_range=(0.1, 500.0))
        push!(β_stars, β_star)
    end

    return mean(β_stars), std(β_stars) / sqrt(length(β_stars))
end

function part4_stored_pattern_analysis(tau_gauss, tau_coeffs, protein_results)
    @info "\n══════════════════════════════════════════════════════════"
    @info "  PART 4: Stored-pattern probe analysis"
    @info "  (Matching what find_entropy_inflection actually computes)"
    @info "══════════════════════════════════════════════════════════\n"

    results4 = Vector{NamedTuple}()

    # Header
    @info @sprintf("  %-12s %4s %4s %8s %8s %8s %8s %8s %8s %8s",
        "Family", "K", "d", "β*_act", "β*_dir", "σ_cross", "σ_rand", "gap", "τ*_st", "τ*_st/σ")
    @info "  " * "-"^100

    for fam in FAMILIES
        sto_path = joinpath(PFAM_DATA_DIR, fam.id, "$(fam.id)_seed.sto")
        if !isfile(sto_path); continue; end

        raw_seqs = parse_stockholm(sto_path)
        char_mat, _ = clean_alignment(raw_seqs)
        X̂, _, _, _ = build_memory_matrix(char_mat; pratio=0.95)
        d, K = size(X̂)

        # Empirical β*
        inflection = find_entropy_inflection(X̂)
        β_star = inflection.β_star

        # Stored-pattern moments
        sp = stored_pattern_moments(X̂)

        # Random-probe moments (for comparison)
        rp = similarity_moments(X̂)

        # Stored-pattern τ*
        τ_stored, τ_stored_se = stored_tau_star(X̂)

        # Direct β* from stored patterns (should closely match find_entropy_inflection)
        β_direct, β_direct_se = stored_beta_star_direct(X̂)

        @info @sprintf("  %-12s %4d %4d %8.3f %8.3f %8.4f %8.4f %8.4f %8.4f %8.3f",
            fam.name, K, d, β_star, β_direct, sp.σ_cross, rp.σ, sp.gap_mean, τ_stored, τ_stored/sp.σ_cross)
        @info @sprintf("    μ_cross=%.5f, σ²_cross=%.6f, σ²_rand=%.6f, ratio=%.3f, kurt_cross=%.3f",
            sp.μ_cross, sp.σ²_cross, rp.σ², sp.σ²_cross/rp.σ², sp.kurt_cross)

        push!(results4, (
            Family = fam.name, K = K, d = d, sqrt_d = sqrt(d),
            β_star = β_star, β_direct = β_direct, β_direct_se = β_direct_se,
            σ_cross = sp.σ_cross, σ²_cross = sp.σ²_cross,
            σ_rand = rp.σ, σ²_rand = rp.σ²,
            μ_cross = sp.μ_cross,
            σ_full = sp.σ_full, σ²_full = sp.σ²_full,
            gap_mean = sp.gap_mean, gap_std = sp.gap_std,
            kurt_cross = sp.kurt_cross,
            τ_stored = τ_stored, τ_stored_se = τ_stored_se,
            σ²_ratio = sp.σ²_cross / rp.σ²,
        ))
    end

    # WW scaling subsamples
    @info "\n  WW scaling subsamples:"
    ww_sto = joinpath(PFAM_DATA_DIR, "PF00397", "PF00397_seed.sto")
    if isfile(ww_sto)
        raw_seqs = parse_stockholm(ww_sto)
        full_char_mat, _ = clean_alignment(raw_seqs)
        K_full = size(full_char_mat, 1)

        for K_sub in [20, 50, 100, 200, 400]
            rng = Random.MersenneTwister(1000 * K_sub + 1)
            K_actual = min(K_sub, K_full)
            indices = sort(Random.randperm(rng, K_full)[1:K_actual])
            char_sub = full_char_mat[indices, :]

            X̂, _, _, _ = build_memory_matrix(char_sub; pratio=0.95)
            d, K = size(X̂)
            inflection = find_entropy_inflection(X̂)
            β_star = inflection.β_star
            sp = stored_pattern_moments(X̂)
            rp = similarity_moments(X̂)
            τ_stored, τ_stored_se = stored_tau_star(X̂)
            β_direct, β_direct_se = stored_beta_star_direct(X̂)

            @info @sprintf("  WW K=%d: d=%d, β*=%.3f, β*_dir=%.3f, σ_cross=%.4f, σ_rand=%.4f, gap=%.4f, ratio=%.3f",
                K_sub, d, β_star, β_direct, sp.σ_cross, rp.σ, sp.gap_mean, sp.σ²_cross/rp.σ²)

            push!(results4, (
                Family = "WW_K$(K_sub)", K = K, d = d, sqrt_d = sqrt(d),
                β_star = β_star, β_direct = β_direct, β_direct_se = β_direct_se,
                σ_cross = sp.σ_cross, σ²_cross = sp.σ²_cross,
                σ_rand = rp.σ, σ²_rand = rp.σ²,
                μ_cross = sp.μ_cross,
                σ_full = sp.σ_full, σ²_full = sp.σ²_full,
                gap_mean = sp.gap_mean, gap_std = sp.gap_std,
                kurt_cross = sp.kurt_cross,
                τ_stored = τ_stored, τ_stored_se = τ_stored_se,
                σ²_ratio = sp.σ²_cross / rp.σ²,
            ))
        end
    end

    # ── Analysis ──────────────────────────────────────────────────────────
    β_actual = [r.β_star for r in results4]
    β_direct = [r.β_direct for r in results4]
    σ_cross = [r.σ_cross for r in results4]
    σ_rand = [r.σ_rand for r in results4]
    τ_stored = [r.τ_stored for r in results4]
    gaps = [r.gap_mean for r in results4]
    σ²_ratio = [r.σ²_ratio for r in results4]

    SS_tot = sum((β_actual .- mean(β_actual)).^2)

    # Test 1: β* = τ*_stored / σ_cross
    β_pred1 = τ_stored ./ σ_cross
    R2_1 = 1 - sum((β_actual .- β_pred1).^2) / SS_tot
    r_1 = cor(β_pred1, β_actual)
    @info @sprintf("\n  β* = τ*_stored / σ_cross:  r = %.3f, R² = %.3f", r_1, R2_1)

    # Test 2: β* ≈ β*_direct (should be ≈1 if our computation matches)
    R2_dir = 1 - sum((β_actual .- β_direct).^2) / SS_tot
    r_dir = cor(β_direct, β_actual)
    @info @sprintf("  β* vs β*_direct:           r = %.3f, R² = %.3f", r_dir, R2_dir)

    # Test 3: Regression β* = a + b/σ_cross
    X_reg = hcat(ones(length(σ_cross)), 1.0 ./ σ_cross)
    c_reg = X_reg \ β_actual
    y_reg = X_reg * c_reg
    R2_reg = 1 - sum((β_actual .- y_reg).^2) / SS_tot
    @info @sprintf("  β* = %.3f + %.4f/σ_cross:  R² = %.3f", c_reg[1], c_reg[2], R2_reg)

    # Test 4: σ²_cross vs σ²_rand — how much bigger?
    @info @sprintf("\n  σ²_cross / σ²_rand:  mean = %.3f, std = %.3f",
        mean(σ²_ratio), std(σ²_ratio))

    # Test 5: Does the gap predict β*?
    r_gap = cor(1.0 ./ gaps, β_actual)
    @info @sprintf("  cor(1/gap, β*) = %.3f", r_gap)

    # ── Generate Part 4 figures ───────────────────────────────────────────
    generate_part4_figures(results4, β_actual, β_direct, β_pred1, σ_cross, σ_rand, gaps)

    # Save
    df4 = DataFrame(results4)
    CSV.write(joinpath(OUTPUT_DIR, "stored_pattern_analysis.csv"), df4)
    @info "Saved: $(joinpath(OUTPUT_DIR, "stored_pattern_analysis.csv"))"

    return results4
end

function generate_part4_figures(results4, β_actual, β_direct, β_pred_stored, σ_cross, σ_rand, gaps)
    mkpath(FIGS_DIR)

    d_vals = [r.d for r in results4]
    fam_names = [r.Family for r in results4]

    # ── Panel A: β*_direct vs β*_actual (validation that our computation matches)
    p1 = scatter(β_direct, β_actual,
                 xlabel="β*_direct (stored-pattern inflection)",
                 ylabel="β* [find_entropy_inflection]",
                 title="(A) Direct Computation Validation", titlefontsize=9,
                 legend=false, markersize=5,
                 markercolor=RGB(0.122, 0.467, 0.706), grid=true)
    lo = min(minimum(β_direct), minimum(β_actual)) * 0.8
    hi = max(maximum(β_direct), maximum(β_actual)) * 1.1
    plot!(p1, [lo, hi], [lo, hi], color=:gray, linestyle=:dash, linewidth=1)
    R2 = 1 - sum((β_actual .- β_direct).^2) / sum((β_actual .- mean(β_actual)).^2)
    annotate!(p1, [(lo + 0.1*(hi-lo), hi - 0.1*(hi-lo),
              text(@sprintf("R² = %.3f", R2), 9, :left))])

    # ── Panel B: σ²_cross vs σ²_rand — how much bigger?
    σ²_cross = [r.σ²_cross for r in results4]
    σ²_rand = [r.σ²_rand for r in results4]
    p2 = scatter(σ²_rand, σ²_cross,
                 xlabel="σ²_random (1/d)", ylabel="σ²_cross (stored probes)",
                 title="(B) Cross vs Random Similarity Variance", titlefontsize=9,
                 legend=false, markersize=5,
                 markercolor=RGB(0.173, 0.627, 0.173), grid=true)
    lo2 = min(minimum(σ²_rand), minimum(σ²_cross)) * 0.5
    hi2 = max(maximum(σ²_rand), maximum(σ²_cross)) * 1.3
    plot!(p2, [lo2, hi2], [lo2, hi2], color=:gray, linestyle=:dash, linewidth=1, label="y=x")

    # ── Panel C: σ²_cross/σ²_rand ratio vs d
    ratios = σ²_cross ./ σ²_rand
    p3 = scatter(d_vals, ratios,
                 xlabel="d_eff", ylabel="σ²_cross / σ²_random",
                 title="(C) Variance Ratio", titlefontsize=9,
                 legend=false, markersize=5,
                 markercolor=RGB(0.839, 0.153, 0.157), grid=true)
    hline!(p3, [1.0], color=:gray, linestyle=:dash, linewidth=1)

    # ── Panel D: β* vs τ*_stored/σ_cross (main test)
    p4 = scatter(β_pred_stored, β_actual,
                 xlabel="τ*_stored / σ_cross", ylabel="β* [empirical]",
                 title="(D) Stored-Pattern Prediction", titlefontsize=9,
                 legend=false, markersize=5,
                 markercolor=RGB(1.0, 0.498, 0.055), grid=true)
    lo4 = min(minimum(β_pred_stored), minimum(β_actual)) * 0.8
    hi4 = max(maximum(β_pred_stored), maximum(β_actual)) * 1.1
    plot!(p4, [lo4, hi4], [lo4, hi4], color=:gray, linestyle=:dash, linewidth=1)
    R2_4 = 1 - sum((β_actual .- β_pred_stored).^2) / sum((β_actual .- mean(β_actual)).^2)
    annotate!(p4, [(lo4 + 0.1*(hi4-lo4), hi4 - 0.1*(hi4-lo4),
              text(@sprintf("R² = %.3f", R2_4), 9, :left))])

    # ── Panel E: Self-similarity gap vs d
    p5 = scatter(d_vals, gaps,
                 xlabel="d_eff", ylabel="Gap = 1 - max(m_j^T m_k)",
                 title="(E) Self-Similarity Gap", titlefontsize=9,
                 legend=false, markersize=5,
                 markercolor=RGB(0.58, 0.404, 0.741), grid=true)

    # ── Panel F: β* vs 1/gap
    inv_gap = 1.0 ./ gaps
    p6 = scatter(inv_gap, β_actual,
                 xlabel="1 / gap", ylabel="β* [empirical]",
                 title="(F) β* vs 1/Gap", titlefontsize=9,
                 legend=false, markersize=5,
                 markercolor=RGB(0.122, 0.467, 0.706), grid=true)
    r_gp = cor(inv_gap, β_actual)
    X_g = hcat(ones(length(inv_gap)), inv_gap)
    c_g = X_g \ β_actual
    x_gl = range(minimum(inv_gap)*0.9, maximum(inv_gap)*1.1, length=100)
    y_gl = c_g[1] .+ c_g[2] .* x_gl
    plot!(p6, x_gl, y_gl, color=:black, linewidth=1.5)
    R2_g = 1 - sum((β_actual .- (X_g * c_g)).^2) / sum((β_actual .- mean(β_actual)).^2)
    annotate!(p6, [(minimum(inv_gap)*1.1, maximum(β_actual)*0.95,
              text(@sprintf("r = %.3f\nR² = %.3f", r_gp, R2_g), 8, :left))])

    composite = plot(p1, p2, p3, p4, p5, p6,
                     layout=(2, 3), size=(1500, 900),
                     left_margin=8Plots.mm, bottom_margin=12Plots.mm)

    savefig(composite, joinpath(FIGS_DIR, "stored_pattern_analysis.pdf"))
    @info "Saved: $(joinpath(FIGS_DIR, "stored_pattern_analysis.pdf"))"

    paper_figs = joinpath(_EXPERIMENT_ROOT, "..", "..", "arxiv", "figs")
    if isdir(paper_figs)
        cp(joinpath(FIGS_DIR, "stored_pattern_analysis.pdf"),
           joinpath(paper_figs, "stored_pattern_analysis.pdf"), force=true)
        @info "Copied to arxiv/figs/"
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# Run everything
# ══════════════════════════════════════════════════════════════════════════════

tau_gauss, tau_coeffs = part1_gaussian_inflection()
protein_results = part2_protein_moments()
part3_combined(tau_gauss, tau_coeffs, protein_results)
part4_results = part4_stored_pattern_analysis(tau_gauss, tau_coeffs, protein_results)

@info "\nDone."
