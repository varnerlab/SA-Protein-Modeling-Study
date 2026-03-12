#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Investigate: Does β* = 1/λ₁(C) from the Jacobian bifurcation condition?
#
# The retrieval map T(ξ) = X softmax(β X^T ξ) has Jacobian at the centroid
# (with uniform attention): J = β C, where C is the sample covariance of the
# unit-norm stored patterns. The fixed point becomes unstable when β λ₁(C) = 1,
# giving β_c = 1/λ₁(C).
#
# This script computes λ₁(C) for each family + scaling study subsample and
# compares β_c with the empirical β* from the entropy inflection.
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
const SCALING_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "..", "scaling-study", "data")
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
# Compute covariance spectrum of unit-norm patterns
# ══════════════════════════════════════════════════════════════════════════════

function analyze_spectrum(X̂::Matrix{Float64})
    d, K = size(X̂)

    # Sample covariance of unit-norm columns
    μ = vec(mean(X̂, dims=2))
    C = (X̂ * X̂') / K - μ * μ'

    eigenvals = sort(eigvals(Symmetric(C)), rev=true)
    λ1 = eigenvals[1]
    tr_C = sum(eigenvals)

    # Bifurcation condition
    β_c = 1.0 / λ1

    # Also compute: second moment matrix spectrum (no centering)
    M2 = (X̂ * X̂') / K
    eigenvals_M2 = sort(eigvals(Symmetric(M2)), rev=true)
    λ1_M2 = eigenvals_M2[1]
    β_c_M2 = 1.0 / λ1_M2

    # Participation ratio (effective dimension of spectrum)
    d_pr = tr_C^2 / sum(eigenvals.^2)

    return (
        d = d, K = K,
        λ1 = λ1, tr_C = tr_C, λ1_ratio = λ1/tr_C,
        β_c = β_c,
        λ1_M2 = λ1_M2, β_c_M2 = β_c_M2,
        d_pr = d_pr,
        eigenvals = eigenvals,
        μ_norm = norm(μ),
    )
end

# ══════════════════════════════════════════════════════════════════════════════
# Main analysis
# ══════════════════════════════════════════════════════════════════════════════

function main()
    mkpath(OUTPUT_DIR)
    mkpath(FIGS_DIR)

    results = Vector{NamedTuple}()

    # ── 7 Pfam families ──────────────────────────────────────────────────
    @info "═══════════════════════════════════════════════════════════"
    @info "  Analyzing 7 Pfam families"
    @info "═══════════════════════════════════════════════════════════"

    for fam in FAMILIES
        sto_path = joinpath(PFAM_DATA_DIR, fam.id, "$(fam.id)_seed.sto")
        if !isfile(sto_path)
            @warn "Missing: $sto_path"
            continue
        end

        raw_seqs = parse_stockholm(sto_path)
        char_mat, seq_names = clean_alignment(raw_seqs)
        X̂, pca_model, L, d_full = build_memory_matrix(char_mat; pratio=0.95)

        # Empirical β*
        inflection = find_entropy_inflection(X̂)
        β_star = inflection.β_star

        # Spectrum analysis
        spec = analyze_spectrum(X̂)

        # PCA eigenvalues (for comparison)
        pca_eigenvals = principalvars(pca_model)
        λ1_pca = pca_eigenvals[1]
        tr_pca = sum(pca_eigenvals)

        @info @sprintf("  %s: d=%d, K=%d, β*=%.3f, β_c=1/λ₁=%.3f, β_c_M2=%.3f, λ₁(C)=%.4f, λ₁(M2)=%.4f, |μ|=%.3f",
            fam.name, spec.d, spec.K, β_star, spec.β_c, spec.β_c_M2,
            spec.λ1, spec.λ1_M2, spec.μ_norm)

        push!(results, (
            Source = "Pfam", Family = fam.name,
            K = spec.K, d = spec.d, sqrt_d = sqrt(spec.d),
            β_star = β_star,
            β_c = spec.β_c,           # 1/λ₁(C)  (centered)
            β_c_M2 = spec.β_c_M2,     # 1/λ₁(M2) (uncentered)
            λ1_C = spec.λ1,
            λ1_M2 = spec.λ1_M2,
            tr_C = spec.tr_C,
            λ1_ratio = spec.λ1_ratio,
            d_pr = spec.d_pr,
            μ_norm = spec.μ_norm,
            λ1_pca = λ1_pca,
            tr_pca = tr_pca,
        ))
    end

    # ── WW scaling study subsamples ──────────────────────────────────────
    @info "\n═══════════════════════════════════════════════════════════"
    @info "  Analyzing WW scaling study subsamples"
    @info "═══════════════════════════════════════════════════════════"

    ww_sto = joinpath(PFAM_DATA_DIR, "PF00397", "PF00397_seed.sto")
    if isfile(ww_sto)
        raw_seqs = parse_stockholm(ww_sto)
        full_char_mat, _ = clean_alignment(raw_seqs)
        K_full = size(full_char_mat, 1)

        for K_sub in [20, 50, 100, 200, 400]
            for rep in 1:5
                rng = Random.MersenneTwister(1000 * K_sub + rep)
                K_actual = min(K_sub, K_full)
                indices = sort(Random.randperm(rng, K_full)[1:K_actual])
                char_sub = full_char_mat[indices, :]

                X̂, pca_model, L, d_full = build_memory_matrix(char_sub; pratio=0.95)
                inflection = find_entropy_inflection(X̂)
                β_star = inflection.β_star
                spec = analyze_spectrum(X̂)

                pca_eigenvals = principalvars(pca_model)

                if rep == 1
                    @info @sprintf("  WW K=%d r=%d: d=%d, β*=%.3f, β_c=%.3f, β_c_M2=%.3f, λ₁(C)=%.4f",
                        K_sub, rep, spec.d, β_star, spec.β_c, spec.β_c_M2, spec.λ1)
                end

                push!(results, (
                    Source = "Scaling", Family = "WW (K=$K_sub, r=$rep)",
                    K = spec.K, d = spec.d, sqrt_d = sqrt(spec.d),
                    β_star = β_star,
                    β_c = spec.β_c,
                    β_c_M2 = spec.β_c_M2,
                    λ1_C = spec.λ1,
                    λ1_M2 = spec.λ1_M2,
                    tr_C = spec.tr_C,
                    λ1_ratio = spec.λ1_ratio,
                    d_pr = spec.d_pr,
                    μ_norm = spec.μ_norm,
                    λ1_pca = pca_eigenvals[1],
                    tr_pca = sum(pca_eigenvals),
                ))
            end
        end
    end

    # ── Save results ─────────────────────────────────────────────────────
    df = DataFrame(results)
    CSV.write(joinpath(OUTPUT_DIR, "bifurcation_analysis.csv"), df)
    @info "\nSaved: $(joinpath(OUTPUT_DIR, "bifurcation_analysis.csv"))"

    # ── Print comparison table ───────────────────────────────────────────
    @info "\n═══════════════════════════════════════════════════════════"
    @info "  Comparison: β* (empirical) vs β_c = 1/λ₁ (bifurcation)"
    @info "═══════════════════════════════════════════════════════════"
    @info @sprintf("  %-20s %5s %5s %8s %8s %8s %8s %8s",
        "Family", "K", "d", "β*", "1/λ₁(C)", "1/λ₁(M2)", "β*/β_c", "β*/β_M2")
    @info "  " * "-"^80

    for r in results
        ratio_c = r.β_star / r.β_c
        ratio_m2 = r.β_star / r.β_c_M2
        @info @sprintf("  %-20s %5d %5d %8.3f %8.3f %8.3f %8.3f %8.3f",
            r.Family, r.K, r.d, r.β_star, r.β_c, r.β_c_M2, ratio_c, ratio_m2)
    end

    # ── Plots ────────────────────────────────────────────────────────────
    generate_plots(df)

    return df
end

function generate_plots(df)
    # Plot 1: β* vs 1/λ₁(C) — does the bifurcation condition hold?
    p1 = scatter(df.β_c, df.β_star,
                 xlabel="β_c = 1/λ₁(C) [bifurcation]",
                 ylabel="β* [empirical inflection]",
                 title="(A) Bifurcation vs Empirical β*",
                 titlefontsize=9, legend=false,
                 markersize=5, markercolor=RGB(0.122, 0.467, 0.706),
                 grid=true)
    # Identity line
    lo = min(minimum(df.β_c), minimum(df.β_star)) * 0.8
    hi = max(maximum(df.β_c), maximum(df.β_star)) * 1.1
    plot!(p1, [lo, hi], [lo, hi], color=:gray, linestyle=:dash, linewidth=1, label=false)

    # Correlation
    r_c = cor(df.β_c, df.β_star)
    annotate!(p1, [(lo + 0.1*(hi-lo), hi - 0.1*(hi-lo),
              text(@sprintf("r = %.3f", r_c), 9, :left))])

    # Plot 2: β* vs 1/λ₁(M2) — uncentered second moment
    p2 = scatter(df.β_c_M2, df.β_star,
                 xlabel="1/λ₁(M₂) [uncentered]",
                 ylabel="β* [empirical inflection]",
                 title="(B) Second Moment vs Empirical β*",
                 titlefontsize=9, legend=false,
                 markersize=5, markercolor=RGB(0.839, 0.153, 0.157),
                 grid=true)
    lo2 = min(minimum(df.β_c_M2), minimum(df.β_star)) * 0.8
    hi2 = max(maximum(df.β_c_M2), maximum(df.β_star)) * 1.1
    plot!(p2, [lo2, hi2], [lo2, hi2], color=:gray, linestyle=:dash, linewidth=1, label=false)
    r_m2 = cor(df.β_c_M2, df.β_star)
    annotate!(p2, [(lo2 + 0.1*(hi2-lo2), hi2 - 0.1*(hi2-lo2),
              text(@sprintf("r = %.3f", r_m2), 9, :left))])

    # Plot 3: β* vs √d — the existing empirical relationship
    p3 = scatter(df.sqrt_d, df.β_star,
                 xlabel="√d_eff",
                 ylabel="β* [empirical inflection]",
                 title="(C) β* vs √d_eff",
                 titlefontsize=9, legend=false,
                 markersize=5, markercolor=RGB(0.173, 0.627, 0.173),
                 grid=true)
    r_d = cor(df.sqrt_d, df.β_star)
    # Fit line
    X_reg = hcat(ones(nrow(df)), df.sqrt_d)
    coeffs = X_reg \ df.β_star
    x_line = range(minimum(df.sqrt_d)*0.9, maximum(df.sqrt_d)*1.1, length=100)
    y_line = coeffs[1] .+ coeffs[2] .* x_line
    plot!(p3, x_line, y_line, color=:black, linewidth=1.5, label=false)
    annotate!(p3, [(minimum(df.sqrt_d)+0.3, maximum(df.β_star)-0.2,
              text(@sprintf("β* = %.2f + %.2f√d\nR² = %.3f", coeffs[1], coeffs[2], r_d^2),
                   8, :left))])

    # Plot 4: λ₁(C) vs d — how does the leading eigenvalue scale?
    p4 = scatter(df.d, df.λ1_C,
                 xlabel="d_eff",
                 ylabel="λ₁(C)",
                 title="(D) Leading Eigenvalue vs d_eff",
                 titlefontsize=9, legend=false,
                 markersize=5, markercolor=RGB(1.0, 0.498, 0.055),
                 grid=true)

    # Plot 5: β* / β_c ratio vs d — is there a systematic relationship?
    ratio = df.β_star ./ df.β_c
    p5 = scatter(df.d, ratio,
                 xlabel="d_eff",
                 ylabel="β* / β_c",
                 title="(E) Ratio β*/β_c vs d_eff",
                 titlefontsize=9, legend=false,
                 markersize=5, markercolor=RGB(0.58, 0.404, 0.741),
                 grid=true)
    hline!(p5, [1.0], color=:gray, linestyle=:dash, linewidth=1, label=false)

    # Plot 6: eigenvalue spectrum for Pfam families
    pfam_rows = filter(r -> r.Source == "Pfam", df)
    p6 = plot(title="(F) Eigenvalue Spectra (Pfam)", titlefontsize=9,
              xlabel="Component index", ylabel="λᵢ / tr(C)",
              legend=:topright, legendfontsize=6, grid=true,
              yscale=:log10)

    # Re-compute eigenvalue spectra for plotting
    for fam in FAMILIES
        sto_path = joinpath(PFAM_DATA_DIR, fam.id, "$(fam.id)_seed.sto")
        if !isfile(sto_path); continue; end
        raw_seqs = parse_stockholm(sto_path)
        char_mat, _ = clean_alignment(raw_seqs)
        X̂, _, _, _ = build_memory_matrix(char_mat; pratio=0.95)
        spec = analyze_spectrum(X̂)
        normalized = spec.eigenvals ./ spec.tr_C
        # Only plot positive eigenvalues
        pos = filter(x -> x > 1e-10, normalized)
        plot!(p6, 1:length(pos), pos, label=fam.name, linewidth=1.5)
    end

    composite = plot(p1, p2, p3, p4, p5, p6,
                     layout=(2, 3), size=(1400, 800),
                     left_margin=8Plots.mm, bottom_margin=8Plots.mm)

    savefig(composite, joinpath(FIGS_DIR, "bifurcation_analysis.pdf"))
    @info "Saved: $(joinpath(FIGS_DIR, "bifurcation_analysis.pdf"))"

    # Copy to paper figs
    paper_figs = joinpath(_EXPERIMENT_ROOT, "..", "..", "paper", "figs")
    if isdir(paper_figs)
        cp(joinpath(FIGS_DIR, "bifurcation_analysis.pdf"),
           joinpath(paper_figs, "bifurcation_analysis.pdf"), force=true)
        @info "Copied to paper/figs/"
    end
end

df = main()
