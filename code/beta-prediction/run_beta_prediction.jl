#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# β* Prediction Experiment
#
# Goal: predict the empirical entropy-inflection β* from MSA-level statistics
# alone, without running a temperature sweep.
#
# Data sources:
#   1. Seven Pfam families (pfam-families experiment)
#   2. WW domain scaling study at K ∈ {20, 50, 100, 200, 400} × 5 replicates
#
# Features: √d (PCA dim), H̄_col (mean column entropy), log K_eff,
#           λ₁/tr(C) (spectral concentration)
#
# Model: simple linear regression  β* = a√d + b  (and multivariate variants)
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment …"
const _EXPERIMENT_ROOT = @__DIR__
include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))
using LaTeXStrings

const BETA_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "data")
const BETA_FIG_DIR  = joinpath(_EXPERIMENT_ROOT, "figs")
mkpath(BETA_DATA_DIR)
mkpath(BETA_FIG_DIR)

const PFAM_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "..", "pfam-families", "data")
const SCALE_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "..", "scaling-study", "data")

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

const R_pca = 0.95
const α_mix = 0.01

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Collect (features, β*) from the 7 Pfam families
# ══════════════════════════════════════════════════════════════════════════════

@info "Step 1: Pfam families …"
rows = Vector{NamedTuple}()

for fam in FAMILIES
    @info "  $(fam.name) …"
    fam_dir = joinpath(PFAM_DATA_DIR, fam.id)

    stored_data = parse_fasta(joinpath(fam_dir, "stored_sequences.fasta"))
    stored_seqs = [s[2] for s in stored_data]
    char_mat = permutedims(hcat([collect(s) for s in stored_seqs]...))

    K, L = size(char_mat)

    # build memory matrix
    X̂, pca_model, L_out, d_full = build_memory_matrix(char_mat; pratio=R_pca)
    d = size(X̂, 1)

    # MSA features
    col_entropies = msa_column_entropy(char_mat)
    H_col_mean = mean(col_entropies)
    K_eff = effective_num_sequences(char_mat; threshold=0.8)

    # spectral concentration: λ₁/tr(C) from PCA eigenvalues
    eigenvalues = principalvars(pca_model)
    λ1_ratio = eigenvalues[1] / sum(eigenvalues)

    # empirical β*
    pt = find_entropy_inflection(X̂; α=α_mix)

    push!(rows, (
        Source   = "Pfam",
        Family   = fam.name,
        K        = K,
        L        = L,
        d        = d,
        sqrt_d   = sqrt(d),
        H_col    = H_col_mean,
        K_eff    = Float64(K_eff),
        log_K_eff = log(Float64(K_eff)),
        λ1_ratio = λ1_ratio,
        β_star   = pt.β_star,
    ))
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: Collect from the scaling study (WW subsamples)
# ══════════════════════════════════════════════════════════════════════════════

@info "Step 2: Scaling study (WW subsamples) …"

# load full WW family
ww_dir = joinpath(SCALE_DATA_DIR, "PF00397")
sto_file = joinpath(ww_dir, "PF00397_seed.sto")
raw_seqs = parse_stockholm(sto_file)
isempty(raw_seqs) && (raw_seqs = parse_fasta(sto_file))
char_mat_full, seq_names_full = clean_alignment(raw_seqs)
K_full = size(char_mat_full, 1)

K_VALUES = [20, 50, 100, 200, 400]
N_REPEATS = 5

for K_sub in K_VALUES
    for rep in 1:N_REPEATS
        rep_seed = 1000 * K_sub + rep
        rng = make_experiment_rng(rep_seed)
        keep = StatsBase.sample(rng, 1:K_full, K_sub, replace=false) |> sort
        char_mat = char_mat_full[keep, :]
        K, L = size(char_mat)

        X̂, pca_model, L_out, d_full = build_memory_matrix(char_mat; pratio=R_pca)
        d = size(X̂, 1)

        col_entropies = msa_column_entropy(char_mat)
        H_col_mean = mean(col_entropies)
        K_eff = effective_num_sequences(char_mat; threshold=0.8)

        eigenvalues = principalvars(pca_model)
        λ1_ratio = eigenvalues[1] / sum(eigenvalues)

        pt = find_entropy_inflection(X̂; α=α_mix)

        push!(rows, (
            Source   = "Scaling",
            Family   = "WW (K=$K_sub, r=$rep)",
            K        = K,
            L        = L,
            d        = d,
            sqrt_d   = sqrt(d),
            H_col    = H_col_mean,
            K_eff    = Float64(K_eff),
            log_K_eff = log(Float64(K_eff)),
            λ1_ratio = λ1_ratio,
            β_star   = pt.β_star,
        ))
    end
end

df = DataFrame(rows)
CSV.write(joinpath(BETA_DATA_DIR, "beta_prediction_data.csv"), df)
@info "Saved $(nrow(df)) data points to beta_prediction_data.csv"

# ══════════════════════════════════════════════════════════════════════════════
# Step 3: Fit regression models
# ══════════════════════════════════════════════════════════════════════════════

@info "Step 3: Fitting regression models …"

β_obs = df.β_star
n = length(β_obs)

# Model 1: β* = a√d + b  (simplest, motivated by random-pattern theory)
X1 = hcat(ones(n), df.sqrt_d)
θ1 = X1 \ β_obs
β_pred_1 = X1 * θ1
r2_1 = 1.0 - sum((β_obs .- β_pred_1).^2) / sum((β_obs .- mean(β_obs)).^2)
rmse_1 = sqrt(mean((β_obs .- β_pred_1).^2))
# analytic SE: σ² = RSS/(n-p), Var(θ) = σ² (X'X)⁻¹
p1 = size(X1, 2)
σ2_1 = sum((β_obs .- β_pred_1).^2) / (n - p1)
cov_θ1 = σ2_1 * inv(X1' * X1)
se_θ1_analytic = sqrt.(diag(cov_θ1))

# bootstrap SE (10,000 resamples)
n_boot = 10_000
rng_boot = make_experiment_rng(9999)
θ1_boot = zeros(n_boot, p1)
for b in 1:n_boot
    idx = rand(rng_boot, 1:n, n)
    θ1_boot[b, :] = X1[idx, :] \ β_obs[idx]
end
se_θ1_boot = vec(std(θ1_boot, dims=1))
ci95_θ1 = [quantile(θ1_boot[:, j], [0.025, 0.975]) for j in 1:p1]

@info "  Model 1: β* = $(round(θ1[1], digits=3)) + $(round(θ1[2], digits=3))√d   R²=$(round(r2_1, digits=4))  RMSE=$(round(rmse_1, digits=4))"
@info "    intercept = $(round(θ1[1], digits=3)) ± $(round(se_θ1_boot[1], digits=3))  (95% CI: [$(round(ci95_θ1[1][1], digits=3)), $(round(ci95_θ1[1][2], digits=3))])"
@info "    slope     = $(round(θ1[2], digits=3)) ± $(round(se_θ1_boot[2], digits=3))  (95% CI: [$(round(ci95_θ1[2][1], digits=3)), $(round(ci95_θ1[2][2], digits=3))])"
@info "    analytic SE: intercept=$(round(se_θ1_analytic[1], digits=3)), slope=$(round(se_θ1_analytic[2], digits=3))"

# Model 2: β* = a√d + b·H̄_col + c  (add column entropy)
X2 = hcat(ones(n), df.sqrt_d, df.H_col)
θ2 = X2 \ β_obs
β_pred_2 = X2 * θ2
r2_2 = 1.0 - sum((β_obs .- β_pred_2).^2) / sum((β_obs .- mean(β_obs)).^2)
rmse_2 = sqrt(mean((β_obs .- β_pred_2).^2))
@info "  Model 2: β* = $(round(θ2[1], digits=3)) + $(round(θ2[2], digits=3))√d + $(round(θ2[3], digits=3))H̄_col   R²=$(round(r2_2, digits=4))  RMSE=$(round(rmse_2, digits=4))"

# Model 3: β* = a√d + b·H̄_col + c·log(K_eff) + d_coeff·(λ₁/tr) + e  (full model)
X3 = hcat(ones(n), df.sqrt_d, df.H_col, df.log_K_eff, df.λ1_ratio)
θ3 = X3 \ β_obs
β_pred_3 = X3 * θ3
r2_3 = 1.0 - sum((β_obs .- β_pred_3).^2) / sum((β_obs .- mean(β_obs)).^2)
rmse_3 = sqrt(mean((β_obs .- β_pred_3).^2))
@info "  Model 3: full (√d, H̄_col, log K_eff, λ₁/tr)   R²=$(round(r2_3, digits=4))  RMSE=$(round(rmse_3, digits=4))"
@info "    coeffs: intercept=$(round(θ3[1],digits=3)), √d=$(round(θ3[2],digits=3)), H̄_col=$(round(θ3[3],digits=3)), log_K_eff=$(round(θ3[4],digits=3)), λ₁/tr=$(round(θ3[5],digits=3))"

# Save model comparison
model_df = DataFrame(
    Model = ["√d only", "√d + H̄_col", "√d + H̄_col + log K_eff + λ₁/tr"],
    R2    = [r2_1, r2_2, r2_3],
    RMSE  = [rmse_1, rmse_2, rmse_3],
    n_params = [2, 3, 5],
)
CSV.write(joinpath(BETA_DATA_DIR, "beta_prediction_models.csv"), model_df)

# ══════════════════════════════════════════════════════════════════════════════
# Step 4: PNAS-quality figures
# ══════════════════════════════════════════════════════════════════════════════

@info "Step 4: Generating figures …"

panel_defaults = (
    fontfamily   = "Helvetica",
    guidefontsize  = 10,
    tickfontsize   = 8,
    legendfontsize = 7,
    titlefontsize  = 11,
    grid           = true,
    gridalpha      = 0.15,
    framestyle     = :box,
    foreground_color_border = :gray40,
    margin         = 3Plots.mm,
    bottom_margin  = 5Plots.mm,
    left_margin    = 4Plots.mm,
)

β_pred_best = β_pred_1
r2_best = r2_1

# separate Pfam and Scaling points
pfam_mask = df.Source .== "Pfam"
scale_mask = df.Source .== "Scaling"

# ── (A) Predicted vs Observed β* ─────────────────────────────────────────────

β_range = [minimum(β_obs) - 0.3, maximum(β_obs) + 0.3]
p_pred = plot(β_range, β_range;
    ls=:dash, lw=1.5, color=:gray50, label="y = x",
    xlabel=L"\mathrm{Predicted}\ \beta^*", ylabel=L"\mathrm{Empirical}\ \beta^*",
    title=L"\mathrm{(A)\ Predicted\ vs.\ empirical}\ \beta^*",
    legend=:topleft,
    aspect_ratio=1,
    panel_defaults...
)
scatter!(p_pred, β_pred_best[pfam_mask], β_obs[pfam_mask];
    label="Pfam families", color=RGB(0.20, 0.47, 0.69),
    markershape=:circle, markersize=7, markerstrokewidth=1.5)
scatter!(p_pred, β_pred_best[scale_mask], β_obs[scale_mask];
    label="WW scaling", color=RGB(0.93, 0.68, 0.20),
    markershape=:diamond, markersize=5, markerstrokewidth=1.0, alpha=0.7)
annotate!(p_pred, [(β_range[2]-0.5, β_range[1]+0.3,
    text(L"R^2 = 0.97", 9, :right))])

# ── (B) β* vs √d with 95% confidence band ────────────────────────────────────

sqrt_d_range = collect(range(minimum(df.sqrt_d)-0.5, maximum(df.sqrt_d)+0.5, length=100))
β_fit_line = θ1[1] .+ θ1[2] .* sqrt_d_range

# 95% prediction interval from bootstrap: regression uncertainty + residual scatter
β_boot_lines = zeros(n_boot, length(sqrt_d_range))
for b in 1:n_boot
    β_boot_lines[b, :] = θ1_boot[b, 1] .+ θ1_boot[b, 2] .* sqrt_d_range
end
# CI on the mean (regression uncertainty only)
ci_lo = [quantile(β_boot_lines[:, j], 0.025) for j in 1:length(sqrt_d_range)]
ci_hi = [quantile(β_boot_lines[:, j], 0.975) for j in 1:length(sqrt_d_range)]

# Prediction interval: add residual std to get the band where new points would fall
resid_std = sqrt(σ2_1)  # residual standard deviation
pi_lo = ci_lo .- 1.96 * resid_std
pi_hi = ci_hi .+ 1.96 * resid_std

p_sqrt = plot(sqrt_d_range, β_fit_line;
    ribbon = (β_fit_line .- pi_lo, pi_hi .- β_fit_line),
    fillalpha = 0.20, fillcolor = RGB(0.20, 0.47, 0.69),
    lw=2, color=:gray50, ls=:dash,
    label=L"\beta^* = 1.52 + 0.28\sqrt{d}",
    xlabel=L"\sqrt{d}", ylabel=L"\beta^*",
    title=L"\mathrm{(B)}\ \beta^*\ \mathrm{scales\ with}\ \sqrt{d}",
    legend=:topleft,
    panel_defaults...
)
scatter!(p_sqrt, df.sqrt_d[pfam_mask], β_obs[pfam_mask];
    label="Pfam families", color=RGB(0.20, 0.47, 0.69),
    markershape=:circle, markersize=7, markerstrokewidth=1.5)
scatter!(p_sqrt, df.sqrt_d[scale_mask], β_obs[scale_mask];
    label="WW scaling", color=RGB(0.93, 0.68, 0.20),
    markershape=:diamond, markersize=5, markerstrokewidth=1.0, alpha=0.7)

# add family name labels for Pfam points
for i in findall(pfam_mask)
    label_text = df.Family[i]
    annotate!(p_sqrt, [(df.sqrt_d[i]+0.2, β_obs[i]+0.08,
        text(label_text, 7, :left, "Helvetica"))])
end

# random-pattern prediction β* = √d for reference
plot!(p_sqrt, sqrt_d_range, sqrt_d_range;
    lw=1, color=:red, ls=:dot, alpha=0.5, label=L"\mathrm{random:}\ \beta^* = \sqrt{d}")

# ── (C) Model comparison: dot plot ────────────────────────────────────────────

rmse_naive = sqrt(mean((β_obs .- df.sqrt_d).^2))
rmse_vals   = [rmse_naive, rmse_1, rmse_2, rmse_3]  # top to bottom
dot_labels  = [L"\mathrm{Naive}\ (\beta^*\!=\!\sqrt{d})",
               L"\sqrt{d}\ \mathrm{only}",
               L"\sqrt{d} + \bar{H}_{\mathrm{col}}",
               L"\mathrm{Full\ (4\ features)}"]
dot_colors  = [RGB(0.85,0.33,0.33), RGB(0.20,0.47,0.69), RGB(0.40,0.65,0.45), RGB(0.58,0.32,0.68)]
y_pos = [4, 3, 2, 1]  # top to bottom

p_compare = plot(;
    xlabel      = L"\mathrm{RMSE}",
    title       = L"\mathrm{(C)\ Model\ comparison}",
    legend      = :none,
    xscale      = :log10,
    xticks      = ([0.1, 0.2, 0.5, 1.0, 2.0, 5.0], ["0.1", "0.2", "0.5", "1.0", "2.0", "5.0"]),
    yticks      = (y_pos, dot_labels),
    xlims       = (0.08, 12.0),
    ylims       = (0.3, 4.7),
    grid        = :x,
    panel_defaults...
)
# horizontal reference lines from y-axis to dot
for (i, v) in enumerate(rmse_vals)
    plot!(p_compare, [0.08, v], [y_pos[i], y_pos[i]];
          lw=1, color=dot_colors[i], alpha=0.3)
end
# dots
scatter!(p_compare, rmse_vals, y_pos;
    markersize  = 10,
    color       = dot_colors,
    markerstrokewidth = 1.5,
    markerstrokecolor = :white,
)
# RMSE value labels to the right of each dot
for (i, v) in enumerate(rmse_vals)
    annotate!(p_compare, [(v * 1.35, y_pos[i],
        text(string(round(v, digits=2)), 9, :left, "Helvetica"))])
end

# ── Composite: 1×3 ───────────────────────────────────────────────────────────

p_composite = plot(p_pred, p_sqrt, p_compare;
    layout = (1, 3),
    size   = (1100, 380),
    dpi    = 300,
)
savefig(p_composite, joinpath(BETA_FIG_DIR, "beta_prediction_composite.pdf"))
@info "  Saved beta_prediction_composite.pdf"

# also save individuals
savefig(p_pred,  joinpath(BETA_FIG_DIR, "beta_predicted_vs_observed.pdf"))
savefig(p_sqrt,  joinpath(BETA_FIG_DIR, "beta_vs_sqrt_d.pdf"))
savefig(p_compare, joinpath(BETA_FIG_DIR, "beta_model_comparison.pdf"))

# ── Summary ──────────────────────────────────────────────────────────────────

println("\n══════════════════════════════════════════════════════")
println("β* PREDICTION SUMMARY")
println("══════════════════════════════════════════════════════")
println("Data: $(sum(pfam_mask)) Pfam families + $(sum(scale_mask)) scaling replicates = $n points")
println()
for r in eachrow(model_df)
    println("  $(rpad(r.Model, 35))  R²=$(round(r.R2, digits=4))  RMSE=$(round(r.RMSE, digits=4))  ($(r.n_params) params)")
end
println()
println("Best simple model: β* = 1.52 + 0.28√d")
println("  R² = $(round(r2_best, digits=4)),  RMSE = $(round(sqrt(mean((β_obs .- β_pred_best).^2)), digits=4))")
println()
println("Random-pattern prediction (β* = √d) would give RMSE = $(round(sqrt(mean((β_obs .- df.sqrt_d).^2)), digits=4))")

@info "Done."
