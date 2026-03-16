#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Leave-One-Family-Out Cross-Validation for beta-star Prediction
#
# Shows that the beta-star ~ sqrt(d) relationship generalizes out-of-sample.
# Two LOOCV schemes:
#   1. Leave-one-family-out on the full 33-point dataset (8 families + 25 WW)
#   2. Leave-one-out on the 8-family-only dataset
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment ..."
const _EXPERIMENT_ROOT = @__DIR__
include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))

const BETA_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "data")
const BETA_FIG_DIR  = joinpath(_EXPERIMENT_ROOT, "figs")
mkpath(BETA_DATA_DIR)
mkpath(BETA_FIG_DIR)

# ══════════════════════════════════════════════════════════════════════════════
# Load the 33-point dataset
# ══════════════════════════════════════════════════════════════════════════════

assert_inputs_exist(
    [joinpath(BETA_DATA_DIR, "beta_prediction_data.csv")];
    context="run_beta_loocv — run run_beta_prediction.jl first"
)
df = CSV.read(joinpath(BETA_DATA_DIR, "beta_prediction_data.csv"), DataFrame)
@info "Loaded $(nrow(df)) data points"

# Extract the 8 unique family names (treating each WW scaling replicate as "WW")
# For leave-one-family-out, we group by the base family name
function base_family(row)
    fam = row.Family
    # WW scaling replicates have names like "WW (K=20, r=1)"
    startswith(fam, "WW") ? "WW" : fam
end

df.BaseFam = [base_family(row) for row in eachrow(df)]
unique_families = unique(df.BaseFam)
@info "Unique families: $(unique_families)"

# ══════════════════════════════════════════════════════════════════════════════
# Scheme 1: Leave-one-family-out on full 33-point dataset
# ══════════════════════════════════════════════════════════════════════════════

@info "Scheme 1: Leave-one-family-out (full dataset, n=33)"

loocv_rows = Vector{NamedTuple}()
predictions_full = fill(NaN, nrow(df))

for held_out_fam in unique_families
    train_mask = df.BaseFam .!= held_out_fam
    test_mask  = df.BaseFam .== held_out_fam

    # fit OLS: beta_star = intercept + slope * sqrt_d
    y_train = df.β_star[train_mask]
    X_train = hcat(ones(sum(train_mask)), df.sqrt_d[train_mask])
    theta = X_train \ y_train

    # predict on held-out
    y_test = df.β_star[test_mask]
    X_test = hcat(ones(sum(test_mask)), df.sqrt_d[test_mask])
    y_pred = X_test * theta

    predictions_full[test_mask] .= y_pred

    # per-family error
    rmse_fam = sqrt(mean((y_test .- y_pred).^2))
    mae_fam  = mean(abs.(y_test .- y_pred))

    push!(loocv_rows, (
        Family    = held_out_fam,
        n_held    = sum(test_mask),
        n_train   = sum(train_mask),
        intercept = theta[1],
        slope     = theta[2],
        RMSE      = rmse_fam,
        MAE       = mae_fam,
        mean_obs  = mean(y_test),
        mean_pred = mean(y_pred),
    ))

    @info "  Held out $held_out_fam ($(sum(test_mask)) pts): " *
          "RMSE=$(round(rmse_fam, digits=4)), " *
          "coefs=[$(round(theta[1], digits=3)), $(round(theta[2], digits=3))]"
end

# overall LOOCV metrics
valid = .!isnan.(predictions_full)
ss_res = sum((df.β_star[valid] .- predictions_full[valid]).^2)
ss_tot = sum((df.β_star[valid] .- mean(df.β_star[valid])).^2)
loocv_r2_full   = 1.0 - ss_res / ss_tot
loocv_rmse_full = sqrt(mean((df.β_star[valid] .- predictions_full[valid]).^2))

@info "Full-dataset LOOCV: R^2 = $(round(loocv_r2_full, digits=4)), RMSE = $(round(loocv_rmse_full, digits=4))"

loocv_df = DataFrame(loocv_rows)
CSV.write(joinpath(BETA_DATA_DIR, "beta_loocv_full.csv"), loocv_df)

# ══════════════════════════════════════════════════════════════════════════════
# Scheme 2: 8-family-only LOO (leave one of 8 points out, fit on 7)
# ══════════════════════════════════════════════════════════════════════════════

@info "Scheme 2: 8-family LOO (one point per family)"

pfam_only = filter(r -> r.Source == "Pfam", df)
n8 = nrow(pfam_only)
predictions_8 = fill(NaN, n8)

loo8_rows = Vector{NamedTuple}()

for i in 1:n8
    train_idx = setdiff(1:n8, i)
    y_train = pfam_only.β_star[train_idx]
    X_train = hcat(ones(n8 - 1), pfam_only.sqrt_d[train_idx])
    theta = X_train \ y_train

    y_pred = theta[1] + theta[2] * pfam_only.sqrt_d[i]
    predictions_8[i] = y_pred

    push!(loo8_rows, (
        Family    = pfam_only.Family[i],
        observed  = pfam_only.β_star[i],
        predicted = y_pred,
        error     = pfam_only.β_star[i] - y_pred,
        intercept = theta[1],
        slope     = theta[2],
    ))
end

ss_res_8 = sum((pfam_only.β_star .- predictions_8).^2)
ss_tot_8 = sum((pfam_only.β_star .- mean(pfam_only.β_star)).^2)
loocv_r2_8   = 1.0 - ss_res_8 / ss_tot_8
loocv_rmse_8 = sqrt(mean((pfam_only.β_star .- predictions_8).^2))

@info "8-family LOO: R^2 = $(round(loocv_r2_8, digits=4)), RMSE = $(round(loocv_rmse_8, digits=4))"

loo8_df = DataFrame(loo8_rows)
CSV.write(joinpath(BETA_DATA_DIR, "beta_loocv_8family.csv"), loo8_df)

# ══════════════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════════════

summary_df = DataFrame(
    Scheme = ["Full 33-pt LOFO", "8-family LOO"],
    R2     = [loocv_r2_full, loocv_r2_8],
    RMSE   = [loocv_rmse_full, loocv_rmse_8],
    n      = [nrow(df), n8],
)
CSV.write(joinpath(BETA_DATA_DIR, "beta_loocv_summary.csv"), summary_df)

println("\n======================================================")
println("BETA-STAR LOOCV SUMMARY")
println("======================================================")
for r in eachrow(summary_df)
    println("  $(rpad(r.Scheme, 25))  R^2=$(round(r.R2, digits=4))  RMSE=$(round(r.RMSE, digits=4))  (n=$(r.n))")
end
println()

# per-family errors for the 8-family scheme
println("8-family LOO detail:")
for r in eachrow(loo8_df)
    println("  $(rpad(r.Family, 15))  obs=$(round(r.observed, digits=2))  pred=$(round(r.predicted, digits=2))  err=$(round(r.error, digits=3))")
end

@info "Done. Results saved to beta_loocv_*.csv"
