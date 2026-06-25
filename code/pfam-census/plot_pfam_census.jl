#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Plot Pfam Family Size Distributions
#
# Reads the census CSV produced by run_pfam_census.jl and generates a
# publication-quality figure showing seed vs. full alignment size distributions.
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment …"
const _EXPERIMENT_ROOT = @__DIR__
include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))

const CENSUS_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "data")
const FIG_DIR = joinpath(dirname(_EXPERIMENT_ROOT), "..", "arxiv", "figs")
mkpath(FIG_DIR)
@info "Environment loaded."

# ══════════════════════════════════════════════════════════════════════════════
# Load data
# ══════════════════════════════════════════════════════════════════════════════

csv_path = joinpath(CENSUS_DATA_DIR, "pfam_family_sizes.csv")
if !isfile(csv_path)
    error("Census CSV not found at $csv_path. Run run_pfam_census.jl first.")
end

df = CSV.read(csv_path, DataFrame)
@info "Loaded $(nrow(df)) families from census."

# ══════════════════════════════════════════════════════════════════════════════
# Figure: Empirical CDF of family sizes (seed vs. full)
# ══════════════════════════════════════════════════════════════════════════════

function plot_family_size_ecdf(df::DataFrame)
    seed = sort(df.n_seed)
    full = sort(skipmissing(df.n_full) |> collect)

    # compute ECDFs
    seed_ecdf = (1:length(seed)) ./ length(seed)
    full_ecdf = (1:length(full)) ./ length(full)

    # medians
    med_seed = median(seed)
    med_full = median(full)

    p = plot(
        seed, seed_ecdf,
        label="Seed alignments (n = $(length(seed)))",
        xlabel="Number of sequences per family",
        ylabel="Cumulative fraction of families",
        xscale=:log10,
        linewidth=2,
        color=:steelblue,
        legend=:bottomright,
        size=(700, 450),
        dpi=300,
        left_margin=5Plots.mm,
        bottom_margin=5Plots.mm,
        framestyle=:box,
        grid=true,
        gridalpha=0.2,
    )

    plot!(p, full, full_ecdf,
        label="Full alignments (n = $(length(full)))",
        linewidth=2,
        color=:coral,
    )

    # median lines
    vline!(p, [med_seed], color=:steelblue, linestyle=:dash, linewidth=1.5,
        label="Seed median = $(Int(med_seed))")
    vline!(p, [med_full], color=:coral, linestyle=:dash, linewidth=1.5,
        label="Full median = $(Int(med_full))")

    # threshold annotations at key values
    for thresh in [50, 100]
        frac_seed = mean(seed .< thresh)
        hline!(p, [frac_seed], color=:gray60, linestyle=:dot, linewidth=0.8, label=false)
        annotate!(p, thresh, frac_seed + 0.02,
            text("$(round(100*frac_seed, digits=1))% < $thresh", 7, :left, :gray40))
    end

    title!(p, "Pfam 36.0 Family Size Distributions")

    return p
end

p = plot_family_size_ecdf(df)
out_path = joinpath(FIG_DIR, "pfam_family_size_distribution.pdf")
savefig(p, out_path)
@info "Saved figure to $out_path"

# ══════════════════════════════════════════════════════════════════════════════
# Print summary table for SI
# ══════════════════════════════════════════════════════════════════════════════

seed = df.n_seed
full = skipmissing(df.n_full) |> collect

println("\nSummary statistics for SI table:")
println("-"^50)
header = ["Statistic", "Seed", "Full"]
data = [
    "Total families"    string(length(seed))           string(length(full));
    "Median"            string(Int(median(seed)))       string(Int(median(full)));
    "Mean"              string(round(mean(seed), digits=1))  string(round(mean(full), digits=1));
    "Q25"               string(Int(quantile(seed, 0.25)))   string(Int(quantile(full, 0.25)));
    "Q75"               string(Int(quantile(seed, 0.75)))   string(Int(quantile(full, 0.75)));
    "< 20 sequences"    "$(sum(seed .< 20)) ($(round(100*mean(seed .< 20), digits=1))%)"   "$(sum(full .< 20)) ($(round(100*mean(full .< 20), digits=1))%)";
    "< 50 sequences"    "$(sum(seed .< 50)) ($(round(100*mean(seed .< 50), digits=1))%)"   "$(sum(full .< 50)) ($(round(100*mean(full .< 50), digits=1))%)";
    "< 100 sequences"   "$(sum(seed .< 100)) ($(round(100*mean(seed .< 100), digits=1))%)" "$(sum(full .< 100)) ($(round(100*mean(full .< 100), digits=1))%)";
]
pretty_table(data, header)
