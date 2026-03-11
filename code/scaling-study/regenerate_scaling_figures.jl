#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Regenerate scaling study figures — PNAS quality with SE error bands
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment …"
const _EXPERIMENT_ROOT = @__DIR__
include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))

const SCALE_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "data")
const SCALE_FIG_DIR  = joinpath(_EXPERIMENT_ROOT, "figs")
mkpath(SCALE_FIG_DIR)

# Load existing results
df = CSV.read(joinpath(SCALE_DATA_DIR, "scaling_results.csv"), DataFrame)
@info "Loaded $(nrow(df)) rows from scaling_results.csv"

n_reps = 5  # N_REPEATS used in the study

# aggregate: mean ± SE over repeats
gdf = groupby(df, [:K, :Method])
agg = combine(gdf,
    :KL_AA    => mean => :KL_mean,     :KL_AA    => (x -> std(x)/sqrt(n_reps)) => :KL_se,
    :Novelty  => mean => :Nov_mean,    :Novelty  => (x -> std(x)/sqrt(n_reps)) => :Nov_se,
    :SeqID    => mean => :SID_mean,    :SeqID    => (x -> std(x)/sqrt(n_reps)) => :SID_se,
    :Diversity => mean => :Div_mean,   :Diversity => (x -> std(x)/sqrt(n_reps)) => :Div_se,
    :β_star   => mean => :β_star_mean, :β_star   => (x -> std(x)/sqrt(n_reps)) => :β_star_se,
)

# ══════════════════════════════════════════════════════════════════════════════
# PNAS-quality plot defaults
# ══════════════════════════════════════════════════════════════════════════════

const PNAS_SIZE    = (520, 380)
const PNAS_DPI     = 300
const PNAS_LW      = 2.5
const PNAS_MS      = 6
const PNAS_ALPHA   = 0.20  # ribbon fill alpha
const PNAS_GUIDE   = :outertopright

method_colors = Dict(
    "SA (generation)" => RGB(0.20, 0.47, 0.69),   # steel blue
    "SA (retrieval)"  => RGB(0.89, 0.44, 0.32),   # coral
    "Bootstrap"       => RGB(0.50, 0.50, 0.50),   # gray
    "Gaussian pert."  => RGB(0.93, 0.68, 0.20),   # amber
    "Convex comb."    => RGB(0.58, 0.32, 0.68),   # purple
)
method_labels = Dict(
    "SA (generation)" => "SA generation",
    "SA (retrieval)"  => "SA retrieval",
    "Bootstrap"       => "Bootstrap",
    "Gaussian pert."  => "Gaussian pert.",
    "Convex comb."    => "Convex comb.",
)
method_markers = Dict(
    "SA (generation)" => :circle,
    "SA (retrieval)"  => :diamond,
    "Bootstrap"       => :square,
    "Gaussian pert."  => :utriangle,
    "Convex comb."    => :star5,
)
method_order = ["SA (generation)", "SA (retrieval)", "Bootstrap", "Gaussian pert.", "Convex comb."]

pnas_defaults = (
    fontfamily   = "Helvetica",
    guidefontsize  = 11,
    tickfontsize   = 9,
    legendfontsize = 8,
    titlefontsize  = 12,
    grid           = :y,
    gridalpha      = 0.15,
    framestyle     = :box,
    foreground_color_border = :gray40,
    size           = PNAS_SIZE,
    dpi            = PNAS_DPI,
    margin         = 4Plots.mm,
    bottom_margin  = 5Plots.mm,
    left_margin    = 5Plots.mm,
)

# ══════════════════════════════════════════════════════════════════════════════
# Fig 3a — KL divergence vs K  (LOG y-axis so CC doesn't crush everything)
# ══════════════════════════════════════════════════════════════════════════════

K_ticks = ([20, 50, 100, 200, 400], ["20", "50", "100", "200", "400"])

p_kl = plot(;
    xlabel = "Family size (K)",
    ylabel = "KL divergence",
    title  = "Amino acid composition fidelity",
    xscale = :log10,
    yscale = :log10,
    xticks = K_ticks,
    yticks = ([0.01, 0.1, 1.0], ["0.01", "0.1", "1.0"]),
    legend = :topright,
    pnas_defaults...
)

for m in method_order
    sub = sort(filter(r -> r.Method == m, agg), :K)
    isempty(sub) && continue
    # on log scale, use explicit error bars not ribbon (ribbon can go negative)
    plot!(p_kl, sub.K, sub.KL_mean,
          label       = method_labels[m],
          lw          = PNAS_LW,
          color       = method_colors[m],
          markershape = method_markers[m],
          markersize  = PNAS_MS,
          markercolor = method_colors[m],
          yerror      = sub.KL_se,
    )
end
savefig(p_kl, joinpath(SCALE_FIG_DIR, "scaling_kl_vs_k.pdf"))
@info "  Saved scaling_kl_vs_k.pdf"

# ══════════════════════════════════════════════════════════════════════════════
# Fig 3b — Novelty vs K
# ══════════════════════════════════════════════════════════════════════════════

p_nov = plot(;
    xlabel = "Family size (K)",
    ylabel = "Novelty (1 − max cosine sim.)",
    title  = "Sample novelty",
    xscale = :log10,
    xticks = K_ticks,
    legend = :right,
    ylims  = (-0.05, 0.85),
    pnas_defaults...
)

for m in method_order
    sub = sort(filter(r -> r.Method == m, agg), :K)
    isempty(sub) && continue
    plot!(p_nov, sub.K, sub.Nov_mean,
          ribbon      = sub.Nov_se,
          fillalpha   = PNAS_ALPHA,
          label       = method_labels[m],
          lw          = PNAS_LW,
          color       = method_colors[m],
          markershape = method_markers[m],
          markersize  = PNAS_MS,
          markercolor = method_colors[m],
    )
end
savefig(p_nov, joinpath(SCALE_FIG_DIR, "scaling_novelty_vs_k.pdf"))
@info "  Saved scaling_novelty_vs_k.pdf"

# ══════════════════════════════════════════════════════════════════════════════
# Fig 3c — Sequence identity vs K
# ══════════════════════════════════════════════════════════════════════════════

p_sid = plot(;
    xlabel = "Family size (K)",
    ylabel = "Nearest sequence identity",
    title  = "Sequence identity to nearest stored pattern",
    xscale = :log10,
    xticks = K_ticks,
    legend = :right,
    pnas_defaults...
)

for m in method_order
    sub = sort(filter(r -> r.Method == m, agg), :K)
    isempty(sub) && continue
    plot!(p_sid, sub.K, sub.SID_mean,
          ribbon      = sub.SID_se,
          fillalpha   = PNAS_ALPHA,
          label       = method_labels[m],
          lw          = PNAS_LW,
          color       = method_colors[m],
          markershape = method_markers[m],
          markersize  = PNAS_MS,
          markercolor = method_colors[m],
    )
end
savefig(p_sid, joinpath(SCALE_FIG_DIR, "scaling_seqid_vs_k.pdf"))
@info "  Saved scaling_seqid_vs_k.pdf"

# ══════════════════════════════════════════════════════════════════════════════
# Fig 3d — β* vs K
# ══════════════════════════════════════════════════════════════════════════════

sa_gen_agg = sort(filter(r -> r.Method == "SA (generation)", agg), :K)
p_beta = plot(sa_gen_agg.K, sa_gen_agg.β_star_mean;
    ribbon      = sa_gen_agg.β_star_se,
    fillalpha   = PNAS_ALPHA,
    xlabel      = "Family size (K)",
    ylabel      = "β*",
    title       = "Critical temperature",
    xscale      = :log10,
    xticks      = K_ticks,
    lw          = PNAS_LW,
    color       = method_colors["SA (generation)"],
    markershape = :circle,
    markersize  = PNAS_MS,
    markercolor = method_colors["SA (generation)"],
    label       = "β* (empirical)",
    legend      = :topleft,
    pnas_defaults...
)
# reference line: √d prediction for comparison
d_vals = sort(filter(r -> r.Method == "SA (generation)", agg), :K)
savefig(p_beta, joinpath(SCALE_FIG_DIR, "scaling_beta_star_vs_k.pdf"))
@info "  Saved scaling_beta_star_vs_k.pdf"

@info "Done — all scaling figures regenerated (PNAS quality, SE error bands, log KL axis)."
