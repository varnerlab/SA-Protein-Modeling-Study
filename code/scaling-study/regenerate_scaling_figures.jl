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
# Build individual panels, then compose into a 2×2 figure
# ══════════════════════════════════════════════════════════════════════════════

K_ticks = ([20, 50, 100, 200, 400], ["20", "50", "100", "200", "400"])

panel_defaults = (
    fontfamily   = "Helvetica",
    guidefontsize  = 10,
    tickfontsize   = 8,
    legendfontsize = 7,
    titlefontsize  = 11,
    grid           = :y,
    gridalpha      = 0.15,
    framestyle     = :box,
    foreground_color_border = :gray40,
    margin         = 3Plots.mm,
    bottom_margin  = 5Plots.mm,
    left_margin    = 4Plots.mm,
)

# helper for line+ribbon panels
function scaling_panel(agg, method_order, method_labels, method_colors, method_markers;
                       val_col::Symbol, se_col::Symbol,
                       ylabel::String, title::String,
                       use_log_y::Bool=false, legend=:none,
                       ylims=:auto, use_yerror::Bool=false)
    kw = Dict{Symbol,Any}(
        :xlabel => "Family size (K)",
        :ylabel => ylabel,
        :title  => title,
        :xscale => :log10,
        :xticks => K_ticks,
        :legend => legend,
    )
    if use_log_y
        kw[:yscale] = :log10
        kw[:yticks] = ([0.01, 0.1, 1.0], ["0.01", "0.1", "1.0"])
    end
    ylims != :auto && (kw[:ylims] = ylims)

    p = plot(; panel_defaults..., kw...)
    for m in method_order
        sub = sort(filter(r -> r.Method == m, agg), :K)
        isempty(sub) && continue
        if use_yerror
            plot!(p, sub.K, getproperty(sub, val_col);
                  label=method_labels[m], lw=PNAS_LW, color=method_colors[m],
                  markershape=method_markers[m], markersize=PNAS_MS,
                  markercolor=method_colors[m],
                  yerror=getproperty(sub, se_col))
        else
            plot!(p, sub.K, getproperty(sub, val_col);
                  ribbon=getproperty(sub, se_col), fillalpha=PNAS_ALPHA,
                  label=method_labels[m], lw=PNAS_LW, color=method_colors[m],
                  markershape=method_markers[m], markersize=PNAS_MS,
                  markercolor=method_colors[m])
        end
    end
    return p
end

# (A) KL — log-log
p_kl = scaling_panel(agg, method_order, method_labels, method_colors, method_markers;
    val_col=:KL_mean, se_col=:KL_se,
    ylabel="KL divergence", title="(A)  AA composition fidelity",
    use_log_y=true, use_yerror=true, legend=:topright)

# (B) Novelty
p_nov = scaling_panel(agg, method_order, method_labels, method_colors, method_markers;
    val_col=:Nov_mean, se_col=:Nov_se,
    ylabel="Novelty", title="(B)  Sample novelty",
    ylims=(-0.05, 0.85))

# (C) SeqID
p_sid = scaling_panel(agg, method_order, method_labels, method_colors, method_markers;
    val_col=:SID_mean, se_col=:SID_se,
    ylabel="Nearest seq. identity", title="(C)  Sequence identity")

# (D) β* vs K — single method
sa_gen_agg = sort(filter(r -> r.Method == "SA (generation)", agg), :K)
p_beta = plot(sa_gen_agg.K, sa_gen_agg.β_star_mean;
    ribbon      = sa_gen_agg.β_star_se,
    fillalpha   = PNAS_ALPHA,
    xlabel      = "Family size (K)",
    ylabel      = "β*",
    title       = "(D)  Critical temperature",
    xscale      = :log10,
    xticks      = K_ticks,
    lw          = PNAS_LW,
    color       = method_colors["SA (generation)"],
    markershape = :circle,
    markersize  = PNAS_MS,
    markercolor = method_colors["SA (generation)"],
    label       = "β* (empirical)",
    legend      = :topleft,
    panel_defaults...
)

# ── Composite: 2×2 ───────────────────────────────────────────────────────────
p_composite = plot(p_kl, p_nov, p_sid, p_beta;
    layout = (2, 2),
    size   = (900, 700),
    dpi    = 300,
)
savefig(p_composite, joinpath(SCALE_FIG_DIR, "scaling_composite.pdf"))
@info "  Saved scaling_composite.pdf"

# also save individual panels for flexibility
savefig(p_kl,   joinpath(SCALE_FIG_DIR, "scaling_kl_vs_k.pdf"))
savefig(p_nov,  joinpath(SCALE_FIG_DIR, "scaling_novelty_vs_k.pdf"))
savefig(p_sid,  joinpath(SCALE_FIG_DIR, "scaling_seqid_vs_k.pdf"))
savefig(p_beta, joinpath(SCALE_FIG_DIR, "scaling_beta_star_vs_k.pdf"))

@info "Done — all scaling figures regenerated (PNAS quality, SE, composite + individual)."
