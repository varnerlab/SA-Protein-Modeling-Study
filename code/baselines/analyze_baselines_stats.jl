#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Baseline Analysis with Statistical Tests + Significance Brackets
#
# Reads per-chain estimates from the baseline comparison, runs Wilcoxon
# rank-sum tests (SA gen vs each baseline), and regenerates figures with
# significance brackets.
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment …"
const _EXPERIMENT_ROOT = @__DIR__

using Pkg
Pkg.activate(joinpath(_EXPERIMENT_ROOT, ".."))

include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))
using LaTeXStrings
using Printf
using HypothesisTests

const BASELINE_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "data")
const PFAM_DATA_DIR     = joinpath(_EXPERIMENT_ROOT, "..", "pfam-families", "data")
const RESULTS_DIR       = joinpath(_EXPERIMENT_ROOT, "results")
const FIG_DIR           = joinpath(_EXPERIMENT_ROOT, "figs")
mkpath(RESULTS_DIR)
mkpath(FIG_DIR)

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

const BASELINE_FILES = [
    ("hmm_sequences.fasta",     "HMM emit"),
    ("evodiff_sequences.fasta", "EvoDiff"),
    ("msat_sequences.fasta",    "MSA Transformer"),
]

const SA_FILES = [
    ("sa_generation_sequences.fasta", "SA generation"),
    ("sa_retrieval_sequences.fasta",  "SA retrieval"),
]

const N_CHAINS = 30
const SPC = 5  # samples per chain

function pval_stars(p::Float64)
    p < 0.001 && return "***"
    p < 0.01  && return "**"
    p < 0.05  && return "*"
    return "n.s."
end

function cohens_d(x::Vector{Float64}, y::Vector{Float64})
    nx, ny = length(x), length(y)
    mx, my = mean(x), mean(y)
    sx, sy = std(x), std(y)
    sp = sqrt(((nx-1)*sx^2 + (ny-1)*sy^2) / (nx + ny - 2))
    return sp > 0 ? (mx - my) / sp : 0.0
end

function load_and_clean_seqs(fasta_path::String, L::Int)
    data = parse_fasta(fasta_path)
    seqs = String[]
    for (_, raw_seq) in data
        cleaned = filter(c -> c in AA_ALPHABET, uppercase(raw_seq))
        if length(cleaned) >= L
            push!(seqs, cleaned[1:L])
        elseif length(cleaned) > 0
            push!(seqs, cleaned * repeat("A", L - length(cleaned)))
        end
    end
    return seqs
end

# ══════════════════════════════════════════════════════════════════════════════
# Compute per-chain metric vectors (for statistical tests)
# ══════════════════════════════════════════════════════════════════════════════

function compute_chain_metrics(gen_seqs::Vector{String}, stored_seqs::Vector{String},
                                X̂::Matrix{Float64}, pca_model, L::Int)
    # Encode into PCA space
    gen_pca = Vector{Vector{Float64}}()
    for seq in gen_seqs
        x_onehot = zeros(Float64, N_AA * L)
        for pos in 1:min(L, length(seq))
            idx = get(AA_TO_IDX, seq[pos], 0)
            if idx > 0
                x_onehot[(pos-1)*N_AA + idx] = 1.0
            end
        end
        ξ_pca = vec(MultivariateStats.transform(pca_model, x_onehot))
        push!(gen_pca, ξ_pca)
    end

    n_total = length(gen_seqs)
    spc = min(SPC, n_total ÷ N_CHAINS)
    spc < 1 && return (Float64[], Float64[], Float64[])

    chain_kls  = Float64[]
    chain_novs = Float64[]
    chain_sids = Float64[]

    for c in 1:N_CHAINS
        i1 = (c-1)*spc + 1
        i2 = c*spc
        i2 > n_total && break
        ch_seqs = gen_seqs[i1:i2]
        ch_pca  = gen_pca[i1:i2]
        push!(chain_kls, aa_composition_kl(ch_seqs, stored_seqs))
        push!(chain_novs, mean(sample_novelty(ξ, X̂) for ξ in ch_pca))
        push!(chain_sids, mean(nearest_sequence_identity(s, stored_seqs) for s in ch_seqs))
    end

    return chain_kls, chain_novs, chain_sids
end

# ══════════════════════════════════════════════════════════════════════════════
# Main: collect chain-level metrics for all methods
# ══════════════════════════════════════════════════════════════════════════════

"""
    bootstrap_kl(gen_seqs, stored_seqs; n_boot=1000) -> (kl_pooled, kl_se)

Compute pooled KL from all sequences and bootstrap SE.
"""
function bootstrap_kl(gen_seqs::Vector{String}, stored_seqs::Vector{String}; n_boot::Int=1000)
    kl_pooled = aa_composition_kl(gen_seqs, stored_seqs)
    n = length(gen_seqs)
    kl_boots = Float64[]
    for _ in 1:n_boot
        boot_idx = rand(1:n, n)
        boot_seqs = gen_seqs[boot_idx]
        push!(kl_boots, aa_composition_kl(boot_seqs, stored_seqs))
    end
    return (kl_pooled, std(kl_boots))
end

function collect_all_chain_metrics()
    # family → method → (chain_kls, chain_novs, chain_sids)
    all_chains = Dict{String, Dict{String, Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}}()
    # family → method → (pooled_kl, kl_se) — correct KL from all sequences
    pooled_kls = Dict{String, Dict{String, Tuple{Float64, Float64}}}()

    for fam in FAMILIES
        @info "════════════════════════════════════════════════════════════"
        @info "  $(fam.name) ($(fam.id))"
        @info "════════════════════════════════════════════════════════════"

        all_chains[fam.name] = Dict{String, Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}()
        pooled_kls[fam.name] = Dict{String, Tuple{Float64, Float64}}()

        stored_fasta = joinpath(PFAM_DATA_DIR, fam.id, "stored_sequences.fasta")
        if !isfile(stored_fasta); @warn "Missing stored — skipping"; continue; end
        stored_data = parse_fasta(stored_fasta)
        stored_seqs = [s[2] for s in stored_data]

        sto_file = joinpath(PFAM_DATA_DIR, fam.id, "$(fam.id)_seed.sto")
        raw_seqs = parse_stockholm(sto_file)
        char_mat, _ = clean_alignment(raw_seqs)
        K, L = size(char_mat)
        X̂, pca_model, L, _ = build_memory_matrix(char_mat; pratio=0.95)

        # helper to process a set of sequences
        function process_method(gen_seqs, method_name)
            chains = compute_chain_metrics(gen_seqs, stored_seqs, X̂, pca_model, L)
            all_chains[fam.name][method_name] = chains
            pooled_kls[fam.name][method_name] = bootstrap_kl(gen_seqs, stored_seqs)
            @info "  $method_name: $(length(chains[1])) chains, pooled KL=$(round(pooled_kls[fam.name][method_name][1], digits=4))"
        end

        # SA methods
        for (fname, method_name) in SA_FILES
            fpath = joinpath(PFAM_DATA_DIR, fam.id, fname)
            if !isfile(fpath); continue; end
            sa_data = parse_fasta(fpath)
            sa_seqs = [s[2] for s in sa_data]
            process_method(sa_seqs, method_name)
        end

        # Bootstrap
        Random.seed!(12345)
        bs_seqs = [stored_seqs[rand(1:length(stored_seqs))] for _ in 1:150]
        process_method(bs_seqs, "Bootstrap")

        # Learned baselines
        for (fname, method_name) in BASELINE_FILES
            fpath = joinpath(BASELINE_DATA_DIR, fam.id, fname)
            if !isfile(fpath); @info "  $method_name: not found"; continue; end
            gen_seqs = load_and_clean_seqs(fpath, L)
            if isempty(gen_seqs); continue; end
            process_method(gen_seqs, method_name)
        end
    end

    return all_chains, pooled_kls
end

# ══════════════════════════════════════════════════════════════════════════════
# Statistical tests: SA generation vs each other method
# ══════════════════════════════════════════════════════════════════════════════

function run_baseline_stats(all_chains)
    @info "\n══════════════════════════════════════════════════════════"
    @info "  Statistical Tests: SA generation vs Baselines"
    @info "══════════════════════════════════════════════════════════\n"

    test_rows = Vector{NamedTuple}()
    comparisons = ["HMM emit", "EvoDiff", "MSA Transformer", "Bootstrap"]

    for fam in FAMILIES
        fam_chains = get(all_chains, fam.name, nothing)
        fam_chains === nothing && continue
        sa_gen = get(fam_chains, "SA generation", nothing)
        sa_gen === nothing && continue

        for comp in comparisons
            comp_chains = get(fam_chains, comp, nothing)
            comp_chains === nothing && continue

            for (metric_idx, metric_name) in enumerate(["KL", "Nov", "SID"])
                sa_vals = sa_gen[metric_idx]
                co_vals = comp_chains[metric_idx]
                if length(sa_vals) >= 5 && length(co_vals) >= 5
                    p = pvalue(MannWhitneyUTest(sa_vals, co_vals))
                    d = cohens_d(sa_vals, co_vals)
                    stars = pval_stars(p)
                    @info @sprintf("  %s %s: SA gen vs %s: p=%.4f (%s), d=%.2f  [SA=%.4f±%.4f, %s=%.4f±%.4f]",
                        fam.name, metric_name, comp, p, stars, d,
                        mean(sa_vals), std(sa_vals)/sqrt(length(sa_vals)),
                        comp, mean(co_vals), std(co_vals)/sqrt(length(co_vals)))

                    push!(test_rows, (
                        Family = fam.name,
                        Metric = metric_name,
                        Method1 = "SA generation",
                        Method2 = comp,
                        Mean1 = mean(sa_vals),
                        SE1 = std(sa_vals)/sqrt(length(sa_vals)),
                        Mean2 = mean(co_vals),
                        SE2 = std(co_vals)/sqrt(length(co_vals)),
                        p_value = p,
                        Cohens_d = d,
                        Stars = stars,
                    ))
                end
            end
        end
    end

    df_tests = DataFrame(test_rows)
    CSV.write(joinpath(RESULTS_DIR, "baseline_statistical_tests.csv"), df_tests)
    @info "Saved baseline_statistical_tests.csv"
    return df_tests
end

# ══════════════════════════════════════════════════════════════════════════════
# Figures with significance brackets
# ══════════════════════════════════════════════════════════════════════════════

function draw_bracket!(p, x1, x2, y, dy, label; color=:black, fontsize=6)
    plot!(p, [x1, x1, x2, x2], [y, y+dy, y+dy, y],
          linecolor=color, linewidth=0.7, label=false)
    annotate!(p, [(mean([x1, x2]), y + dy + 0.4*dy, text(label, fontsize, :center, color))])
end

function generate_figures_with_brackets(all_chains, df_tests, pooled_kls)
    family_names = [f.name for f in FAMILIES]

    # Determine which methods are available
    all_methods_seen = Set{String}()
    for (_, fam_chains) in all_chains
        for (m, _) in fam_chains
            push!(all_methods_seen, m)
        end
    end

    # Show only the four generative methods (SA gen + 3 learned baselines).
    # Simple baselines (bootstrap, Gaussian perturbation, convex combination)
    # and SA retrieval are reported in SI Table S2.
    method_order_full = ["SA generation", "HMM emit", "EvoDiff", "MSA Transformer"]
    method_order = [m for m in method_order_full if m in all_methods_seen]

    method_colors = Dict(
        "SA generation"    => RGB(0.20, 0.47, 0.69),
        "HMM emit"         => RGB(0.30, 0.69, 0.29),
        "EvoDiff"          => RGB(0.93, 0.68, 0.20),
        "MSA Transformer"  => RGB(0.58, 0.32, 0.68),
    )

    n_methods = length(method_order)
    bar_width = 0.8 / n_methods
    offsets = Dict(m => -(n_methods-1)/2 * bar_width + (i-1) * bar_width
                   for (i, m) in enumerate(method_order))

    panel_defaults = (
        fontfamily   = "Helvetica",
        guidefontsize  = 10,
        tickfontsize   = 8,
        legendfontsize = 6,
        titlefontsize  = 11,
        grid           = :y,
        gridalpha      = 0.15,
        framestyle     = :box,
        foreground_color_border = :gray40,
        margin         = 3Plots.mm,
        bottom_margin  = 5Plots.mm,
        left_margin    = 4Plots.mm,
    )

    # Helper to get mean/se from chain data (for novelty and SID)
    # For KL (metric_idx=1), use pooled values instead of per-chain means
    function get_mean_se(fam_name, method, metric_idx)
        if metric_idx == 1
            # Use pooled KL (computed from all 150 sequences) with bootstrap SE
            fam_kls = get(pooled_kls, fam_name, nothing)
            fam_kls === nothing && return (NaN, 0.0)
            kl_data = get(fam_kls, method, nothing)
            kl_data === nothing && return (NaN, 0.0)
            return kl_data  # (pooled_kl, bootstrap_se)
        end
        fam_chains = get(all_chains, fam_name, nothing)
        fam_chains === nothing && return (NaN, 0.0)
        chains = get(fam_chains, method, nothing)
        chains === nothing && return (NaN, 0.0)
        vals = chains[metric_idx]
        isempty(vals) && return (NaN, 0.0)
        return (mean(vals), std(vals)/sqrt(length(vals)))
    end

    # For each metric, find the comparison to annotate (SA gen vs best learned baseline)
    # We'll annotate SA gen vs HMM emit as the most direct comparison
    # (both use only the seed alignment)

    metrics = [
        (idx=1, ylabel=L"\mathrm{KL\ divergence}",        title="(A)  AA composition fidelity (lower is better)",    ylims_base=(0, 0.08)),
        (idx=2, ylabel=L"\mathrm{Novelty}",                title="(B)  Sample novelty (higher is better)",              ylims_base=(0, 0.85)),
        (idx=3, ylabel=L"\mathrm{Nearest\ seq.\ identity}", title="(C)  Sequence identity vs. natural nearest-neighbor envelope", ylims_base=(0, 1.15)),
    ]

    panels = []

    for (mi, met) in enumerate(metrics)
        # Compute max bar height for ylim
        max_val = 0.0
        for fam_name in family_names
            for method in method_order
                m, s = get_mean_se(fam_name, method, met.idx)
                !isnan(m) && (max_val = max(max_val, m + s))
            end
        end

        # Add headroom for up to 3 stacked brackets
        ylim_top = max_val * 1.55
        ylim_top = max(ylim_top, met.ylims_base[2])

        p = plot(; xlabel="", ylabel=met.ylabel, title=met.title,
                 legend=(mi == 1 ? :topright : :none),
                 ylims=(met.ylims_base[1], ylim_top),
                 panel_defaults...)

        for method in method_order
            xs = Float64[]; ys = Float64[]; es = Float64[]
            for (fi, fam_name) in enumerate(family_names)
                m, s = get_mean_se(fam_name, method, met.idx)
                push!(xs, fi + offsets[method])
                push!(ys, isnan(m) ? 0.0 : m)
                push!(es, s)
            end
            bar!(p, xs, ys; bar_width=bar_width, yerror=es,
                 label=method, color=method_colors[method], alpha=0.85, lw=0.5, linecolor=:gray30)
        end
        plot!(p, xticks=(1:length(family_names), family_names), xrotation=0)

        # Add KL = 0.05 reference line on panel A only
        if met.idx == 1
            hline!(p, [0.05]; color=:black, linestyle=:dash, lw=1.0, label="KL = 0.05")
        end

        # Per-family natural within-family nearest-neighbor identity envelope on the
        # SID panel, replacing the old flat 30-70% band. Dark box = p10-p90; light box
        # = full [min, max] range. Data: analysis/data/review_identity_band.csv.
        if met.idx == 3
            env = CSV.read(joinpath(_EXPERIMENT_ROOT, "..", "analysis", "data", "review_identity_band.csv"), DataFrame)
            env_by = Dict(String(r.Family) => r for r in eachrow(env))
            hw = 0.46
            firstlab = true
            for (fi, fam_name) in enumerate(family_names)
                haskey(env_by, fam_name) || continue
                r = env_by[fam_name]
                outer = Shape([fi-hw, fi+hw, fi+hw, fi-hw], [r.NN_min, r.NN_min, r.NN_max, r.NN_max])
                inner = Shape([fi-hw, fi+hw, fi+hw, fi-hw], [r.NN_p10, r.NN_p10, r.NN_p90, r.NN_p90])
                plot!(p, outer; fillalpha=0.06, color=:darkgreen, lw=0,
                      label=(firstlab ? "natural NN range (min-max)" : ""))
                plot!(p, inner; fillalpha=0.15, color=:darkgreen, lw=0,
                      label=(firstlab ? "natural NN envelope (p10-p90)" : ""))
                firstlab = false
            end
        end

        # Add significance brackets: SA gen vs each baseline, stacked by height
        metric_name = ["KL", "Nov", "SID"][met.idx]
        bracket_comps = ["HMM emit", "EvoDiff", "MSA Transformer"]
        dy = max_val * 0.02

        for (fi, fam_name) in enumerate(family_names)
            # compute the base height (top of tallest bar in this family group)
            base_top = 0.0
            for method in method_order
                m, s = get_mean_se(fam_name, method, met.idx)
                !isnan(m) && (base_top = max(base_top, m + s))
            end

            for (bi, bracket_comp) in enumerate(bracket_comps)
                row = filter(r -> r.Family == fam_name && r.Metric == metric_name &&
                                  r.Method2 == bracket_comp, df_tests)
                nrow(row) == 0 && continue
                stars = row.Stars[1]
                stars == "n.s." && continue   # skip non-significant
                x1 = fi + offsets["SA generation"]
                x2 = fi + offsets[bracket_comp]
                top_bar = base_top + (bi * max_val * 0.08)
                draw_bracket!(p, x1, x2, top_bar, dy, stars)
            end
        end

        push!(panels, p)
    end

    composite = plot(panels...; layout=(3, 1), size=(800, 1050), dpi=300)

    # Save as PNG first (GR PDF can timeout on complex plots), then convert
    png_path = joinpath(FIG_DIR, "baseline_comparison_composite.png")
    pdf_path = joinpath(FIG_DIR, "baseline_comparison_composite.pdf")
    savefig(composite, png_path)
    @info "Saved baseline_comparison_composite.png"
    # Also try PDF with individual panels
    try
        savefig(composite, pdf_path)
        @info "Saved baseline_comparison_composite.pdf"
    catch e
        @warn "PDF save failed ($(e)), using PNG. Convert manually if needed."
    end

    # Mirror into both manuscript figure directories. pnas/figs is a real
    # directory (no longer a symlink to paper/figs), so it must be written
    # to explicitly in addition to paper/figs.
    for dest in (joinpath(_EXPERIMENT_ROOT, "..", "..", "paper", "figs"),
                 joinpath(_EXPERIMENT_ROOT, "..", "..", "pnas", "figs"))
        if isdir(dest)
            cp(joinpath(FIG_DIR, "baseline_comparison_composite.pdf"),
               joinpath(dest, "baseline_comparison_composite.pdf"), force=true)
            @info "Copied baseline_comparison_composite.pdf -> $dest"
        end
    end

    return composite
end

# ══════════════════════════════════════════════════════════════════════════════
# Run
# ══════════════════════════════════════════════════════════════════════════════

(all_chains, pooled_kls) = collect_all_chain_metrics()
df_tests = run_baseline_stats(all_chains)
generate_figures_with_brackets(all_chains, df_tests, pooled_kls)

@info "\nDone."
