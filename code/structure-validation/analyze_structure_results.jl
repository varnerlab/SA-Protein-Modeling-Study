#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Structure Validation Analysis: Statistical tests + publication figures
#
# Reads cached PDB files from run_structure_validation.jl, extracts per-sequence
# pLDDT and TM-scores, runs Wilcoxon rank-sum tests (SA gen vs Stored),
# computes effect sizes (Cohen's d), and produces figures with significance
# brackets.
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment …"
const _EXPERIMENT_ROOT = @__DIR__

using Pkg
Pkg.activate(joinpath(_EXPERIMENT_ROOT, ".."))

include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))

using Printf
using HypothesisTests

const STRUCT_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "data")
const RESULTS_DIR = joinpath(_EXPERIMENT_ROOT, "results")
const FIGS_DIR = joinpath(_EXPERIMENT_ROOT, "figs")
const TMALIGN_BIN = joinpath(_EXPERIMENT_ROOT, "bin", "TMalign")

const FAMILIES = [
    (id="PF00076", name="RRM",     pdb="1FXL", chain="A"),
    (id="PF00018", name="SH3",     pdb="1SHG", chain="A"),
    (id="PF00397", name="WW",      pdb="1PIN", chain="A"),
    (id="PF00014", name="Kunitz",  pdb="1BPI", chain="A"),
    (id="PF00096", name="zf-C2H2", pdb="1ZAA", chain="C"),
    (id="PF00595", name="PDZ",     pdb="1BE9", chain="A"),
    (id="PF00069", name="Pkinase", pdb="1ATP", chain="E"),
    (id="PF00711", name="Defensin_beta", pdb="1E4S", chain="A"),
]

const METHODS = ["SA generation", "SA retrieval", "Bootstrap", "Stored"]
const METHOD_DIRS = ["sa_generation", "sa_retrieval", "bootstrap", "stored"]

# ══════════════════════════════════════════════════════════════════════════════
# Helper functions
# ══════════════════════════════════════════════════════════════════════════════

function extract_mean_plddt(pdb_string::String)
    plddts = Float64[]
    for line in split(pdb_string, '\n')
        if startswith(line, "ATOM") && length(line) >= 66
            try
                plddt = parse(Float64, strip(line[61:66]))
                push!(plddts, plddt)
            catch; continue; end
        end
    end
    if isempty(plddts)
        return NaN
    end
    m = mean(plddts)
    return m < 1.5 ? m * 100.0 : m
end

function run_tmalign(query_pdb::String, reference_pdb::String)
    try
        output = read(`$TMALIGN_BIN $query_pdb $reference_pdb`, String)
        for line in split(output, '\n')
            if contains(line, "TM-score") && contains(line, "Chain_2")
                m = match(r"TM-score=\s*([\d.]+)", line)
                m !== nothing && return parse(Float64, m.captures[1])
            end
        end
        for line in split(output, '\n')
            if contains(line, "TM-score")
                m = match(r"TM-score=\s*([\d.]+)", line)
                m !== nothing && return parse(Float64, m.captures[1])
            end
        end
    catch e
        @warn "TMalign failed: $e"
    end
    return NaN
end

function cohens_d(x::Vector{Float64}, y::Vector{Float64})
    nx, ny = length(x), length(y)
    mx, my = mean(x), mean(y)
    sx, sy = std(x), std(y)
    sp = sqrt(((nx-1)*sx^2 + (ny-1)*sy^2) / (nx + ny - 2))
    return sp > 0 ? (mx - my) / sp : 0.0
end

function pval_stars(p::Float64)
    p < 0.001 && return "***"
    p < 0.01  && return "**"
    p < 0.05  && return "*"
    return "n.s."
end

# ══════════════════════════════════════════════════════════════════════════════
# Extract raw per-sequence metrics from cached PDB files
# ══════════════════════════════════════════════════════════════════════════════

function collect_raw_metrics()
    # Store raw per-sequence values: family → method → vector of values
    plddt_data = Dict{String, Dict{String, Vector{Float64}}}()
    tm_data    = Dict{String, Dict{String, Vector{Float64}}}()

    for fam in FAMILIES
        @info "Collecting metrics for $(fam.name) ($(fam.id)) …"
        plddt_data[fam.name] = Dict{String, Vector{Float64}}()
        tm_data[fam.name]    = Dict{String, Vector{Float64}}()

        ref_pdb = joinpath(STRUCT_DATA_DIR, fam.id, "reference", "$(fam.pdb)_$(fam.chain).pdb")
        if !isfile(ref_pdb)
            @warn "  Missing reference PDB: $ref_pdb"
            continue
        end

        for (method, method_dir) in zip(METHODS, METHOD_DIRS)
            struct_dir = joinpath(STRUCT_DATA_DIR, fam.id, "structures", method_dir)
            if !isdir(struct_dir)
                @warn "  Missing structure dir: $struct_dir"
                plddt_data[fam.name][method] = Float64[]
                tm_data[fam.name][method]    = Float64[]
                continue
            end

            pdb_files = filter(f -> endswith(f, ".pdb"), readdir(struct_dir, join=true))
            plddts = Float64[]
            tms = Float64[]

            for pdb_file in pdb_files
                pdb_content = read(pdb_file, String)
                plddt = extract_mean_plddt(pdb_content)
                !isnan(plddt) && push!(plddts, plddt)

                tm = run_tmalign(pdb_file, ref_pdb)
                !isnan(tm) && push!(tms, tm)
            end

            plddt_data[fam.name][method] = plddts
            tm_data[fam.name][method]    = tms
            @info @sprintf("  %s: n=%d, pLDDT=%.1f±%.1f, TM=%.3f±%.3f",
                method, length(plddts),
                isempty(plddts) ? 0.0 : mean(plddts),
                isempty(plddts) ? 0.0 : std(plddts)/sqrt(length(plddts)),
                isempty(tms) ? 0.0 : mean(tms),
                isempty(tms) ? 0.0 : std(tms)/sqrt(length(tms)))
        end
    end

    return plddt_data, tm_data
end

# ══════════════════════════════════════════════════════════════════════════════
# Statistical tests: SA generation vs Stored (primary comparison)
# ══════════════════════════════════════════════════════════════════════════════

function run_statistical_tests(plddt_data, tm_data)
    @info "\n══════════════════════════════════════════════════════════"
    @info "  Statistical Tests: SA generation vs Stored"
    @info "══════════════════════════════════════════════════════════\n"

    test_rows = Vector{NamedTuple}()

    for fam in FAMILIES
        fam_name = fam.name

        # pLDDT comparison
        gen_plddt = get(plddt_data[fam_name], "SA generation", Float64[])
        sto_plddt = get(plddt_data[fam_name], "Stored", Float64[])

        if length(gen_plddt) >= 5 && length(sto_plddt) >= 5
            mw_plddt = MannWhitneyUTest(gen_plddt, sto_plddt)
            p_plddt = pvalue(mw_plddt)
            d_plddt = cohens_d(gen_plddt, sto_plddt)

            @info @sprintf("  %s pLDDT: SA gen=%.1f±%.1f, Stored=%.1f±%.1f, p=%.4f (%s), d=%.2f",
                fam_name, mean(gen_plddt), std(gen_plddt)/sqrt(length(gen_plddt)),
                mean(sto_plddt), std(sto_plddt)/sqrt(length(sto_plddt)),
                p_plddt, pval_stars(p_plddt), d_plddt)
        else
            p_plddt = NaN; d_plddt = NaN
        end

        # TM-score comparison
        gen_tm = get(tm_data[fam_name], "SA generation", Float64[])
        sto_tm = get(tm_data[fam_name], "Stored", Float64[])

        if length(gen_tm) >= 5 && length(sto_tm) >= 5
            mw_tm = MannWhitneyUTest(gen_tm, sto_tm)
            p_tm = pvalue(mw_tm)
            d_tm = cohens_d(gen_tm, sto_tm)

            @info @sprintf("  %s TM:    SA gen=%.3f±%.3f, Stored=%.3f±%.3f, p=%.4f (%s), d=%.2f",
                fam_name, mean(gen_tm), std(gen_tm)/sqrt(length(gen_tm)),
                mean(sto_tm), std(sto_tm)/sqrt(length(sto_tm)),
                p_tm, pval_stars(p_tm), d_tm)
        else
            p_tm = NaN; d_tm = NaN
        end

        push!(test_rows, (
            Family = fam_name,
            pLDDT_gen_mean = isempty(gen_plddt) ? NaN : mean(gen_plddt),
            pLDDT_gen_se   = isempty(gen_plddt) ? NaN : std(gen_plddt)/sqrt(length(gen_plddt)),
            pLDDT_sto_mean = isempty(sto_plddt) ? NaN : mean(sto_plddt),
            pLDDT_sto_se   = isempty(sto_plddt) ? NaN : std(sto_plddt)/sqrt(length(sto_plddt)),
            pLDDT_p        = p_plddt,
            pLDDT_d        = d_plddt,
            TM_gen_mean    = isempty(gen_tm) ? NaN : mean(gen_tm),
            TM_gen_se      = isempty(gen_tm) ? NaN : std(gen_tm)/sqrt(length(gen_tm)),
            TM_sto_mean    = isempty(sto_tm) ? NaN : mean(sto_tm),
            TM_sto_se      = isempty(sto_tm) ? NaN : std(sto_tm)/sqrt(length(sto_tm)),
            TM_p           = p_tm,
            TM_d           = d_tm,
        ))
    end

    # Also test SA generation vs Bootstrap for each family
    @info "\n══════════════════════════════════════════════════════════"
    @info "  Statistical Tests: SA generation vs Bootstrap"
    @info "══════════════════════════════════════════════════════════\n"

    for fam in FAMILIES
        fam_name = fam.name
        gen_plddt = get(plddt_data[fam_name], "SA generation", Float64[])
        bs_plddt  = get(plddt_data[fam_name], "Bootstrap", Float64[])
        gen_tm    = get(tm_data[fam_name], "SA generation", Float64[])
        bs_tm     = get(tm_data[fam_name], "Bootstrap", Float64[])

        if length(gen_plddt) >= 5 && length(bs_plddt) >= 5
            p = pvalue(MannWhitneyUTest(gen_plddt, bs_plddt))
            d = cohens_d(gen_plddt, bs_plddt)
            @info @sprintf("  %s pLDDT: SA gen vs BS: p=%.4f (%s), d=%.2f", fam_name, p, pval_stars(p), d)
        end
        if length(gen_tm) >= 5 && length(bs_tm) >= 5
            p = pvalue(MannWhitneyUTest(gen_tm, bs_tm))
            d = cohens_d(gen_tm, bs_tm)
            @info @sprintf("  %s TM:    SA gen vs BS: p=%.4f (%s), d=%.2f", fam_name, p, pval_stars(p), d)
        end
    end

    return DataFrame(test_rows)
end

# ══════════════════════════════════════════════════════════════════════════════
# Figures with significance brackets
# ══════════════════════════════════════════════════════════════════════════════

function draw_significance_bracket!(p, x1, x2, y, dy, label; color=:black, fontsize=7)
    # Draw bracket: two vertical ticks + horizontal bar + label
    plot!(p, [x1, x1, x2, x2], [y, y+dy, y+dy, y],
          linecolor=color, linewidth=0.8, label=false)
    annotate!(p, [(mean([x1, x2]), y + dy + 0.3*dy, text(label, fontsize, :center, color))])
end

function generate_figures_with_stats(plddt_data, tm_data)
    families = [fam.name for fam in FAMILIES]
    method_colors = [
        RGB(0.20, 0.47, 0.69),   # SA gen - steel blue (canonical)
        RGB(0.89, 0.44, 0.32),   # SA ret - coral (canonical)
        RGB(0.50, 0.50, 0.50),   # Bootstrap - gray (canonical)
        RGB(0.20, 0.63, 0.17),   # Stored - green (canonical)
    ]

    n_fam = length(families)
    bar_width = 0.18
    offsets = [-1.5, -0.5, 0.5, 1.5] .* bar_width

    # ── Panel 1: pLDDT with significance brackets ───────────────────────
    p1 = plot(size=(560, 380), ylabel="Mean pLDDT",
              xticks=(1:n_fam, families), xrotation=30,
              legend=:bottomright, legendfontsize=6,
              ylims=(0, 115), grid=:y,
              left_margin=5Plots.mm, bottom_margin=10Plots.mm,
              top_margin=5Plots.mm)
    hline!(p1, [70], color=:gray, linestyle=:dash, linewidth=0.8, label="pLDDT = 70")

    for (j, method) in enumerate(METHODS)
        xs = Float64[]; ys = Float64[]; es = Float64[]
        for (i, fam) in enumerate(families)
            vals = get(plddt_data[fam], method, Float64[])
            if !isempty(vals)
                push!(xs, i + offsets[j])
                push!(ys, mean(vals))
                push!(es, std(vals)/sqrt(length(vals)))
            end
        end
        bar!(p1, xs, ys, yerror=es, bar_width=bar_width,
             color=method_colors[j], label=method, linewidth=0.5)
    end

    # Add significance brackets: SA gen (offset[1]) vs Stored (offset[4])
    for (i, fam) in enumerate(families)
        gen_vals = get(plddt_data[fam], "SA generation", Float64[])
        sto_vals = get(plddt_data[fam], "Stored", Float64[])
        if length(gen_vals) >= 5 && length(sto_vals) >= 5
            p = pvalue(MannWhitneyUTest(gen_vals, sto_vals))
            label = pval_stars(p)
            max_y = max(mean(gen_vals) + std(gen_vals)/sqrt(length(gen_vals)),
                        mean(sto_vals) + std(sto_vals)/sqrt(length(sto_vals)))
            bracket_y = max_y + 2.0
            x1 = i + offsets[1]
            x2 = i + offsets[4]
            draw_significance_bracket!(p1, x1, x2, bracket_y, 1.5, label; fontsize=6)
        end
    end

    # ── Panel 2: TM-score with significance brackets ────────────────────
    p2 = plot(size=(560, 380), ylabel="TM-score",
              xticks=(1:n_fam, families), xrotation=30,
              legend=:bottomright, legendfontsize=6,
              ylims=(0, 1.15), grid=:y,
              left_margin=5Plots.mm, bottom_margin=10Plots.mm,
              top_margin=5Plots.mm)
    hline!(p2, [0.5], color=:gray, linestyle=:dash, linewidth=0.8, label="TM = 0.5")

    for (j, method) in enumerate(METHODS)
        xs = Float64[]; ys = Float64[]; es = Float64[]
        for (i, fam) in enumerate(families)
            vals = get(tm_data[fam], method, Float64[])
            if !isempty(vals)
                push!(xs, i + offsets[j])
                push!(ys, mean(vals))
                push!(es, std(vals)/sqrt(length(vals)))
            end
        end
        bar!(p2, xs, ys, yerror=es, bar_width=bar_width,
             color=method_colors[j], label=method, linewidth=0.5)
    end

    # Add significance brackets for TM-score
    for (i, fam) in enumerate(families)
        gen_vals = get(tm_data[fam], "SA generation", Float64[])
        sto_vals = get(tm_data[fam], "Stored", Float64[])
        if length(gen_vals) >= 5 && length(sto_vals) >= 5
            p = pvalue(MannWhitneyUTest(gen_vals, sto_vals))
            label = pval_stars(p)
            max_y = max(mean(gen_vals) + std(gen_vals)/sqrt(length(gen_vals)),
                        mean(sto_vals) + std(sto_vals)/sqrt(length(sto_vals)))
            bracket_y = max_y + 0.02
            x1 = i + offsets[1]
            x2 = i + offsets[4]
            draw_significance_bracket!(p2, x1, x2, bracket_y, 0.015, label; fontsize=6)
        end
    end

    # ── Composite ───────────────────────────────────────────────────────
    composite = plot(p1, p2, layout=(1, 2), size=(1100, 420),
                     title=["(A) ESMFold Confidence (pLDDT, higher is better)" "(B) Structural Similarity (TM-score, higher is better)"],
                     titlefontsize=10)

    mkpath(FIGS_DIR)
    savefig(composite, joinpath(FIGS_DIR, "structure_validation_composite.pdf"))
    @info "Figure saved: $(joinpath(FIGS_DIR, "structure_validation_composite.pdf"))"

    # Copy to paper
    paper_figs = joinpath(_EXPERIMENT_ROOT, "..", "..", "arxiv", "figs")
    if isdir(paper_figs)
        cp(joinpath(FIGS_DIR, "structure_validation_composite.pdf"),
           joinpath(paper_figs, "structure_validation_composite.pdf"), force=true)
        @info "Copied to arxiv/figs/"
    end

    return composite
end

# ══════════════════════════════════════════════════════════════════════════════
# Generate LaTeX-ready results table
# ══════════════════════════════════════════════════════════════════════════════

function generate_latex_table(plddt_data, tm_data)
    @info "\n══════════════════════════════════════════════════════════"
    @info "  LaTeX Table: Structure Validation Results"
    @info "══════════════════════════════════════════════════════════\n"

    for fam in FAMILIES
        for method in METHODS
            pv = get(plddt_data[fam.name], method, Float64[])
            tv = get(tm_data[fam.name], method, Float64[])
            if isempty(pv); continue; end

            pm = mean(pv); pse = std(pv)/sqrt(length(pv))
            p70 = mean(x -> x > 70 ? 1.0 : 0.0, pv)

            if !isempty(tv)
                tm = mean(tv); tse = std(tv)/sqrt(length(tv))
                t50 = mean(x -> x > 0.5 ? 1.0 : 0.0, tv)
                @info @sprintf("  %s & %s & \$%.1f \\pm %.1f\$ & %.0f\\%% & \$%.3f \\pm %.3f\$ & %.0f\\%%",
                    fam.name, method, pm, pse, 100*p70, tm, tse, 100*t50)
            else
                @info @sprintf("  %s & %s & \$%.1f \\pm %.1f\$ & %.0f\\%% & --- & ---",
                    fam.name, method, pm, pse, 100*p70)
            end
        end
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

@info "Collecting raw per-sequence metrics from cached PDB files …"
plddt_data, tm_data = collect_raw_metrics()

@info "\nRunning statistical tests …"
df_tests = run_statistical_tests(plddt_data, tm_data)
mkpath(RESULTS_DIR)
CSV.write(joinpath(RESULTS_DIR, "statistical_tests.csv"), df_tests)
@info "Statistical tests saved to results/statistical_tests.csv"

@info "\nGenerating figures with significance brackets …"
generate_figures_with_stats(plddt_data, tm_data)

generate_latex_table(plddt_data, tm_data)

@info "\nDone."
