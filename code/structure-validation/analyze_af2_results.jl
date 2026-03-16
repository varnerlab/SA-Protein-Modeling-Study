#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# AlphaFold2 (ColabFold) Structure Validation Analysis
#
# Reads ColabFold output PDB files, extracts pLDDT and TM-scores, runs
# statistical tests comparing SA gen vs Stored, and produces figures that
# compare both ESMFold and AlphaFold2 results side by side.
#
# Prerequisites:
#   - Run prepare_colabfold_input.jl to generate FASTA files
#   - Run colabfold_predict.py on Google Colab to generate AF2 predictions
#   - Download colabfold_output/ directory back to this location
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment …"
const _EXPERIMENT_ROOT = @__DIR__

using Pkg
Pkg.activate(joinpath(_EXPERIMENT_ROOT, ".."))

include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))

using Printf
using HypothesisTests

const STRUCT_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "data")          # ESMFold results
const AF2_OUTPUT_DIR  = joinpath(_EXPERIMENT_ROOT, "colabfold_output")  # AF2 results
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

const CATEGORIES = ["sa_generation", "stored"]
const CATEGORY_LABELS = ["SA generation", "Stored"]

# ══════════════════════════════════════════════════════════════════════════════
# Helper functions
# ══════════════════════════════════════════════════════════════════════════════

function extract_mean_plddt(pdb_path::String)
    plddts = Float64[]
    for line in eachline(pdb_path)
        if startswith(line, "ATOM") && length(line) >= 66
            try
                plddt = parse(Float64, strip(line[61:66]))
                push!(plddts, plddt)
            catch; continue; end
        end
    end
    isempty(plddts) && return NaN
    m = mean(plddts)
    # ColabFold outputs pLDDT in [0,100]; ESMFold API outputs [0,1]
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

function draw_significance_bracket!(p, x1, x2, y, dy, label; color=:black, fontsize=7)
    plot!(p, [x1, x1, x2, x2], [y, y+dy, y+dy, y],
          linecolor=color, linewidth=0.8, label=false)
    annotate!(p, [(mean([x1, x2]), y + dy + 0.3*dy, text(label, fontsize, :center, color))])
end

# ══════════════════════════════════════════════════════════════════════════════
# Find ColabFold output PDB for a given sequence ID
# ══════════════════════════════════════════════════════════════════════════════

function find_af2_pdb(output_dir::String, seq_id::String)
    # ColabFold outputs: {id}_unrelaxed_rank_001_alphafold2_ptm_model_*_seed_000.pdb
    if !isdir(output_dir)
        return nothing
    end

    # Try exact match first
    for f in readdir(output_dir)
        if startswith(f, seq_id) && contains(f, "rank_001") && endswith(f, ".pdb")
            return joinpath(output_dir, f)
        end
    end

    # Try with _unrelaxed pattern
    for f in readdir(output_dir)
        if startswith(f, seq_id) && endswith(f, ".pdb") && !contains(f, "rank_002")
            return joinpath(output_dir, f)
        end
    end

    return nothing
end

# ══════════════════════════════════════════════════════════════════════════════
# Collect AF2 metrics
# ══════════════════════════════════════════════════════════════════════════════

function collect_af2_metrics()
    plddt_data = Dict{String, Dict{String, Vector{Float64}}}()
    tm_data    = Dict{String, Dict{String, Vector{Float64}}}()

    for fam in FAMILIES
        @info "Collecting AF2 metrics for $(fam.name) ($(fam.id)) …"
        plddt_data[fam.name] = Dict{String, Vector{Float64}}()
        tm_data[fam.name]    = Dict{String, Vector{Float64}}()

        # Reference PDB (already downloaded by ESMFold experiment)
        ref_pdb = joinpath(STRUCT_DATA_DIR, fam.id, "reference", "$(fam.pdb)_$(fam.chain).pdb")
        if !isfile(ref_pdb)
            @warn "  Missing reference PDB: $ref_pdb"
            continue
        end

        for (category, label) in zip(CATEGORIES, CATEGORY_LABELS)
            af2_dir = joinpath(AF2_OUTPUT_DIR, fam.id, category)
            if !isdir(af2_dir)
                @warn "  Missing AF2 output dir: $af2_dir"
                plddt_data[fam.name][label] = Float64[]
                tm_data[fam.name][label]    = Float64[]
                continue
            end

            plddts = Float64[]
            tms = Float64[]

            # Find all rank_001 PDB files in the output directory
            pdb_files = filter(readdir(af2_dir, join=true)) do f
                endswith(f, ".pdb") && contains(basename(f), "rank_001")
            end

            # If no rank files found, try all PDB files
            if isempty(pdb_files)
                pdb_files = filter(f -> endswith(f, ".pdb"), readdir(af2_dir, join=true))
            end

            for pdb_file in pdb_files
                plddt = extract_mean_plddt(pdb_file)
                !isnan(plddt) && push!(plddts, plddt)

                tm = run_tmalign(pdb_file, ref_pdb)
                !isnan(tm) && push!(tms, tm)
            end

            plddt_data[fam.name][label] = plddts
            tm_data[fam.name][label]    = tms

            @info @sprintf("  %s (%s): n=%d, pLDDT=%.1f±%.1f, TM=%.3f±%.3f",
                label, category, length(plddts),
                isempty(plddts) ? 0.0 : mean(plddts),
                isempty(plddts) ? 0.0 : std(plddts)/sqrt(length(plddts)),
                isempty(tms) ? 0.0 : mean(tms),
                isempty(tms) ? 0.0 : std(tms)/sqrt(length(tms)))
        end
    end

    return plddt_data, tm_data
end

# ══════════════════════════════════════════════════════════════════════════════
# Collect ESMFold metrics (same as analyze_structure_results.jl)
# ══════════════════════════════════════════════════════════════════════════════

function collect_esmfold_metrics()
    plddt_data = Dict{String, Dict{String, Vector{Float64}}}()
    tm_data    = Dict{String, Dict{String, Vector{Float64}}}()

    esmfold_dirs = ["sa_generation", "stored"]
    esmfold_labels = ["SA generation", "Stored"]

    for fam in FAMILIES
        plddt_data[fam.name] = Dict{String, Vector{Float64}}()
        tm_data[fam.name]    = Dict{String, Vector{Float64}}()

        ref_pdb = joinpath(STRUCT_DATA_DIR, fam.id, "reference", "$(fam.pdb)_$(fam.chain).pdb")
        if !isfile(ref_pdb); continue; end

        for (method_dir, label) in zip(esmfold_dirs, esmfold_labels)
            struct_dir = joinpath(STRUCT_DATA_DIR, fam.id, "structures", method_dir)
            if !isdir(struct_dir)
                plddt_data[fam.name][label] = Float64[]
                tm_data[fam.name][label]    = Float64[]
                continue
            end

            pdb_files = filter(f -> endswith(f, ".pdb"), readdir(struct_dir, join=true))
            plddts = Float64[]
            tms = Float64[]

            for pdb_file in pdb_files
                pdb_content = read(pdb_file, String)
                plddt = let
                    ps = Float64[]
                    for line in split(pdb_content, '\n')
                        if startswith(line, "ATOM") && length(line) >= 66
                            try push!(ps, parse(Float64, strip(line[61:66]))); catch; end
                        end
                    end
                    isempty(ps) ? NaN : (let m=mean(ps); m < 1.5 ? m*100.0 : m; end)
                end
                !isnan(plddt) && push!(plddts, plddt)

                tm = run_tmalign(pdb_file, ref_pdb)
                !isnan(tm) && push!(tms, tm)
            end

            plddt_data[fam.name][label] = plddts
            tm_data[fam.name][label]    = tms
        end
    end

    return plddt_data, tm_data
end

# ══════════════════════════════════════════════════════════════════════════════
# Statistical tests for AF2 results
# ══════════════════════════════════════════════════════════════════════════════

function run_af2_statistical_tests(plddt_data, tm_data)
    @info "\n══════════════════════════════════════════════════════════"
    @info "  AlphaFold2 Statistical Tests: SA generation vs Stored"
    @info "══════════════════════════════════════════════════════════\n"

    test_rows = Vector{NamedTuple}()

    for fam in FAMILIES
        gen_plddt = get(plddt_data[fam.name], "SA generation", Float64[])
        sto_plddt = get(plddt_data[fam.name], "Stored", Float64[])
        gen_tm    = get(tm_data[fam.name], "SA generation", Float64[])
        sto_tm    = get(tm_data[fam.name], "Stored", Float64[])

        p_plddt = NaN; d_plddt = NaN
        p_tm = NaN; d_tm = NaN

        if length(gen_plddt) >= 5 && length(sto_plddt) >= 5
            p_plddt = pvalue(MannWhitneyUTest(gen_plddt, sto_plddt))
            d_plddt = cohens_d(gen_plddt, sto_plddt)
            @info @sprintf("  %s pLDDT: SA gen=%.1f±%.1f, Stored=%.1f±%.1f, p=%.4f (%s), d=%.2f",
                fam.name, mean(gen_plddt), std(gen_plddt)/sqrt(length(gen_plddt)),
                mean(sto_plddt), std(sto_plddt)/sqrt(length(sto_plddt)),
                p_plddt, pval_stars(p_plddt), d_plddt)
        end

        if length(gen_tm) >= 5 && length(sto_tm) >= 5
            p_tm = pvalue(MannWhitneyUTest(gen_tm, sto_tm))
            d_tm = cohens_d(gen_tm, sto_tm)
            @info @sprintf("  %s TM:    SA gen=%.3f±%.3f, Stored=%.3f±%.3f, p=%.4f (%s), d=%.2f",
                fam.name, mean(gen_tm), std(gen_tm)/sqrt(length(gen_tm)),
                mean(sto_tm), std(sto_tm)/sqrt(length(sto_tm)),
                p_tm, pval_stars(p_tm), d_tm)
        end

        push!(test_rows, (
            Family = fam.name, Predictor = "AlphaFold2",
            pLDDT_gen_mean = isempty(gen_plddt) ? NaN : mean(gen_plddt),
            pLDDT_sto_mean = isempty(sto_plddt) ? NaN : mean(sto_plddt),
            pLDDT_p = p_plddt, pLDDT_d = d_plddt, pLDDT_stars = pval_stars(p_plddt),
            TM_gen_mean = isempty(gen_tm) ? NaN : mean(gen_tm),
            TM_sto_mean = isempty(sto_tm) ? NaN : mean(sto_tm),
            TM_p = p_tm, TM_d = d_tm, TM_stars = pval_stars(p_tm),
        ))
    end

    return DataFrame(test_rows)
end

# ══════════════════════════════════════════════════════════════════════════════
# Combined ESMFold + AF2 figure (4-panel: ESMFold pLDDT, AF2 pLDDT,
#                                         ESMFold TM, AF2 TM)
# ══════════════════════════════════════════════════════════════════════════════

function generate_combined_figure(esm_plddt, esm_tm, af2_plddt, af2_tm)
    families = [fam.name for fam in FAMILIES]
    n_fam = length(families)

    colors_gen = RGB(0.20, 0.47, 0.69)   # SA gen - steel blue (canonical)
    colors_sto = RGB(0.20, 0.63, 0.17)   # Stored - green (canonical)

    bar_width = 0.3
    offsets = [-0.5, 0.5] .* bar_width
    labels = ["SA generation", "Stored"]

    function make_panel(data, ylabel, ylims_range, threshold, threshold_label, title_str)
        p = plot(size=(500, 350), ylabel=ylabel,
                 xticks=(1:n_fam, families), xrotation=30,
                 legend=:bottomright, legendfontsize=6,
                 ylims=ylims_range, grid=:y,
                 left_margin=5Plots.mm, bottom_margin=10Plots.mm,
                 top_margin=5Plots.mm, title=title_str, titlefontsize=9)

        if !isnan(threshold)
            hline!(p, [threshold], color=:gray, linestyle=:dash, linewidth=0.8,
                   label=threshold_label)
        end

        for (j, (label, color)) in enumerate(zip(labels, [colors_gen, colors_sto]))
            xs = Float64[]; ys = Float64[]; es = Float64[]
            for (i, fam) in enumerate(families)
                vals = get(data[fam], label, Float64[])
                if !isempty(vals)
                    push!(xs, i + offsets[j])
                    push!(ys, mean(vals))
                    push!(es, std(vals)/sqrt(length(vals)))
                end
            end
            bar!(p, xs, ys, yerror=es, bar_width=bar_width,
                 color=color, label=label, linewidth=0.5)
        end

        # Significance brackets
        bracket_dy = (ylims_range[2] - ylims_range[1]) * 0.015
        for (i, fam) in enumerate(families)
            gen_vals = get(data[fam], "SA generation", Float64[])
            sto_vals = get(data[fam], "Stored", Float64[])
            if length(gen_vals) >= 5 && length(sto_vals) >= 5
                pv = pvalue(MannWhitneyUTest(gen_vals, sto_vals))
                label_str = pval_stars(pv)
                max_y = max(
                    mean(gen_vals) + std(gen_vals)/sqrt(length(gen_vals)),
                    mean(sto_vals) + std(sto_vals)/sqrt(length(sto_vals))
                )
                bracket_y = max_y + bracket_dy * 2
                x1 = i + offsets[1]; x2 = i + offsets[2]
                draw_significance_bracket!(p, x1, x2, bracket_y, bracket_dy, label_str; fontsize=6)
            end
        end

        return p
    end

    p1 = make_panel(esm_plddt, "Mean pLDDT", (0, 115), 70.0, "pLDDT = 70",
                    "(A) ESMFold pLDDT")
    p2 = make_panel(af2_plddt, "Mean pLDDT", (0, 115), 70.0, "pLDDT = 70",
                    "(B) AlphaFold2 pLDDT")
    p3 = make_panel(esm_tm, "TM-score", (0, 1.15), 0.5, "TM = 0.5",
                    "(C) ESMFold TM-score")
    p4 = make_panel(af2_tm, "TM-score", (0, 1.15), 0.5, "TM = 0.5",
                    "(D) AlphaFold2 TM-score")

    composite = plot(p1, p2, p3, p4, layout=(2, 2), size=(1100, 800),
                     left_margin=5Plots.mm, bottom_margin=5Plots.mm)

    mkpath(FIGS_DIR)
    savefig(composite, joinpath(FIGS_DIR, "structure_validation_esmfold_af2.pdf"))
    @info "Saved: $(joinpath(FIGS_DIR, "structure_validation_esmfold_af2.pdf"))"

    # Copy to paper
    paper_figs = joinpath(_EXPERIMENT_ROOT, "..", "..", "paper", "figs")
    if isdir(paper_figs)
        cp(joinpath(FIGS_DIR, "structure_validation_esmfold_af2.pdf"),
           joinpath(paper_figs, "structure_validation_esmfold_af2.pdf"), force=true)
        @info "Copied to paper/figs/"
    end

    return composite
end

# ══════════════════════════════════════════════════════════════════════════════
# Correlation / concordance between ESMFold and AF2
# ══════════════════════════════════════════════════════════════════════════════

function generate_concordance_figure(esm_plddt, esm_tm, af2_plddt, af2_tm)
    families = [fam.name for fam in FAMILIES]

    # Collect per-family mean pLDDT and TM for both predictors
    esm_plddt_means = Float64[]
    af2_plddt_means = Float64[]
    esm_tm_means = Float64[]
    af2_tm_means = Float64[]
    fam_labels = String[]

    for fam in families
        for label in ["SA generation", "Stored"]
            ep = get(esm_plddt[fam], label, Float64[])
            ap = get(af2_plddt[fam], label, Float64[])
            et = get(esm_tm[fam], label, Float64[])
            at = get(af2_tm[fam], label, Float64[])

            if !isempty(ep) && !isempty(ap)
                push!(esm_plddt_means, mean(ep))
                push!(af2_plddt_means, mean(ap))
                push!(fam_labels, "$fam ($label)")
            end
            if !isempty(et) && !isempty(at)
                push!(esm_tm_means, mean(et))
                push!(af2_tm_means, mean(at))
            end
        end
    end

    p1 = scatter(esm_plddt_means, af2_plddt_means,
                 xlabel="ESMFold pLDDT", ylabel="AlphaFold2 pLDDT",
                 title="(A) pLDDT Concordance", titlefontsize=9,
                 legend=false, markersize=6, markercolor=RGB(0.20, 0.47, 0.69),
                 xlims=(40, 100), ylims=(40, 100), aspect_ratio=:equal, grid=true)
    plot!(p1, [40, 100], [40, 100], color=:gray, linestyle=:dash, linewidth=1, label=false)

    if length(esm_plddt_means) >= 3
        r = cor(esm_plddt_means, af2_plddt_means)
        annotate!(p1, [(50, 95, text(@sprintf("r = %.3f", r), 9, :left))])
    end

    p2 = scatter(esm_tm_means, af2_tm_means,
                 xlabel="ESMFold TM-score", ylabel="AlphaFold2 TM-score",
                 title="(B) TM-score Concordance", titlefontsize=9,
                 legend=false, markersize=6, markercolor=RGB(0.20, 0.63, 0.17),
                 xlims=(0.2, 1.0), ylims=(0.2, 1.0), aspect_ratio=:equal, grid=true)
    plot!(p2, [0.2, 1.0], [0.2, 1.0], color=:gray, linestyle=:dash, linewidth=1, label=false)

    if length(esm_tm_means) >= 3
        r = cor(esm_tm_means, af2_tm_means)
        annotate!(p2, [(0.3, 0.95, text(@sprintf("r = %.3f", r), 9, :left))])
    end

    composite = plot(p1, p2, layout=(1, 2), size=(900, 400),
                     left_margin=5Plots.mm, bottom_margin=8Plots.mm)

    savefig(composite, joinpath(FIGS_DIR, "esmfold_af2_concordance.pdf"))
    @info "Saved: $(joinpath(FIGS_DIR, "esmfold_af2_concordance.pdf"))"

    paper_figs = joinpath(_EXPERIMENT_ROOT, "..", "..", "paper", "figs")
    if isdir(paper_figs)
        cp(joinpath(FIGS_DIR, "esmfold_af2_concordance.pdf"),
           joinpath(paper_figs, "esmfold_af2_concordance.pdf"), force=true)
        @info "Copied to paper/figs/"
    end

    return composite
end

# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

function main()
    if !isdir(AF2_OUTPUT_DIR)
        @error """
        ColabFold output directory not found: $AF2_OUTPUT_DIR

        To generate AF2 predictions:
          1. Run: julia prepare_colabfold_input.jl
          2. Upload colabfold_input/ to Google Colab
          3. Run: python colabfold_predict.py
          4. Download colabfold_output/ back here
        """
        return
    end

    mkpath(RESULTS_DIR)
    mkpath(FIGS_DIR)

    # Collect AF2 metrics
    @info "Collecting AlphaFold2 metrics …"
    af2_plddt, af2_tm = collect_af2_metrics()

    # Statistical tests
    @info "\nRunning AF2 statistical tests …"
    df_af2 = run_af2_statistical_tests(af2_plddt, af2_tm)
    CSV.write(joinpath(RESULTS_DIR, "af2_statistical_tests.csv"), df_af2)
    @info "AF2 statistical tests saved."

    # Collect ESMFold metrics for comparison
    @info "\nCollecting ESMFold metrics for comparison …"
    esm_plddt, esm_tm = collect_esmfold_metrics()

    # Combined 4-panel figure
    @info "\nGenerating combined ESMFold + AF2 figure …"
    generate_combined_figure(esm_plddt, esm_tm, af2_plddt, af2_tm)

    # Concordance scatter plots
    @info "\nGenerating concordance figure …"
    generate_concordance_figure(esm_plddt, esm_tm, af2_plddt, af2_tm)

    @info "\nDone. Key outputs:"
    @info "  results/af2_statistical_tests.csv"
    @info "  figs/structure_validation_esmfold_af2.pdf"
    @info "  figs/esmfold_af2_concordance.pdf"
end

main()
