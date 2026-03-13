#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Evaluate all baselines: read generated FASTAs, compute metrics, compare to SA
#
# This script reads generated sequences from each baseline method (HMM, EvoDiff,
# MSA Transformer) and the existing SA results, then produces a unified comparison
# table and composite figure.
#
# Run after: run_hmm_baseline.jl (and optionally run_evodiff_baseline.py,
#            run_msa_transformer_baseline.py)
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment …"
const _EXPERIMENT_ROOT = @__DIR__
include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))
using LaTeXStrings
using Printf

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

# Baseline FASTA file names and method labels
const BASELINE_FILES = [
    ("hmm_sequences.fasta",     "HMM emit"),
    ("evodiff_sequences.fasta", "EvoDiff"),
    ("msat_sequences.fasta",    "MSA Transformer"),
]

# SA results from the pfam-families experiment (we reload and re-evaluate
# so everything uses the same evaluation code)
const SA_FILES = [
    ("sa_generation_sequences.fasta", "SA generation"),
    ("sa_retrieval_sequences.fasta",  "SA retrieval"),
]

const n_chains_eval = 30
const spc_eval = 5

# ══════════════════════════════════════════════════════════════════════════════
# Evaluation helpers
# ══════════════════════════════════════════════════════════════════════════════

"""
Evaluate a set of generated sequences against stored sequences.
Returns a NamedTuple with mean and SE for KL, novelty, and SeqID.
"""
function evaluate_sequences(gen_seqs::Vector{String}, stored_seqs::Vector{String},
                             X̂::Matrix{Float64}, pca_model, L::Int)

    # Encode generated sequences into PCA space for novelty
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

    # Full-sample point estimates
    kl_full = aa_composition_kl(gen_seqs, stored_seqs)
    nov_full = mean(sample_novelty(ξ, X̂) for ξ in gen_pca)
    sid_full = mean(nearest_sequence_identity(s, stored_seqs) for s in gen_seqs)

    # Per-chain SE
    n_total = length(gen_seqs)
    spc = min(spc_eval, n_total ÷ n_chains_eval)
    spc < 1 && return (KL_mean=kl_full, KL_se=0.0, Nov_mean=nov_full, Nov_se=0.0,
                        SID_mean=sid_full, SID_se=0.0, ValidAA=mean(valid_residue_fraction(s) for s in gen_seqs))

    chain_kls = Float64[]
    chain_novs = Float64[]
    chain_sids = Float64[]

    for c in 1:n_chains_eval
        i1 = (c-1)*spc + 1
        i2 = c*spc
        i2 > n_total && break
        ch_seqs = gen_seqs[i1:i2]
        ch_pca  = gen_pca[i1:i2]
        push!(chain_kls, aa_composition_kl(ch_seqs, stored_seqs))
        push!(chain_novs, mean(sample_novelty(ξ, X̂) for ξ in ch_pca))
        push!(chain_sids, mean(nearest_sequence_identity(s, stored_seqs) for s in ch_seqs))
    end

    n_ch = length(chain_kls)
    return (
        KL_mean  = kl_full,
        KL_se    = n_ch > 1 ? std(chain_kls)  / sqrt(n_ch) : 0.0,
        Nov_mean = nov_full,
        Nov_se   = n_ch > 1 ? std(chain_novs) / sqrt(n_ch) : 0.0,
        SID_mean = sid_full,
        SID_se   = n_ch > 1 ? std(chain_sids) / sqrt(n_ch) : 0.0,
        ValidAA  = mean(valid_residue_fraction(s) for s in gen_seqs),
    )
end

"""
Read a FASTA file and clean sequences to length L (truncate/pad, keep only standard AAs).
"""
function load_and_clean_seqs(fasta_path::String, L::Int)
    data = parse_fasta(fasta_path)
    seqs = String[]
    for (_, raw_seq) in data
        cleaned = filter(c -> c in AA_ALPHABET, uppercase(raw_seq))
        if length(cleaned) >= L
            push!(seqs, cleaned[1:L])
        elseif length(cleaned) > 0
            # pad with most common residue placeholder
            push!(seqs, cleaned * repeat("A", L - length(cleaned)))
        end
    end
    return seqs
end

# ══════════════════════════════════════════════════════════════════════════════
# Main evaluation
# ══════════════════════════════════════════════════════════════════════════════

function evaluate_all_baselines()
    rows = Vector{NamedTuple}()

    for fam in FAMILIES
        @info "════════════════════════════════════════════════════════════"
        @info "  $(fam.name) ($(fam.id))"
        @info "════════════════════════════════════════════════════════════"

        # Load stored sequences
        stored_fasta = joinpath(PFAM_DATA_DIR, fam.id, "stored_sequences.fasta")
        if !isfile(stored_fasta)
            @warn "  Missing stored sequences — skipping"
            continue
        end
        stored_data = parse_fasta(stored_fasta)
        stored_seqs = [s[2] for s in stored_data]

        # Build memory matrix for PCA-space novelty
        sto_file = joinpath(PFAM_DATA_DIR, fam.id, "$(fam.id)_seed.sto")
        raw_seqs = parse_stockholm(sto_file)
        char_mat, _ = clean_alignment(raw_seqs)
        K, L = size(char_mat)
        X̂, pca_model, L, d_full = build_memory_matrix(char_mat; pratio=0.95)

        # ── Evaluate SA methods (from pfam-families data) ────────────────
        for (fname, method_name) in SA_FILES
            fpath = joinpath(PFAM_DATA_DIR, fam.id, fname)
            if !isfile(fpath)
                @info "  $method_name: not found, skipping"
                continue
            end
            sa_data = parse_fasta(fpath)
            sa_seqs = [s[2] for s in sa_data]
            # SA sequences are already the right length (decoded from PCA)
            metrics = evaluate_sequences(sa_seqs, stored_seqs, X̂, pca_model, L)
            @info "  $method_name: KL=$(round(metrics.KL_mean, digits=4)), Nov=$(round(metrics.Nov_mean, digits=4)), SeqID=$(round(metrics.SID_mean, digits=4))"
            push!(rows, (Family=fam.name, PfamID=fam.id, Method=method_name, metrics...))
        end

        # ── Evaluate bootstrap baseline ──────────────────────────────────
        # Regenerate bootstrap (same as pfam-families experiment)
        Random.seed!(12345)
        bs_seqs = [stored_seqs[rand(1:length(stored_seqs))] for _ in 1:150]
        metrics = evaluate_sequences(bs_seqs, stored_seqs, X̂, pca_model, L)
        @info "  Bootstrap: KL=$(round(metrics.KL_mean, digits=4)), Nov=$(round(metrics.Nov_mean, digits=4)), SeqID=$(round(metrics.SID_mean, digits=4))"
        push!(rows, (Family=fam.name, PfamID=fam.id, Method="Bootstrap", metrics...))

        # ── Evaluate each baseline method ────────────────────────────────
        for (fname, method_name) in BASELINE_FILES
            fpath = joinpath(BASELINE_DATA_DIR, fam.id, fname)
            if !isfile(fpath)
                @info "  $method_name: not found (run the corresponding script first)"
                continue
            end
            gen_seqs = load_and_clean_seqs(fpath, L)
            if isempty(gen_seqs)
                @warn "  $method_name: no valid sequences after cleaning"
                continue
            end
            metrics = evaluate_sequences(gen_seqs, stored_seqs, X̂, pca_model, L)
            @info "  $method_name: KL=$(round(metrics.KL_mean, digits=4)), Nov=$(round(metrics.Nov_mean, digits=4)), SeqID=$(round(metrics.SID_mean, digits=4))"
            push!(rows, (Family=fam.name, PfamID=fam.id, Method=method_name, metrics...))
        end
    end

    df = DataFrame(rows)
    CSV.write(joinpath(RESULTS_DIR, "baseline_comparison.csv"), df)
    @info "Results saved to results/baseline_comparison.csv"

    # Print summary
    pretty_table(df[:, [:Family, :Method, :KL_mean, :KL_se, :Nov_mean, :Nov_se, :SID_mean, :SID_se]])

    return df
end

# ══════════════════════════════════════════════════════════════════════════════
# Figures: PNAS-quality comparison
# ══════════════════════════════════════════════════════════════════════════════

function generate_comparison_figures(df::DataFrame)
    family_names = [f.name for f in FAMILIES]

    # Only plot methods that have data
    available_methods = unique(df.Method)
    method_order_full = ["SA generation", "SA retrieval", "HMM emit",
                         "EvoDiff", "MSA Transformer", "Bootstrap"]
    method_order = [m for m in method_order_full if m in available_methods]

    method_colors = Dict(
        "SA generation"    => RGB(0.20, 0.47, 0.69),   # steel blue
        "SA retrieval"     => RGB(0.89, 0.44, 0.32),   # coral
        "HMM emit"         => RGB(0.30, 0.69, 0.29),   # green
        "EvoDiff"          => RGB(0.93, 0.68, 0.20),   # amber
        "MSA Transformer"  => RGB(0.58, 0.32, 0.68),   # purple
        "Bootstrap"        => RGB(0.50, 0.50, 0.50),   # gray
    )

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

    n_methods = length(method_order)
    bar_width = 0.8 / n_methods
    offsets = [-(n_methods-1)/2 * bar_width + (i-1) * bar_width for i in 1:n_methods]

    function make_bar_panel(df, family_names, method_order, method_colors, offsets, bar_width;
                            val_col::Symbol, err_col::Symbol, ylabel, title,
                            legend=:none, ylims=(0, Inf))
        p = plot(; xlabel="", ylabel=ylabel, title=title, legend=legend, ylims=ylims, panel_defaults...)
        for (mi, m) in enumerate(method_order)
            sub = filter(r -> r.Method == m, df)
            vals = Float64[]
            errs = Float64[]
            for fn in family_names
                row = filter(r -> r.Family == fn, sub)
                if nrow(row) > 0
                    push!(vals, getproperty(row, val_col)[1])
                    push!(errs, getproperty(row, err_col)[1])
                else
                    push!(vals, NaN)
                    push!(errs, 0.0)
                end
            end
            x_pos = collect(1:length(family_names)) .+ offsets[mi]
            bar!(p, x_pos, vals; bar_width=bar_width, yerror=errs,
                 label=m, color=method_colors[m], alpha=0.85, lw=0.5, linecolor=:gray30)
        end
        plot!(p, xticks=(1:length(family_names), family_names), xrotation=0)
        return p
    end

    # (A) KL divergence
    p_kl = make_bar_panel(df, family_names, method_order, method_colors, offsets, bar_width;
        val_col=:KL_mean, err_col=:KL_se,
        ylabel=L"\mathrm{KL\ divergence}", title="(A)  AA composition fidelity",
        legend=:topright)

    # (B) Novelty
    p_nov = make_bar_panel(df, family_names, method_order, method_colors, offsets, bar_width;
        val_col=:Nov_mean, err_col=:Nov_se,
        ylabel=L"\mathrm{Novelty}", title="(B)  Sample novelty",
        legend=:none)

    # (C) Sequence identity
    p_sid = make_bar_panel(df, family_names, method_order, method_colors, offsets, bar_width;
        val_col=:SID_mean, err_col=:SID_se,
        ylabel=L"\mathrm{Nearest\ seq.\ identity}", title="(C)  Sequence identity",
        legend=:none, ylims=(0, 1.05))

    # Composite: 3 rows × 1 column
    p_composite = plot(p_kl, p_nov, p_sid;
        layout = (3, 1),
        size   = (800, 1000),
        dpi    = 300,
    )

    savefig(p_composite, joinpath(FIG_DIR, "baseline_comparison_composite.pdf"))
    @info "  Saved baseline_comparison_composite.pdf"

    # Individual panels
    savefig(p_kl,  joinpath(FIG_DIR, "baseline_kl.pdf"))
    savefig(p_nov, joinpath(FIG_DIR, "baseline_novelty.pdf"))
    savefig(p_sid, joinpath(FIG_DIR, "baseline_seqid.pdf"))

    return p_composite
end

# ══════════════════════════════════════════════════════════════════════════════
# LaTeX table output for paper
# ══════════════════════════════════════════════════════════════════════════════

function generate_latex_table(df::DataFrame)
    available_methods = unique(df.Method)
    method_order_full = ["SA generation", "SA retrieval", "HMM emit",
                         "EvoDiff", "MSA Transformer", "Bootstrap"]
    method_order = [m for m in method_order_full if m in available_methods]
    method_abbrev = Dict(
        "SA generation"   => "SA gen",
        "SA retrieval"    => "SA ret",
        "HMM emit"        => "HMM",
        "EvoDiff"         => "EvoDiff",
        "MSA Transformer" => "MSA-T",
        "Bootstrap"       => "BS",
    )

    println("\n% ── LaTeX table rows ──")
    for fam in FAMILIES
        fam_df = filter(r -> r.Family == fam.name, df)
        isempty(fam_df) && continue
        for (mi, m) in enumerate(method_order)
            row = filter(r -> r.Method == m, fam_df)
            nrow(row) == 0 && continue
            r = row[1, :]
            prefix = mi == 1 ? "\\multirow{$(length(method_order))}{*}{$(fam.name)}" : ""
            kl_str  = @sprintf("\$%.3f \\pm %.3f\$", r.KL_mean, r.KL_se)
            nov_str = @sprintf("\$%.3f \\pm %.3f\$", r.Nov_mean, r.Nov_se)
            sid_str = @sprintf("\$%.3f \\pm %.3f\$", r.SID_mean, r.SID_se)
            println("$prefix & $(method_abbrev[m]) & $kl_str & $nov_str & $sid_str \\\\")
        end
        println("\\midrule")
    end
end

# ── Run ──────────────────────────────────────────────────────────────────────

df = evaluate_all_baselines()
generate_comparison_figures(df)
generate_latex_table(df)

@info "Done — all baseline evaluations complete."
