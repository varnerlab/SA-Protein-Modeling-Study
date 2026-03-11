#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Regenerate cross-family figures — PNAS quality with per-chain SE error bars
# Uses per-chain bootstrap: 30 chains × 5 samples → SE over chain-level metrics
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment …"
const _EXPERIMENT_ROOT = @__DIR__
include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))

const PFAM_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "data")
const PFAM_FIG_DIR  = joinpath(_EXPERIMENT_ROOT, "figs")

const FAMILIES = [
    (id="PF00076", name="RRM"),
    (id="PF00018", name="SH3"),
    (id="PF00397", name="WW"),
    (id="PF00014", name="Kunitz"),
    (id="PF00096", name="zf-C2H2"),
    (id="PF00595", name="PDZ"),
    (id="PF00069", name="Pkinase"),
]

const n_chains = 30
const samples_per_chain = 5
const R_pca = 0.95

@info "Computing per-chain error bars …"

# collect per-family, per-method results with SE
rows = Vector{NamedTuple}()

for fam in FAMILIES
    @info "  $(fam.name) …"
    fam_dir = joinpath(PFAM_DATA_DIR, fam.id)

    # load stored sequences
    stored_data = parse_fasta(joinpath(fam_dir, "stored_sequences.fasta"))
    stored_seqs = [s[2] for s in stored_data]

    # load SA retrieval and generation sequences
    for (label, method_name) in [("sa_retrieval", "SA retrieval"),
                                  ("sa_generation", "SA generation")]
        fasta_file = joinpath(fam_dir, "$(label)_sequences.fasta")
        isfile(fasta_file) || continue
        gen_data = parse_fasta(fasta_file)
        gen_seqs = [s[2] for s in gen_data]

        # split into 30 chains of 5 samples
        n_total = length(gen_seqs)
        spc = min(samples_per_chain, n_total ÷ n_chains)
        chain_kls = Float64[]
        chain_novelties = Float64[]
        chain_seqids = Float64[]

        for c in 1:n_chains
            start_idx = (c - 1) * spc + 1
            end_idx = c * spc
            end_idx > n_total && break
            chain_seqs = gen_seqs[start_idx:end_idx]

            push!(chain_kls, aa_composition_kl(chain_seqs, stored_seqs))
            push!(chain_seqids, mean(nearest_sequence_identity(s, stored_seqs) for s in chain_seqs))
            push!(chain_novelties, mean(1.0 - nearest_sequence_identity(s, stored_seqs) for s in chain_seqs))
        end

        n_ch = length(chain_kls)
        push!(rows, (
            Family     = fam.name,
            Method     = method_name,
            KL_mean    = mean(chain_kls),
            KL_se      = std(chain_kls) / sqrt(n_ch),
            Nov_mean   = mean(chain_novelties),
            Nov_se     = std(chain_novelties) / sqrt(n_ch),
            SID_mean   = mean(chain_seqids),
            SID_se     = std(chain_seqids) / sqrt(n_ch),
        ))
    end

    # Bootstrap baseline: resample stored sequences into 30 groups of 5
    gen_seqs = begin
        Random.seed!(12345)
        [stored_seqs[rand(1:length(stored_seqs))] for _ in 1:150]
    end
    chain_kls = Float64[]
    chain_seqids = Float64[]
    chain_novelties = Float64[]
    spc = 5

    for c in 1:n_chains
        start_idx = (c - 1) * spc + 1
        end_idx = c * spc
        chain_seqs = gen_seqs[start_idx:end_idx]
        push!(chain_kls, aa_composition_kl(chain_seqs, stored_seqs))
        push!(chain_seqids, mean(nearest_sequence_identity(s, stored_seqs) for s in chain_seqs))
        push!(chain_novelties, mean(1.0 - nearest_sequence_identity(s, stored_seqs) for s in chain_seqs))
    end

    n_ch = length(chain_kls)
    push!(rows, (
        Family     = fam.name,
        Method     = "Bootstrap",
        KL_mean    = mean(chain_kls),
        KL_se      = std(chain_kls) / sqrt(n_ch),
        Nov_mean   = mean(chain_novelties),
        Nov_se     = std(chain_novelties) / sqrt(n_ch),
        SID_mean   = mean(chain_seqids),
        SID_se     = std(chain_seqids) / sqrt(n_ch),
    ))
end

df = DataFrame(rows)
CSV.write(joinpath(PFAM_DATA_DIR, "pfam_results_with_se.csv"), df)

# ══════════════════════════════════════════════════════════════════════════════
# PNAS-quality composite figure: 3-panel (KL, Novelty, SeqID)
# ══════════════════════════════════════════════════════════════════════════════

@info "Generating PNAS-quality composite figure …"

family_names = [f.name for f in FAMILIES]
method_order = ["SA generation", "SA retrieval", "Bootstrap"]

method_colors = Dict(
    "SA generation" => RGB(0.20, 0.47, 0.69),   # steel blue
    "SA retrieval"  => RGB(0.89, 0.44, 0.32),   # coral
    "Bootstrap"     => RGB(0.50, 0.50, 0.50),   # gray
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
bar_width = 0.25
offsets = [-(n_methods-1)/2 * bar_width + (i-1) * bar_width for i in 1:n_methods]

# helper: build a grouped bar panel from df
function make_bar_panel(df, family_names, method_order, method_colors, offsets, bar_width;
                        val_col::Symbol, err_col::Symbol, ylabel::String, title::String,
                        legend=:none, ylims=(0, Inf), panel_defaults...)
    p = plot(; xlabel="", ylabel=ylabel, title=title, legend=legend, ylims=ylims, panel_defaults...)
    for (mi, m) in enumerate(method_order)
        sub = filter(r -> r.Method == m, df)
        vals = Float64[]
        errs = Float64[]
        for fam_name in family_names
            row = filter(r -> r.Family == fam_name, sub)
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

p_kl = make_bar_panel(df, family_names, method_order, method_colors, offsets, bar_width;
    val_col=:KL_mean, err_col=:KL_se,
    ylabel="KL divergence", title="(A)  AA composition fidelity",
    legend=:topright, panel_defaults...)

p_nov = make_bar_panel(df, family_names, method_order, method_colors, offsets, bar_width;
    val_col=:Nov_mean, err_col=:Nov_se,
    ylabel="Novelty (1 − max seq. identity)", title="(B)  Sample novelty",
    legend=:none, panel_defaults...)

p_sid = make_bar_panel(df, family_names, method_order, method_colors, offsets, bar_width;
    val_col=:SID_mean, err_col=:SID_se,
    ylabel="Nearest sequence identity", title="(C)  Sequence identity",
    legend=:none, ylims=(0, 1.05), panel_defaults...)

# composite: 3 rows × 1 column
p_composite = plot(p_kl, p_nov, p_sid;
    layout = (3, 1),
    size   = (700, 900),
    dpi    = 300,
)

savefig(p_composite, joinpath(PFAM_FIG_DIR, "cross_family_composite.pdf"))
@info "  Saved cross_family_composite.pdf"

# also save individual panels for flexibility
savefig(p_kl,  joinpath(PFAM_FIG_DIR, "cross_family_kl.pdf"))
savefig(p_nov, joinpath(PFAM_FIG_DIR, "cross_family_novelty.pdf"))
savefig(p_sid, joinpath(PFAM_FIG_DIR, "cross_family_seqid.pdf"))

@info "Done."
