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
# PNAS-quality figures
# ══════════════════════════════════════════════════════════════════════════════

@info "Generating PNAS-quality figures …"

family_names = [f.name for f in FAMILIES]
method_order = ["SA generation", "SA retrieval", "Bootstrap"]

method_colors = Dict(
    "SA generation" => RGB(0.20, 0.47, 0.69),   # steel blue
    "SA retrieval"  => RGB(0.89, 0.44, 0.32),   # coral
    "Bootstrap"     => RGB(0.50, 0.50, 0.50),   # gray
)

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
    size           = (600, 380),
    dpi            = 300,
    margin         = 4Plots.mm,
    bottom_margin  = 6Plots.mm,
    left_margin    = 5Plots.mm,
)

n_methods = length(method_order)
bar_width = 0.25
offsets = [-(n_methods-1)/2 * bar_width + (i-1) * bar_width for i in 1:n_methods]

# ── Cross-family KL ──────────────────────────────────────────────────────────

p_kl = plot(;
    xlabel = "",
    ylabel = "KL divergence (AA composition)",
    title  = "Amino acid composition fidelity across families",
    legend = :topright,
    ylims  = (0, Inf),
    pnas_defaults...
)

for (mi, m) in enumerate(method_order)
    sub = filter(r -> r.Method == m, df)
    kl_vals = Float64[]
    kl_errs = Float64[]
    for fam_name in family_names
        row = filter(r -> r.Family == fam_name, sub)
        if nrow(row) > 0
            push!(kl_vals, row.KL_mean[1])
            push!(kl_errs, row.KL_se[1])
        else
            push!(kl_vals, NaN)
            push!(kl_errs, 0.0)
        end
    end

    x_pos = collect(1:length(family_names)) .+ offsets[mi]
    bar!(p_kl, x_pos, kl_vals;
         bar_width = bar_width,
         yerror    = kl_errs,
         label     = m,
         color     = method_colors[m],
         alpha     = 0.85,
         lw        = 0.5,
         linecolor = :gray30,
    )
end
plot!(p_kl, xticks=(1:length(family_names), family_names), xrotation=0)
savefig(p_kl, joinpath(PFAM_FIG_DIR, "cross_family_kl.pdf"))
@info "  Saved cross_family_kl.pdf"

# ── Cross-family Novelty ─────────────────────────────────────────────────────

p_nov = plot(;
    xlabel = "",
    ylabel = "Novelty (1 − max seq. identity)",
    title  = "Sample novelty across families",
    legend = :topleft,
    ylims  = (0, Inf),
    pnas_defaults...
)

for (mi, m) in enumerate(method_order)
    sub = filter(r -> r.Method == m, df)
    nov_vals = Float64[]
    nov_errs = Float64[]
    for fam_name in family_names
        row = filter(r -> r.Family == fam_name, sub)
        if nrow(row) > 0
            push!(nov_vals, row.Nov_mean[1])
            push!(nov_errs, row.Nov_se[1])
        else
            push!(nov_vals, NaN)
            push!(nov_errs, 0.0)
        end
    end

    x_pos = collect(1:length(family_names)) .+ offsets[mi]
    bar!(p_nov, x_pos, nov_vals;
         bar_width = bar_width,
         yerror    = nov_errs,
         label     = m,
         color     = method_colors[m],
         alpha     = 0.85,
         lw        = 0.5,
         linecolor = :gray30,
    )
end
plot!(p_nov, xticks=(1:length(family_names), family_names), xrotation=0)
savefig(p_nov, joinpath(PFAM_FIG_DIR, "cross_family_novelty.pdf"))
@info "  Saved cross_family_novelty.pdf"

# ── Cross-family SeqID ───────────────────────────────────────────────────────

p_sid = plot(;
    xlabel = "",
    ylabel = "Nearest sequence identity",
    title  = "Sequence identity to nearest stored pattern",
    legend = :topright,
    ylims  = (0, 1.05),
    pnas_defaults...
)

for (mi, m) in enumerate(method_order)
    sub = filter(r -> r.Method == m, df)
    sid_vals = Float64[]
    sid_errs = Float64[]
    for fam_name in family_names
        row = filter(r -> r.Family == fam_name, sub)
        if nrow(row) > 0
            push!(sid_vals, row.SID_mean[1])
            push!(sid_errs, row.SID_se[1])
        else
            push!(sid_vals, NaN)
            push!(sid_errs, 0.0)
        end
    end

    x_pos = collect(1:length(family_names)) .+ offsets[mi]
    bar!(p_sid, x_pos, sid_vals;
         bar_width = bar_width,
         yerror    = sid_errs,
         label     = m,
         color     = method_colors[m],
         alpha     = 0.85,
         lw        = 0.5,
         linecolor = :gray30,
    )
end
plot!(p_sid, xticks=(1:length(family_names), family_names), xrotation=0)
savefig(p_sid, joinpath(PFAM_FIG_DIR, "cross_family_seqid.pdf"))
@info "  Saved cross_family_seqid.pdf"

@info "Done."
