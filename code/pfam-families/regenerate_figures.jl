#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Regenerate cross-family figures — PNAS quality
# KL: pooled over all sequences (matching Table 2), SE via 1000-fold bootstrap
# Novelty: PCA-space cosine novelty = 1 - max_k cos(xi, m_k) (matching Table 2)
# SeqID: per-chain means, SE over 30 chains
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
    (id="PF00711", name="Defensin_beta"),
]

const n_chains = 30
const samples_per_chain = 5
const n_bootstrap_kl = 1000
const R_pca = 0.95

"""
    project_to_pca(seq::String, pca_model, L::Int) -> Vector{Float64}

One-hot encode a single sequence, project into PCA space, and unit-normalize.
"""
function project_to_pca(seq::String, pca_model, L::Int)
    AA_ALPHABET = collect("ACDEFGHIKLMNPQRSTVWY")
    n_aa = length(AA_ALPHABET)
    aa_to_idx = Dict(aa => i for (i, aa) in enumerate(AA_ALPHABET))

    oh = zeros(Float64, n_aa * L)
    for (pos, ch) in enumerate(seq)
        pos > L && break
        idx = get(aa_to_idx, ch, 0)
        idx > 0 && (oh[(pos-1)*n_aa + idx] = 1.0)
    end

    v = MultivariateStats.transform(pca_model, oh)
    nv = norm(v)
    nv > 1e-12 ? v ./ nv : v
end

@info "Computing metrics …"

# collect per-family, per-method results with SE
rows = Vector{NamedTuple}()

for fam in FAMILIES
    @info "  $(fam.name) …"
    fam_dir = joinpath(PFAM_DATA_DIR, fam.id)

    # load stored sequences
    stored_data = parse_fasta(joinpath(fam_dir, "stored_sequences.fasta"))
    stored_seqs = [s[2] for s in stored_data]

    # rebuild memory matrix for PCA-space novelty (matching Table 2)
    char_mat = permutedims(hcat([collect(s) for s in stored_seqs]...))
    X̂, pca_model, L_out, d_full = build_memory_matrix(char_mat; pratio=R_pca)
    L = size(char_mat, 2)
    d, K = size(X̂)

    # load SA retrieval and generation sequences
    for (label, method_name) in [("sa_retrieval", "SA retrieval"),
                                  ("sa_generation", "SA generation")]
        fasta_file = joinpath(fam_dir, "$(label)_sequences.fasta")
        isfile(fasta_file) || continue
        gen_data = parse_fasta(fasta_file)
        gen_seqs = [s[2] for s in gen_data]

        # ── KL: pooled over ALL sequences, SE via bootstrap ──
        kl_pooled = aa_composition_kl(gen_seqs, stored_seqs)
        n_gen = length(gen_seqs)
        boot_kls = Float64[]
        for _ in 1:n_bootstrap_kl
            idx = rand(1:n_gen, n_gen)
            push!(boot_kls, aa_composition_kl(gen_seqs[idx], stored_seqs))
        end
        kl_se = std(boot_kls)

        # ── PCA-space novelty and AA-space SeqID: per-chain means ──
        gen_pca = [project_to_pca(s, pca_model, L) for s in gen_seqs]
        n_total = length(gen_seqs)
        spc = min(samples_per_chain, n_total ÷ n_chains)
        chain_novelties = Float64[]
        chain_seqids = Float64[]

        for c in 1:n_chains
            si = (c - 1) * spc + 1
            ei = c * spc
            ei > n_total && break
            # novelty: 1 - max cosine similarity in PCA space
            push!(chain_novelties, mean(1.0 - nearest_cosine_similarity(gen_pca[j], X̂) for j in si:ei))
            # seqid: nearest sequence identity in AA space
            push!(chain_seqids, mean(nearest_sequence_identity(gen_seqs[j], stored_seqs) for j in si:ei))
        end

        n_ch = length(chain_novelties)
        push!(rows, (
            Family     = fam.name,
            Method     = method_name,
            KL_mean    = kl_pooled,
            KL_se      = kl_se,
            Nov_mean   = mean(chain_novelties),
            Nov_se     = std(chain_novelties) / sqrt(n_ch),
            SID_mean   = mean(chain_seqids),
            SID_se     = std(chain_seqids) / sqrt(n_ch),
        ))
    end

    # Bootstrap baseline: resample stored patterns in PCA space
    Random.seed!(12345)
    bs_pca = [copy(X̂[:, rand(1:K)]) for _ in 1:150]
    bs_seqs = [decode_sample(ξ, pca_model, L) for ξ in bs_pca]

    # ── KL: pooled, SE via bootstrap ──
    kl_pooled = aa_composition_kl(bs_seqs, stored_seqs)
    n_gen = length(bs_seqs)
    boot_kls = Float64[]
    for _ in 1:n_bootstrap_kl
        idx = rand(1:n_gen, n_gen)
        push!(boot_kls, aa_composition_kl(bs_seqs[idx], stored_seqs))
    end
    kl_se = std(boot_kls)

    # ── PCA-space novelty and AA-space SeqID: per-chain ──
    chain_novelties = Float64[]
    chain_seqids = Float64[]
    spc = 5
    for c in 1:n_chains
        si = (c - 1) * spc + 1
        ei = c * spc
        ei > length(bs_pca) && break
        push!(chain_novelties, mean(1.0 - nearest_cosine_similarity(bs_pca[j], X̂) for j in si:ei))
        push!(chain_seqids, mean(nearest_sequence_identity(bs_seqs[j], stored_seqs) for j in si:ei))
    end

    n_ch = length(chain_novelties)
    push!(rows, (
        Family     = fam.name,
        Method     = "Bootstrap",
        KL_mean    = kl_pooled,
        KL_se      = kl_se,
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
    ylabel="Novelty (1 - max cosine sim.)", title="(B)  Sample novelty",
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
