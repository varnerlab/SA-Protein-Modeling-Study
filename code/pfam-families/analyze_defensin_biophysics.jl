#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Biophysical property analysis for Defensin_beta (PF00711)
#
# Computes AMP-relevant physicochemical properties for stored vs. generated
# sequences: net charge, hydrophobic ratio, cysteine count, mean hydrophobicity.
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment …"
const _EXPERIMENT_ROOT = @__DIR__
include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))
using LaTeXStrings, Printf

const DATA_DIR = joinpath(_EXPERIMENT_ROOT, "data", "PF00711")
const BASELINE_DIR = joinpath(_EXPERIMENT_ROOT, "..", "baselines", "data", "PF00711")
const FIG_DIR = joinpath(_EXPERIMENT_ROOT, "figs")
mkpath(FIG_DIR)

# ══════════════════════════════════════════════════════════════════════════════
# Biophysical property functions
# ══════════════════════════════════════════════════════════════════════════════

# Kyte-Doolittle hydrophobicity scale
const KD_HYDRO = Dict(
    'A'=> 1.8, 'R'=>-4.5, 'N'=>-3.5, 'D'=>-3.5, 'C'=> 2.5,
    'Q'=>-3.5, 'E'=>-3.5, 'G'=>-0.4, 'H'=>-3.2, 'I'=> 4.5,
    'L'=> 3.8, 'K'=>-3.9, 'M'=> 1.9, 'F'=> 2.8, 'P'=>-1.6,
    'S'=>-0.8, 'T'=>-0.7, 'W'=>-0.9, 'Y'=>-1.3, 'V'=> 4.2,
)

const HYDROPHOBIC_SET = Set(['A','V','I','L','M','F','W','P'])

"""Net charge at pH 7 (simplified: K,R = +1; D,E = -1; H = +0.5)."""
function net_charge(seq::String)
    q = 0.0
    for c in seq
        if c == 'K' || c == 'R'
            q += 1.0
        elseif c == 'D' || c == 'E'
            q -= 1.0
        elseif c == 'H'
            q += 0.5
        end
    end
    return q
end

"""Fraction of residues that are hydrophobic."""
function hydrophobic_ratio(seq::String)
    n_hydro = count(c -> c in HYDROPHOBIC_SET, seq)
    return n_hydro / length(seq)
end

"""Mean Kyte-Doolittle hydrophobicity per residue."""
function mean_hydrophobicity(seq::String)
    vals = [get(KD_HYDRO, c, 0.0) for c in seq]
    return mean(vals)
end

"""Count of cysteine residues (critical for defensin disulfide scaffold)."""
function cysteine_count(seq::String)
    return count(c -> c == 'C', seq)
end

"""Compute all biophysical properties for a sequence."""
function biophysical_properties(seq::String)
    return (
        charge = net_charge(seq),
        hydro_ratio = hydrophobic_ratio(seq),
        mean_hydro = mean_hydrophobicity(seq),
        n_cys = cysteine_count(seq),
        length = length(seq),
    )
end

# ══════════════════════════════════════════════════════════════════════════════
# Load sequences
# ══════════════════════════════════════════════════════════════════════════════

function load_seqs(fasta_path::String)
    data = parse_fasta(fasta_path)
    return [s[2] for s in data]
end

function load_and_clean(fasta_path::String, L::Int)
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

@info "Loading sequences …"

stored_seqs = load_seqs(joinpath(DATA_DIR, "stored_sequences.fasta"))
sa_gen_seqs = load_seqs(joinpath(DATA_DIR, "sa_generation_sequences.fasta"))
L = length(stored_seqs[1])

# Load baselines if available
methods = Dict{String, Vector{String}}()
methods["Stored"] = stored_seqs
methods["SA generation"] = sa_gen_seqs

hmm_fasta = joinpath(BASELINE_DIR, "hmm_sequences.fasta")
if isfile(hmm_fasta)
    methods["HMM emit"] = load_and_clean(hmm_fasta, L)
end

evodiff_fasta = joinpath(BASELINE_DIR, "evodiff_sequences.fasta")
if isfile(evodiff_fasta)
    methods["EvoDiff"] = load_and_clean(evodiff_fasta, L)
end

msat_fasta = joinpath(BASELINE_DIR, "msat_sequences.fasta")
if isfile(msat_fasta)
    methods["MSA Transformer"] = load_and_clean(msat_fasta, L)
end

# ══════════════════════════════════════════════════════════════════════════════
# Compute properties for all methods
# ══════════════════════════════════════════════════════════════════════════════

@info "Computing biophysical properties …"

results = Dict{String, DataFrame}()
summary_rows = Vector{NamedTuple}()

for (method_name, seqs) in methods
    props = [biophysical_properties(s) for s in seqs]
    df = DataFrame(props)
    df.method .= method_name
    results[method_name] = df

    # Summary statistics
    charges = [p.charge for p in props]
    hydro_ratios = [p.hydro_ratio for p in props]
    n_cys_vals = [p.n_cys for p in props]
    mean_hydros = [p.mean_hydro for p in props]

    # AMP plausibility filter: charge >= +2, hydro_ratio in [0.3, 0.7], n_cys >= 4
    n_pass = count(p -> p.charge >= 2.0 && 0.3 <= p.hydro_ratio <= 0.7 && p.n_cys >= 4, props)
    pass_rate = n_pass / length(props)

    push!(summary_rows, (
        Method = method_name,
        N = length(seqs),
        Charge_mean = mean(charges),
        Charge_std = std(charges),
        HydroRatio_mean = mean(hydro_ratios),
        HydroRatio_std = std(hydro_ratios),
        MeanHydro_mean = mean(mean_hydros),
        MeanHydro_std = std(mean_hydros),
        Cys_mean = mean(n_cys_vals),
        Cys_std = std(n_cys_vals),
        AMP_pass_rate = pass_rate,
    ))
end

df_summary = DataFrame(summary_rows)
println("\n══════════════════════════════════════════════════════")
println("DEFENSIN BIOPHYSICAL PROPERTIES")
println("══════════════════════════════════════════════════════")
pretty_table(df_summary[:, [:Method, :N, :Charge_mean, :Charge_std, :HydroRatio_mean, :Cys_mean, :AMP_pass_rate]])

CSV.write(joinpath(DATA_DIR, "biophysical_properties.csv"), df_summary)
@info "Saved biophysical_properties.csv"

# ══════════════════════════════════════════════════════════════════════════════
# Figures
# ══════════════════════════════════════════════════════════════════════════════

@info "Generating figures …"

method_colors = Dict(
    "Stored"           => RGB(0.20, 0.20, 0.20),   # dark gray
    "SA generation"    => RGB(0.20, 0.47, 0.69),   # steel blue
    "HMM emit"         => RGB(0.30, 0.69, 0.29),   # green
    "EvoDiff"          => RGB(0.93, 0.68, 0.20),   # amber
    "MSA Transformer"  => RGB(0.58, 0.32, 0.68),   # purple
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
)

method_order = ["Stored", "SA generation", "HMM emit", "EvoDiff", "MSA Transformer"]
available_methods = [m for m in method_order if haskey(methods, m)]

# (A) Charge vs Hydrophobic Ratio scatter
p_scatter = plot(xlabel="Net charge (pH 7)", ylabel="Hydrophobic ratio",
                 title="(A)  Charge vs. hydrophobic content", legend=:topright;
                 panel_defaults...)
for m in available_methods
    df = results[m]
    scatter!(p_scatter, df.charge, df.hydro_ratio;
             label=m, color=get(method_colors, m, :gray),
             alpha=0.6, markersize=4, markerstrokewidth=0.5)
end
# Add AMP "plausible zone"
vline!(p_scatter, [2.0], ls=:dash, color=:gray60, label="", lw=0.8)
hline!(p_scatter, [0.3, 0.7], ls=:dash, color=:gray60, label="", lw=0.8)

# (B) Cysteine count distribution
p_cys = plot(xlabel="Cysteine count", ylabel="Frequency",
             title="(B)  Cysteine residues (disulfide scaffold)", legend=:topright;
             panel_defaults...)
for m in available_methods
    df = results[m]
    histogram!(p_cys, df.n_cys; bins=0:1:12, alpha=0.5, label=m,
               color=get(method_colors, m, :gray), normalize=:probability)
end
vline!(p_cys, [6], ls=:dash, color=:red, label="Expected (6 Cys)", lw=1.5)

# (C) Charge distribution
p_charge = plot(xlabel="Net charge (pH 7)", ylabel="Frequency",
                title="(C)  Net charge distribution", legend=:topright;
                panel_defaults...)
for m in available_methods
    df = results[m]
    histogram!(p_charge, df.charge; bins=-5:1:15, alpha=0.5, label=m,
               color=get(method_colors, m, :gray), normalize=:probability)
end

# (D) AMP pass rate bar chart
pass_rates = [df_summary[df_summary.Method .== m, :AMP_pass_rate][1] for m in available_methods]
colors_bar = [get(method_colors, m, :gray) for m in available_methods]
p_pass = bar(available_methods, pass_rates;
             ylabel="AMP plausibility pass rate",
             title="(D)  Fraction passing biophysical filter",
             color=colors_bar, alpha=0.85, legend=false,
             xrotation=15, ylims=(0, 1.05),
             panel_defaults...)

# Composite
p_composite = plot(p_scatter, p_cys, p_charge, p_pass;
                   layout=(2, 2), size=(900, 700), dpi=300)
savefig(p_composite, joinpath(FIG_DIR, "defensin_biophysics.pdf"))
@info "Saved defensin_biophysics.pdf"

# Also save individual panels
savefig(p_scatter, joinpath(FIG_DIR, "defensin_charge_hydro.pdf"))
savefig(p_cys, joinpath(FIG_DIR, "defensin_cysteine.pdf"))

@info "Done!"
