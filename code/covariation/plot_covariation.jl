#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Plot covariation analysis figures
# ──────────────────────────────────────────────────────────────────────────────

using Statistics, LinearAlgebra, Printf
using CairoMakie

include(joinpath(@__DIR__, "..", "src", "Protein.jl"))

# ── MI functions (shared with run_covariation_analysis.jl) ────────────────────

function compute_mi_matrix(char_mat::Matrix{Char}; pseudocount::Float64=1.0)
    K, L = size(char_mat)
    q = N_AA
    MI = zeros(L, L)
    freq_single = Matrix{Float64}(undef, q, L)
    for j in 1:L
        counts = fill(pseudocount, q)
        for i in 1:K
            idx = get(AA_TO_IDX, char_mat[i, j], 0)
            idx > 0 && (counts[idx] += 1.0)
        end
        freq_single[:, j] = counts ./ sum(counts)
    end
    for a in 1:L
        for b in (a+1):L
            joint = fill(pseudocount / q, q, q)
            for i in 1:K
                ia = get(AA_TO_IDX, char_mat[i, a], 0)
                ib = get(AA_TO_IDX, char_mat[i, b], 0)
                (ia > 0 && ib > 0) && (joint[ia, ib] += 1.0)
            end
            joint ./= sum(joint)
            mi = 0.0
            for ia in 1:q, ib in 1:q
                pab = joint[ia, ib]
                pa = freq_single[ia, a]
                pb = freq_single[ib, b]
                (pab > 0 && pa > 0 && pb > 0) && (mi += pab * log(pab / (pa * pb)))
            end
            MI[a, b] = max(mi, 0.0)
            MI[b, a] = MI[a, b]
        end
    end
    return MI
end

function seqs_to_char_mat(seqs::Vector{String}, L::Int)
    K = length(seqs)
    mat = fill('-', K, L)
    for (i, seq) in enumerate(seqs)
        for (j, c) in enumerate(seq)
            j <= L && (mat[i, j] = c)
        end
    end
    return mat
end

function upper_triangle(M::Matrix{Float64})
    L = size(M, 1)
    vals = Float64[]
    for i in 1:L, j in (i+1):L
        push!(vals, M[i, j])
    end
    return vals
end

function spearman_corr(x::Vector{Float64}, y::Vector{Float64})
    rx = sortperm(sortperm(x)) .|> Float64
    ry = sortperm(sortperm(y)) .|> Float64
    return cor(rx, ry)
end

const THREE_TO_ONE = Dict(
    "ALA" => 'A', "CYS" => 'C', "ASP" => 'D', "GLU" => 'E', "PHE" => 'F',
    "GLY" => 'G', "HIS" => 'H', "ILE" => 'I', "LYS" => 'K', "LEU" => 'L',
    "MET" => 'M', "ASN" => 'N', "PRO" => 'P', "GLN" => 'Q', "ARG" => 'R',
    "SER" => 'S', "THR" => 'T', "VAL" => 'V', "TRP" => 'W', "TYR" => 'Y'
)

function extract_ca_from_pdb(pdb_path::String)
    resnums = Int[]; aas = Char[]; coords = Float64[]
    for line in eachline(pdb_path)
        startswith(line, "ATOM") || continue
        length(line) < 54 && continue
        strip(line[13:16]) == "CA" || continue
        aa = get(THREE_TO_ONE, strip(line[18:20]), 'X')
        push!(resnums, parse(Int, strip(line[23:26])))
        push!(aas, aa)
        append!(coords, [parse(Float64, strip(line[31:38])),
                         parse(Float64, strip(line[39:46])),
                         parse(Float64, strip(line[47:54]))])
    end
    N = length(resnums)
    return resnums, aas, Matrix(reshape(coords, 3, N)')
end

function contact_pairs(coords::Matrix{Float64}; threshold::Float64=8.0)
    N = size(coords, 1)
    pairs = Tuple{Int,Int}[]
    for i in 1:N, j in (i+1):N
        abs(i - j) < 5 && continue
        sqrt(sum((coords[i, :] .- coords[j, :]).^2)) < threshold && push!(pairs, (i, j))
    end
    return pairs
end

function align_pdb_to_msa(pdb_aas::Vector{Char}, char_mat::Matrix{Char})
    K, L = size(char_mat)
    pdb_len = length(pdb_aas)
    consensus = [AA_ALPHABET[argmax([count(==(AA_ALPHABET[a]), char_mat[:, j]) for a in 1:N_AA])] for j in 1:L]
    best_score = -1; best_offset = 0
    for offset in 0:max(0, pdb_len - L)
        score = sum(pdb_aas[offset + j] == consensus[j] for j in 1:min(L, pdb_len - offset))
        if score > best_score; best_score = score; best_offset = offset; end
    end
    if pdb_len <= L
        for offset in 0:(L - pdb_len)
            score = sum(pdb_aas[j] == consensus[offset + j] for j in 1:pdb_len)
            if score > best_score; best_score = score; best_offset = -offset; end
        end
    end
    best_score / min(L, pdb_len) < 0.3 && return nothing
    mapping = zeros(Int, L)
    if best_offset >= 0
        for j in 1:L; idx = best_offset + j; 1 <= idx <= pdb_len && (mapping[j] = idx); end
    else
        off = -best_offset
        for j in 1:L; idx = j - off; 1 <= idx <= pdb_len && (mapping[j] = idx); end
    end
    return mapping
end

# ── Configuration ─────────────────────────────────────────────────────────────

base_dir = joinpath(@__DIR__, "..")
data_dir = joinpath(base_dir, "pfam-families", "data")
baselines_dir = joinpath(base_dir, "baselines", "data")
structure_dir = joinpath(base_dir, "structure-validation", "data")
fig_dir = joinpath(@__DIR__, "..", "..", "arxiv", "figs")
mkpath(fig_dir)

# representative families for scatter plots
scatter_families = [
    ("PF00014", "Kunitz"),
    ("PF00018", "SH3"),
    ("PF00711", "Defensin_beta"),
]

# ══════════════════════════════════════════════════════════════════════════════
# Figure: MI scatter plots for representative families
# ══════════════════════════════════════════════════════════════════════════════

fig = Figure(size=(1200, 400), fontsize=12)

for (idx, (pfam_id, family_name)) in enumerate(scatter_families)
    # load data
    sto_file = joinpath(data_dir, pfam_id, "$(pfam_id)_seed.sto")
    raw_seqs = parse_stockholm(sto_file)
    char_mat_stored, _ = clean_alignment(raw_seqs)
    K_stored, L = size(char_mat_stored)

    sa_gen_file = joinpath(data_dir, pfam_id, "sa_generation_sequences.fasta")
    sa_gen_raw = parse_fasta(sa_gen_file)
    sa_gen_seqs = [s[2] for s in sa_gen_raw]
    char_mat_sa_gen = seqs_to_char_mat(sa_gen_seqs, L)

    MI_stored = compute_mi_matrix(char_mat_stored)
    MI_sa_gen = compute_mi_matrix(char_mat_sa_gen)

    v_stored = upper_triangle(MI_stored)
    v_sa_gen = upper_triangle(MI_sa_gen)

    r = cor(v_stored, v_sa_gen)
    ρ = spearman_corr(v_stored, v_sa_gen)

    ax = Axis(fig[1, idx],
        xlabel="MI (stored sequences)",
        ylabel="MI (SA generation)",
        title="$family_name (r=$(round(r, digits=2)), ρ=$(round(ρ, digits=2)))",
        aspect=1)

    # check for Level 2 contacts
    pdb_refs_local = Dict(
        "PF00711" => "1E4S_A.pdb",
        "PF00014" => "1BPI_A.pdb",
        "PF00018" => "1SHG_A.pdb",
    )

    has_contacts = false
    contact_aln_pairs = Set{Tuple{Int,Int}}()

    if haskey(pdb_refs_local, pfam_id)
        pdb_file = joinpath(structure_dir, pfam_id, "reference", pdb_refs_local[pfam_id])
        if isfile(pdb_file)
            resnums, pdb_aas, coords = extract_ca_from_pdb(pdb_file)
            mapping = align_pdb_to_msa(pdb_aas, char_mat_stored)
            if mapping !== nothing
                contacts = contact_pairs(coords; threshold=8.0)
                for (pi, pj) in contacts
                    ai = findfirst(==(pi), mapping)
                    aj = findfirst(==(pj), mapping)
                    if ai !== nothing && aj !== nothing
                        a, b = minmax(ai, aj)
                        push!(contact_aln_pairs, (a, b))
                    end
                end
                has_contacts = length(contact_aln_pairs) > 0
            end
        end
    end

    if has_contacts
        # separate contact and non-contact pairs
        contact_stored = Float64[]
        contact_gen = Float64[]
        noncontact_stored = Float64[]
        noncontact_gen = Float64[]

        k = 0
        for i in 1:L
            for j in (i+1):L
                k += 1
                if (i, j) in contact_aln_pairs
                    push!(contact_stored, v_stored[k])
                    push!(contact_gen, v_sa_gen[k])
                else
                    push!(noncontact_stored, v_stored[k])
                    push!(noncontact_gen, v_sa_gen[k])
                end
            end
        end

        scatter!(ax, noncontact_stored, noncontact_gen,
            color=(:gray, 0.2), markersize=3, label="Non-contact")
        scatter!(ax, contact_stored, contact_gen,
            color=(:red, 0.6), markersize=5, label="Contact (<8Å)")
    else
        scatter!(ax, v_stored, v_sa_gen,
            color=(RGBf(0.20, 0.47, 0.69), 0.15), markersize=3)
    end

    # identity line
    max_val = max(maximum(v_stored), maximum(v_sa_gen))
    lines!(ax, [0, max_val], [0, max_val], color=:black, linestyle=:dash, linewidth=1)

    if has_contacts
        axislegend(ax, position=:rb, framevisible=false, labelsize=9)
    end
end

# add panel labels
for (i, label) in enumerate(["A", "B", "C"])
    Label(fig[1, i, TopLeft()], label, fontsize=16, font=:bold, padding=(0, 5, 5, 0))
end

save(joinpath(fig_dir, "mi_covariation_scatter.pdf"), fig)
@info "Saved figure to $(joinpath(fig_dir, "mi_covariation_scatter.pdf"))"

# ══════════════════════════════════════════════════════════════════════════════
# Figure: Bar chart comparing methods across families
# ══════════════════════════════════════════════════════════════════════════════

fig2 = Figure(size=(1100, 450), fontsize=13)
ax2 = Axis(fig2[1, 1],
    ylabel="Pearson r (MI correlation with stored)",
    xticks=(1:8, ["RRM", "SH3", "WW", "Kunitz", "zf-C2H2", "PDZ", "Pkinase", "Def_β"]),
    xticklabelrotation=π/6,
    title="Pairwise MI preservation across families and methods",
    limits=(0.5, 8.5, -0.2, 1.05))

all_families = [
    ("PF00076", "RRM"),
    ("PF00018", "SH3"),
    ("PF00397", "WW"),
    ("PF00014", "Kunitz"),
    ("PF00096", "zf-C2H2"),
    ("PF00595", "PDZ"),
    ("PF00069", "Pkinase"),
    ("PF00711", "Defensin_beta"),
]

# 5 methods: SA gen, HMM, EvoDiff, MSAT, Potts (last only where available)
methods = ["SA gen", "HMM", "EvoDiff", "MSAT", "Potts"]
method_colors = [
    RGBf(0.20, 0.47, 0.69),   # SA gen - steel blue (canonical)
    RGBf(0.30, 0.69, 0.29),   # HMM - green (canonical)
    RGBf(0.93, 0.68, 0.20),   # EvoDiff - amber (canonical)
    RGBf(0.58, 0.32, 0.68),   # MSAT - purple (canonical)
    RGBf(0.84, 0.15, 0.16),   # Potts - crimson (canonical)
]
bar_width = 0.15
offsets = [-2.0, -1.0, 0.0, 1.0, 2.0] .* bar_width

for (fam_idx, (pfam_id, family_name)) in enumerate(all_families)
    sto_file = joinpath(data_dir, pfam_id, "$(pfam_id)_seed.sto")
    raw_seqs = parse_stockholm(sto_file)
    char_mat_stored, _ = clean_alignment(raw_seqs)
    K_stored, L = size(char_mat_stored)
    MI_stored = compute_mi_matrix(char_mat_stored)
    v_stored = upper_triangle(MI_stored)

    files = [
        joinpath(data_dir, pfam_id, "sa_generation_sequences.fasta"),
        joinpath(baselines_dir, pfam_id, "hmm_sequences.fasta"),
        joinpath(baselines_dir, pfam_id, "evodiff_sequences.fasta"),
        joinpath(baselines_dir, pfam_id, "msat_sequences.fasta"),
    ]

    for (m_idx, fpath) in enumerate(files)
        if isfile(fpath)
            raw = parse_fasta(fpath)
            seqs = [s[2] for s in raw]
            cmat = seqs_to_char_mat(seqs, L)
            mi = compute_mi_matrix(cmat)
            v = upper_triangle(mi)
            r_val = cor(v_stored, v)
            barplot!(ax2, [fam_idx + offsets[m_idx]], [r_val],
                color=method_colors[m_idx], width=bar_width,
                label=fam_idx == 1 ? methods[m_idx] : nothing)
        end
    end

    # Potts as 5th bar (only where available)
    potts_file = joinpath(baselines_dir, pfam_id, "potts_sequences.fasta")
    if isfile(potts_file)
        raw = parse_fasta(potts_file)
        seqs = [s[2] for s in raw]
        cmat = seqs_to_char_mat(seqs, L)
        mi = compute_mi_matrix(cmat)
        v = upper_triangle(mi)
        r_val = cor(v_stored, v)
        barplot!(ax2, [fam_idx + offsets[5]], [r_val],
            color=method_colors[5], width=bar_width,
            label=pfam_id == "PF00018" ? methods[5] : nothing)
    end
end

hlines!(ax2, [0], color=:gray, linestyle=:dash, linewidth=0.5)
Legend(fig2[1, 2], ax2, framevisible=false, labelsize=12)

save(joinpath(fig_dir, "mi_covariation_barplot.pdf"), fig2)
@info "Saved figure to $(joinpath(fig_dir, "mi_covariation_barplot.pdf"))"

@info "All figures generated."
