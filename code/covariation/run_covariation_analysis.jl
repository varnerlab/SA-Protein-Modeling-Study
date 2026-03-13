#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Covariation analysis: MI matrices for stored vs generated sequences
#
# Computes pairwise mutual information (MI) for stored and SA-generated
# sequences, correlates them, and optionally maps to structural contacts.
# Also evaluates Potts and HMM baselines where available.
# ──────────────────────────────────────────────────────────────────────────────

using Statistics, LinearAlgebra, Printf, DelimitedFiles

# ── load shared protein pipeline ──────────────────────────────────────────────
include(joinpath(@__DIR__, "..", "src", "Protein.jl"))

# ══════════════════════════════════════════════════════════════════════════════
# MI computation
# ══════════════════════════════════════════════════════════════════════════════

"""
    compute_mi_matrix(char_mat; pseudocount=1.0) -> Matrix{Float64}

Compute the L×L pairwise mutual information matrix from a K×L character matrix.
Uses additive pseudocount correction to avoid log(0).
MI(i,j) = sum_{a,b} p(a,b) * log(p(a,b) / (p_i(a) * p_j(b)))
"""
function compute_mi_matrix(char_mat::Matrix{Char}; pseudocount::Float64=1.0)
    K, L = size(char_mat)
    q = N_AA  # 20
    MI = zeros(L, L)

    # precompute per-position frequency vectors (with pseudocount)
    freq_single = Matrix{Float64}(undef, q, L)
    for j in 1:L
        counts = fill(pseudocount, q)
        for i in 1:K
            idx = get(AA_TO_IDX, char_mat[i, j], 0)
            idx > 0 && (counts[idx] += 1.0)
        end
        freq_single[:, j] = counts ./ sum(counts)
    end

    # compute MI for each pair (upper triangle)
    for a in 1:L
        for b in (a+1):L
            # joint frequency table q×q
            joint = fill(pseudocount / q, q, q)  # spread pseudocount across joint cells
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
                if pab > 0 && pa > 0 && pb > 0
                    mi += pab * log(pab / (pa * pb))
                end
            end
            MI[a, b] = max(mi, 0.0)  # MI is non-negative; clamp numerical noise
            MI[b, a] = MI[a, b]
        end
    end

    return MI
end

"""
    seqs_to_char_mat(seqs, L) -> Matrix{Char}

Convert a vector of sequence strings to a K×L character matrix.
Sequences shorter than L are padded with '-'; longer are truncated.
"""
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

"""
    upper_triangle(M) -> Vector{Float64}

Extract upper-triangle elements (i < j) of a symmetric matrix.
"""
function upper_triangle(M::Matrix{Float64})
    L = size(M, 1)
    vals = Float64[]
    for i in 1:L
        for j in (i+1):L
            push!(vals, M[i, j])
        end
    end
    return vals
end

# ══════════════════════════════════════════════════════════════════════════════
# Correlation helpers
# ══════════════════════════════════════════════════════════════════════════════

"""
    spearman_corr(x, y) -> Float64

Spearman rank correlation.
"""
function spearman_corr(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    rx = sortperm(sortperm(x)) .|> Float64
    ry = sortperm(sortperm(y)) .|> Float64
    return cor(rx, ry)
end

"""
    top_k_overlap(v1, v2, k) -> Float64

Fraction of top-k pairs in v1 that are also in top-k of v2.
v1 and v2 are upper-triangle MI vectors (same ordering).
"""
function top_k_overlap(v1::Vector{Float64}, v2::Vector{Float64}, k::Int)
    k = min(k, length(v1))
    top1 = Set(partialsortperm(v1, 1:k, rev=true))
    top2 = Set(partialsortperm(v2, 1:k, rev=true))
    return length(intersect(top1, top2)) / k
end

# ══════════════════════════════════════════════════════════════════════════════
# PDB contact extraction (Level 2)
# ══════════════════════════════════════════════════════════════════════════════

const THREE_TO_ONE = Dict(
    "ALA" => 'A', "CYS" => 'C', "ASP" => 'D', "GLU" => 'E', "PHE" => 'F',
    "GLY" => 'G', "HIS" => 'H', "ILE" => 'I', "LYS" => 'K', "LEU" => 'L',
    "MET" => 'M', "ASN" => 'N', "PRO" => 'P', "GLN" => 'Q', "ARG" => 'R',
    "SER" => 'S', "THR" => 'T', "VAL" => 'V', "TRP" => 'W', "TYR" => 'Y'
)

"""
    extract_ca_from_pdb(pdb_path) -> (Vector{Int}, Vector{Char}, Matrix{Float64})

Extract Cα atoms from a PDB file.
Returns (residue_numbers, aa_sequence, coords) where coords is N×3.
"""
function extract_ca_from_pdb(pdb_path::String)
    resnums = Int[]
    aas = Char[]
    coords = Float64[]

    for line in eachline(pdb_path)
        startswith(line, "ATOM") || continue
        length(line) < 54 && continue
        atom_name = strip(line[13:16])
        atom_name == "CA" || continue

        resname = strip(line[18:20])
        aa = get(THREE_TO_ONE, resname, 'X')
        resnum = parse(Int, strip(line[23:26]))
        x = parse(Float64, strip(line[31:38]))
        y = parse(Float64, strip(line[39:46]))
        z = parse(Float64, strip(line[47:54]))

        push!(resnums, resnum)
        push!(aas, aa)
        append!(coords, [x, y, z])
    end

    N = length(resnums)
    coords_mat = reshape(coords, 3, N)'  # N×3
    return resnums, aas, Matrix(coords_mat)
end

"""
    contact_pairs(coords; threshold=8.0) -> Vector{Tuple{Int,Int}}

Find pairs of residues with Cα-Cα distance < threshold (Å).
Returns 1-based indices into the coords array.
"""
function contact_pairs(coords::Matrix{Float64}; threshold::Float64=8.0)
    N = size(coords, 1)
    pairs = Tuple{Int,Int}[]
    for i in 1:N
        for j in (i+1):N
            abs(i - j) < 5 && continue  # skip short-range (< 5 residues apart)
            d = sqrt(sum((coords[i, :] .- coords[j, :]).^2))
            d < threshold && push!(pairs, (i, j))
        end
    end
    return pairs
end

"""
    align_pdb_to_msa(pdb_aas, char_mat) -> Union{Vector{Int}, Nothing}

Find the best offset alignment of PDB sequence to MSA consensus.
Returns a mapping: alignment_position -> PDB_index (1-based), or nothing if no good match.
Only works when PDB length is close to MSA length.
"""
function align_pdb_to_msa(pdb_aas::Vector{Char}, char_mat::Matrix{Char})
    K, L = size(char_mat)
    pdb_len = length(pdb_aas)

    # compute MSA consensus
    consensus = Char[]
    for j in 1:L
        counts = zeros(Int, N_AA)
        for i in 1:K
            idx = get(AA_TO_IDX, char_mat[i, j], 0)
            idx > 0 && (counts[idx] += 1)
        end
        push!(consensus, AA_ALPHABET[argmax(counts)])
    end

    # try all offsets of MSA within PDB sequence
    best_score = -1
    best_offset = 0
    for offset in 0:(pdb_len - L)
        score = sum(pdb_aas[offset + j] == consensus[j] for j in 1:L)
        if score > best_score
            best_score = score
            best_offset = offset
        end
    end

    # also try PDB within MSA (if PDB shorter)
    best_msa_offset = 0
    if pdb_len <= L
        for offset in 0:(L - pdb_len)
            score = sum(pdb_aas[j] == consensus[offset + j] for j in 1:pdb_len)
            if score > best_score
                best_score = score
                best_offset = -offset  # negative means PDB starts before MSA
                best_msa_offset = offset
            end
        end
    end

    match_frac = best_score / min(L, pdb_len)
    if match_frac < 0.3
        @warn "  Poor PDB-MSA alignment: $(round(match_frac*100, digits=1))% match"
        return nothing
    end

    @info "  PDB-MSA alignment: offset=$best_offset, match=$(round(match_frac*100, digits=1))%"

    # build mapping: alignment position -> PDB index
    mapping = zeros(Int, L)
    if best_offset >= 0
        for j in 1:L
            pdb_idx = best_offset + j
            if 1 <= pdb_idx <= pdb_len
                mapping[j] = pdb_idx
            end
        end
    else
        msa_off = -best_offset
        for j in 1:L
            pdb_idx = j - msa_off
            if 1 <= pdb_idx <= pdb_len
                mapping[j] = pdb_idx
            end
        end
    end

    return mapping
end

# ══════════════════════════════════════════════════════════════════════════════
# Main analysis
# ══════════════════════════════════════════════════════════════════════════════

# family configuration
families = [
    ("PF00076", "RRM"),
    ("PF00018", "SH3"),
    ("PF00397", "WW"),
    ("PF00014", "Kunitz"),
    ("PF00096", "zf-C2H2"),
    ("PF00595", "PDZ"),
    ("PF00069", "Pkinase"),
    ("PF00711", "Defensin_beta"),
]

# PDB reference files for Level 2
pdb_refs = Dict(
    "PF00711" => "1E4S_A.pdb",
    "PF00014" => "1BPI_A.pdb",
    "PF00018" => "1SHG_A.pdb",
    "PF00397" => "1PIN_A.pdb",
    "PF00595" => "1BE9_A.pdb",
    "PF00076" => "1FXL_A.pdb",
    "PF00096" => "1ZAA_C.pdb",
    "PF00069" => "1ATP_E.pdb",
)

# families suitable for Level 2 contact analysis (PDB size close to L)
level2_families = Set(["PF00711", "PF00014", "PF00018"])

base_dir = joinpath(@__DIR__, "..")
data_dir = joinpath(base_dir, "pfam-families", "data")
baselines_dir = joinpath(base_dir, "baselines", "data")
structure_dir = joinpath(base_dir, "structure-validation", "data")
output_dir = joinpath(@__DIR__, "results")
mkpath(output_dir)

# results storage
results = []

for (pfam_id, family_name) in families
    @info "═══ Processing $family_name ($pfam_id) ═══"

    # ── load stored alignment ─────────────────────────────────────────────
    sto_file = joinpath(data_dir, pfam_id, "$(pfam_id)_seed.sto")
    raw_seqs = parse_stockholm(sto_file)
    char_mat_stored, _ = clean_alignment(raw_seqs)
    K_stored, L = size(char_mat_stored)

    # ── load SA generation sequences ──────────────────────────────────────
    sa_gen_file = joinpath(data_dir, pfam_id, "sa_generation_sequences.fasta")
    sa_gen_raw = parse_fasta(sa_gen_file)
    sa_gen_seqs = [s[2] for s in sa_gen_raw]
    char_mat_sa_gen = seqs_to_char_mat(sa_gen_seqs, L)

    # ── compute MI matrices ───────────────────────────────────────────────
    @info "  Computing MI matrices (L=$L, K_stored=$K_stored, K_gen=$(length(sa_gen_seqs)))..."
    MI_stored = compute_mi_matrix(char_mat_stored)
    MI_sa_gen = compute_mi_matrix(char_mat_sa_gen)

    # ── correlations ──────────────────────────────────────────────────────
    v_stored = upper_triangle(MI_stored)
    v_sa_gen = upper_triangle(MI_sa_gen)

    r_pearson = cor(v_stored, v_sa_gen)
    r_spearman = spearman_corr(v_stored, v_sa_gen)

    # top-N overlap
    n_pairs = length(v_stored)
    top_50 = top_k_overlap(v_stored, v_sa_gen, 50)
    top_100 = top_k_overlap(v_stored, v_sa_gen, min(100, n_pairs))

    @info "  SA generation: Pearson r=$(round(r_pearson, digits=3)), Spearman ρ=$(round(r_spearman, digits=3))"
    @info "  Top-50 overlap: $(round(top_50*100, digits=1))%, Top-100 overlap: $(round(top_100*100, digits=1))%"

    row = Dict(
        "family" => family_name,
        "pfam_id" => pfam_id,
        "L" => L,
        "K_stored" => K_stored,
        "n_pairs" => n_pairs,
        "sa_gen_pearson" => r_pearson,
        "sa_gen_spearman" => r_spearman,
        "sa_gen_top50" => top_50,
        "sa_gen_top100" => top_100,
    )

    # ── SA retrieval ──────────────────────────────────────────────────────
    sa_ret_file = joinpath(data_dir, pfam_id, "sa_retrieval_sequences.fasta")
    if isfile(sa_ret_file)
        sa_ret_raw = parse_fasta(sa_ret_file)
        sa_ret_seqs = [s[2] for s in sa_ret_raw]
        char_mat_sa_ret = seqs_to_char_mat(sa_ret_seqs, L)
        MI_sa_ret = compute_mi_matrix(char_mat_sa_ret)
        v_sa_ret = upper_triangle(MI_sa_ret)
        row["sa_ret_pearson"] = cor(v_stored, v_sa_ret)
        row["sa_ret_spearman"] = spearman_corr(v_stored, v_sa_ret)
        @info "  SA retrieval: Pearson r=$(round(row["sa_ret_pearson"], digits=3)), Spearman ρ=$(round(row["sa_ret_spearman"], digits=3))"
    end

    # ── Potts baseline (if available) ─────────────────────────────────────
    potts_file = joinpath(baselines_dir, pfam_id, "potts_sequences.fasta")
    if isfile(potts_file)
        potts_raw = parse_fasta(potts_file)
        potts_seqs = [s[2] for s in potts_raw]
        char_mat_potts = seqs_to_char_mat(potts_seqs, L)
        MI_potts = compute_mi_matrix(char_mat_potts)
        v_potts = upper_triangle(MI_potts)
        row["potts_pearson"] = cor(v_stored, v_potts)
        row["potts_spearman"] = spearman_corr(v_stored, v_potts)
        row["potts_top50"] = top_k_overlap(v_stored, v_potts, 50)
        @info "  Potts: Pearson r=$(round(row["potts_pearson"], digits=3)), Spearman ρ=$(round(row["potts_spearman"], digits=3))"
    end

    # ── HMM baseline ──────────────────────────────────────────────────────
    hmm_file = joinpath(baselines_dir, pfam_id, "hmm_sequences.fasta")
    if isfile(hmm_file)
        hmm_raw = parse_fasta(hmm_file)
        hmm_seqs = [s[2] for s in hmm_raw]
        char_mat_hmm = seqs_to_char_mat(hmm_seqs, L)
        MI_hmm = compute_mi_matrix(char_mat_hmm)
        v_hmm = upper_triangle(MI_hmm)
        row["hmm_pearson"] = cor(v_stored, v_hmm)
        row["hmm_spearman"] = spearman_corr(v_stored, v_hmm)
        row["hmm_top50"] = top_k_overlap(v_stored, v_hmm, 50)
        @info "  HMM: Pearson r=$(round(row["hmm_pearson"], digits=3)), Spearman ρ=$(round(row["hmm_spearman"], digits=3))"
    end

    # ── EvoDiff baseline ──────────────────────────────────────────────────
    evodiff_file = joinpath(baselines_dir, pfam_id, "evodiff_sequences.fasta")
    if isfile(evodiff_file)
        evodiff_raw = parse_fasta(evodiff_file)
        evodiff_seqs = [s[2] for s in evodiff_raw]
        char_mat_evodiff = seqs_to_char_mat(evodiff_seqs, L)
        MI_evodiff = compute_mi_matrix(char_mat_evodiff)
        v_evodiff = upper_triangle(MI_evodiff)
        row["evodiff_pearson"] = cor(v_stored, v_evodiff)
        row["evodiff_spearman"] = spearman_corr(v_stored, v_evodiff)
        row["evodiff_top50"] = top_k_overlap(v_stored, v_evodiff, 50)
        @info "  EvoDiff: Pearson r=$(round(row["evodiff_pearson"], digits=3)), Spearman ρ=$(round(row["evodiff_spearman"], digits=3))"
    end

    # ── MSA Transformer baseline ──────────────────────────────────────────
    msat_file = joinpath(baselines_dir, pfam_id, "msat_sequences.fasta")
    if isfile(msat_file)
        msat_raw = parse_fasta(msat_file)
        msat_seqs = [s[2] for s in msat_raw]
        char_mat_msat = seqs_to_char_mat(msat_seqs, L)
        MI_msat = compute_mi_matrix(char_mat_msat)
        v_msat = upper_triangle(MI_msat)
        row["msat_pearson"] = cor(v_stored, v_msat)
        row["msat_spearman"] = spearman_corr(v_stored, v_msat)
        row["msat_top50"] = top_k_overlap(v_stored, v_msat, 50)
        @info "  MSAT: Pearson r=$(round(row["msat_pearson"], digits=3)), Spearman ρ=$(round(row["msat_spearman"], digits=3))"
    end

    # ── Level 2: structural contact analysis ──────────────────────────────
    if pfam_id in level2_families
        pdb_file = joinpath(structure_dir, pfam_id, "reference", pdb_refs[pfam_id])
        if isfile(pdb_file)
            @info "  Level 2: extracting contacts from $(pdb_refs[pfam_id])..."
            resnums, pdb_aas, coords = extract_ca_from_pdb(pdb_file)
            mapping = align_pdb_to_msa(pdb_aas, char_mat_stored)

            if mapping !== nothing
                contacts = contact_pairs(coords; threshold=8.0)
                @info "  Found $(length(contacts)) structural contacts (>5 residues apart, <8Å)"

                # map PDB contacts to alignment positions
                contact_mi_stored = Float64[]
                contact_mi_gen = Float64[]
                noncontact_mi_stored = Float64[]
                noncontact_mi_gen = Float64[]

                # build set of alignment-position contact pairs
                contact_aln_pairs = Set{Tuple{Int,Int}}()
                for (pi, pj) in contacts
                    # find alignment positions mapping to these PDB indices
                    ai = findfirst(==(pi), mapping)
                    aj = findfirst(==(pj), mapping)
                    if ai !== nothing && aj !== nothing
                        a, b = minmax(ai, aj)
                        push!(contact_aln_pairs, (a, b))
                    end
                end
                @info "  Mapped $(length(contact_aln_pairs)) contacts to alignment positions"

                for i in 1:L
                    for j in (i+1):L
                        if (i, j) in contact_aln_pairs
                            push!(contact_mi_stored, MI_stored[i, j])
                            push!(contact_mi_gen, MI_sa_gen[i, j])
                        else
                            push!(noncontact_mi_stored, MI_stored[i, j])
                            push!(noncontact_mi_gen, MI_sa_gen[i, j])
                        end
                    end
                end

                if length(contact_mi_stored) > 0
                    row["n_contacts"] = length(contact_mi_stored)
                    row["contact_mi_stored_mean"] = mean(contact_mi_stored)
                    row["contact_mi_gen_mean"] = mean(contact_mi_gen)
                    row["noncontact_mi_stored_mean"] = mean(noncontact_mi_stored)
                    row["noncontact_mi_gen_mean"] = mean(noncontact_mi_gen)
                    row["contact_enrichment_stored"] = mean(contact_mi_stored) / max(mean(noncontact_mi_stored), 1e-10)
                    row["contact_enrichment_gen"] = mean(contact_mi_gen) / max(mean(noncontact_mi_gen), 1e-10)
                    row["contact_mi_corr"] = cor(contact_mi_stored, contact_mi_gen)

                    @info "  Contact MI (stored): $(round(row["contact_mi_stored_mean"], digits=4)) vs non-contact: $(round(row["noncontact_mi_stored_mean"], digits=4))"
                    @info "  Contact MI (SA gen): $(round(row["contact_mi_gen_mean"], digits=4)) vs non-contact: $(round(row["noncontact_mi_gen_mean"], digits=4))"
                    @info "  Contact MI correlation (stored vs gen): $(round(row["contact_mi_corr"], digits=3))"
                end
            end
        end
    end

    push!(results, row)
    println()
end

# ══════════════════════════════════════════════════════════════════════════════
# Write CSV output
# ══════════════════════════════════════════════════════════════════════════════

csv_file = joinpath(output_dir, "mi_correlation_results.csv")
@info "Writing results to $csv_file"

# collect all keys across all rows
all_keys = String[]
for row in results
    for k in keys(row)
        k in all_keys || push!(all_keys, k)
    end
end
sort!(all_keys)

open(csv_file, "w") do f
    println(f, join(all_keys, ","))
    for row in results
        vals = [get(row, k, "") for k in all_keys]
        println(f, join(string.(vals), ","))
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# Print summary table
# ══════════════════════════════════════════════════════════════════════════════

println("\n" * "="^80)
println("COVARIATION ANALYSIS SUMMARY")
println("="^80)

# Level 1: MI correlation across all families
println("\nLevel 1: Pairwise MI correlation (stored vs generated)")
println("-"^90)
@printf("%-15s %5s %6s  %8s %8s  %8s %8s  %6s\n",
    "Family", "L", "K", "SA_r", "SA_ρ", "Potts_r", "HMM_r", "Top50")
println("-"^90)

sa_pearson_vals = Float64[]
sa_spearman_vals = Float64[]

for row in results
    potts_r = get(row, "potts_pearson", "")
    hmm_r = get(row, "hmm_pearson", "")
    @printf("%-15s %5d %6d  %8.3f %8.3f  %8s %8s  %5.1f%%\n",
        row["family"], row["L"], row["K_stored"],
        row["sa_gen_pearson"], row["sa_gen_spearman"],
        potts_r == "" ? "---" : @sprintf("%.3f", potts_r),
        hmm_r == "" ? "---" : @sprintf("%.3f", hmm_r),
        row["sa_gen_top50"] * 100)
    push!(sa_pearson_vals, row["sa_gen_pearson"])
    push!(sa_spearman_vals, row["sa_gen_spearman"])
end

println("-"^90)
@printf("%-15s %5s %6s  %8.3f %8.3f\n", "Mean", "", "",
    mean(sa_pearson_vals), mean(sa_spearman_vals))
@printf("%-15s %5s %6s  %8.3f %8.3f\n", "Range", "", "",
    minimum(sa_pearson_vals), minimum(sa_spearman_vals))
@printf("%-15s %5s %6s  %8.3f %8.3f\n", "  to", "", "",
    maximum(sa_pearson_vals), maximum(sa_spearman_vals))

# Level 2: Contact analysis
println("\n\nLevel 2: MI at structural contacts")
println("-"^70)
has_contacts = [row for row in results if haskey(row, "n_contacts")]
if !isempty(has_contacts)
    @printf("%-15s %6s  %12s %12s  %12s\n",
        "Family", "N_cont", "MI_cont(sto)", "MI_cont(gen)", "Contact_r")
    println("-"^70)
    for row in has_contacts
        @printf("%-15s %6d  %12.4f %12.4f  %12.3f\n",
            row["family"], row["n_contacts"],
            row["contact_mi_stored_mean"], row["contact_mi_gen_mean"],
            row["contact_mi_corr"])
    end
end

# Full baseline comparison table
println("\n\nFull method comparison: Pearson r of MI (stored vs method)")
println("-"^85)
@printf("%-15s  %8s %8s %8s %8s %8s %8s\n",
    "Family", "SA_gen", "SA_ret", "Potts", "HMM", "EvoDiff", "MSAT")
println("-"^85)
for row in results
    @printf("%-15s  %8s %8s %8s %8s %8s %8s\n",
        row["family"],
        @sprintf("%.3f", row["sa_gen_pearson"]),
        haskey(row, "sa_ret_pearson") ? @sprintf("%.3f", row["sa_ret_pearson"]) : "---",
        haskey(row, "potts_pearson") ? @sprintf("%.3f", row["potts_pearson"]) : "---",
        haskey(row, "hmm_pearson") ? @sprintf("%.3f", row["hmm_pearson"]) : "---",
        haskey(row, "evodiff_pearson") ? @sprintf("%.3f", row["evodiff_pearson"]) : "---",
        haskey(row, "msat_pearson") ? @sprintf("%.3f", row["msat_pearson"]) : "---")
end

@info "\nDone. Results saved to $csv_file"
