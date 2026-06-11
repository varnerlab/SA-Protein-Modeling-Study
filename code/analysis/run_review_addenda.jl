#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Review-addenda analysis (no new generation; reuses existing FASTA outputs).
#
# Produces, per family, using the paper's exact sequence_identity definition:
#   (1) Empirical WITHIN-FAMILY nearest-neighbor identity band: for each stored
#       sequence, max identity to any OTHER stored sequence (leave-one-out NN).
#       This is the apples-to-apples comparator to the generated "SeqID"
#       (nearest-neighbor identity), unlike the MEAN pairwise "Stored ID".
#   (2) Sequence-space diversity among the 150 SA-generated sequences:
#       mean pairwise identity, mean pairwise Hamming-distance fraction, and
#       number of distinct decoded sequences (replaces high-dim cosine claim).
#   (3) Per-SA-generated-sequence identity to the family consensus (most
#       frequent residue per column among stored), written to CSV for the
#       consensus-vs-pLDDT regression.
# ──────────────────────────────────────────────────────────────────────────────

using Statistics
using Printf

const PFAM_DATA_DIR = joinpath(@__DIR__, "..", "pfam-families", "data")
const OUT_DIR = joinpath(@__DIR__, "data")
mkpath(OUT_DIR)

const AA_ALPHABET = collect("ACDEFGHIKLMNPQRSTVWY")
const AA_TO_IDX = Dict(aa => i for (i, aa) in enumerate(AA_ALPHABET))

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

# ---- exact copies of the paper's parse_fasta and sequence_identity ----
function parse_fasta(filepath::String)
    sequences = Tuple{String,String}[]
    current_name = ""
    current_seq = IOBuffer()
    for line in eachline(filepath)
        if startswith(line, ">")
            if !isempty(current_name)
                push!(sequences, (current_name, String(take!(current_seq))))
            end
            current_name = strip(line[2:end])
            current_seq = IOBuffer()
        else
            write(current_seq, uppercase(strip(line)))
        end
    end
    if !isempty(current_name)
        push!(sequences, (current_name, String(take!(current_seq))))
    end
    return sequences
end

function sequence_identity(seq1::String, seq2::String)
    L = min(length(seq1), length(seq2))
    matches = 0; compared = 0
    for i in 1:L
        (seq1[i] == '-' || seq2[i] == '-') && continue
        compared += 1
        seq1[i] == seq2[i] && (matches += 1)
    end
    return compared > 0 ? matches / compared : 0.0
end

# consensus = most frequent standard residue per column among stored seqs
function consensus_sequence(stored::Vector{String})
    L = maximum(length.(stored))
    cons = Char[]
    for pos in 1:L
        counts = zeros(Int, length(AA_ALPHABET))
        for s in stored
            pos > length(s) && continue
            idx = get(AA_TO_IDX, s[pos], 0)
            idx > 0 && (counts[idx] += 1)
        end
        push!(cons, sum(counts) == 0 ? '-' : AA_ALPHABET[argmax(counts)])
    end
    return String(cons)
end

qstr(x) = @sprintf("%.3f", x)

println("\nFamily        K     L | NN-band(stored): mean  sd  p10  p50  p90  min  max | gen: pairID  Hamming  nDistinct | consID mean")
println("-"^150)

band_rows = String[]
push!(band_rows, "Family,K,L,NN_mean,NN_sd,NN_p05,NN_p10,NN_p25,NN_p50,NN_p75,NN_p90,NN_p95,NN_min,NN_max,gen_meanPairID,gen_meanHamming,gen_nDistinct,gen_consID_mean,gen_consID_sd")

consid_rows = String[]
push!(consid_rows, "Family,SeqIdx,ConsensusID")

for fam in FAMILIES
    fam_dir = joinpath(PFAM_DATA_DIR, fam.id)
    stored = [s[2] for s in parse_fasta(joinpath(fam_dir, "stored_sequences.fasta"))]
    gen    = [s[2] for s in parse_fasta(joinpath(fam_dir, "sa_generation_sequences.fasta"))]
    K = length(stored); L = length(stored[1]); S = length(gen)

    # (1) leave-one-out nearest-neighbor identity per stored sequence
    nn = Float64[]
    for i in 1:K
        best = 0.0
        for j in 1:K
            i == j && continue
            best = max(best, sequence_identity(stored[i], stored[j]))
        end
        push!(nn, best)
    end
    qs = quantile(nn, [0.05,0.10,0.25,0.50,0.75,0.90,0.95])

    # (2) generated diversity (pairwise among generated)
    pid = Float64[]; ham = Float64[]
    for i in 1:(S-1), j in (i+1):S
        id = sequence_identity(gen[i], gen[j])
        push!(pid, id); push!(ham, 1 - id)
    end
    n_distinct = length(unique(gen))

    # (3) consensus identity per generated sequence
    cons = consensus_sequence(stored)
    consid = [sequence_identity(g, cons) for g in gen]
    for (i, c) in enumerate(consid)
        push!(consid_rows, "$(fam.name),$i,$(qstr(c))")
    end

    @printf("%-12s %4d %4d | %5s %5s %5s %5s %5s %5s %5s | %7s %7s %6d | %6s\n",
        fam.name, K, L,
        qstr(mean(nn)), qstr(std(nn)), qstr(qs[2]), qstr(qs[4]), qstr(qs[6]),
        qstr(minimum(nn)), qstr(maximum(nn)),
        qstr(mean(pid)), qstr(mean(ham)), n_distinct, qstr(mean(consid)))

    push!(band_rows, join([fam.name, K, L,
        qstr(mean(nn)), qstr(std(nn)), qstr(qs[1]), qstr(qs[2]), qstr(qs[3]),
        qstr(qs[4]), qstr(qs[5]), qstr(qs[6]), qstr(qs[7]),
        qstr(minimum(nn)), qstr(maximum(nn)),
        qstr(mean(pid)), qstr(mean(ham)), n_distinct,
        qstr(mean(consid)), qstr(std(consid))], ","))
end

write(joinpath(OUT_DIR, "review_identity_band.csv"), join(band_rows, "\n") * "\n")
write(joinpath(OUT_DIR, "review_consensus_identity.csv"), join(consid_rows, "\n") * "\n")
println("\nWrote review_identity_band.csv and review_consensus_identity.csv to $OUT_DIR")
