#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Near-Duplicate / Identity Distribution Analysis
#
# For each family, compute the full pairwise identity between every generated
# sequence and every stored sequence. Report:
#   - fraction of generated seqs with >95%, >90%, >80% identity to any stored
#   - fraction within 1, 2, 5 substitutions of a stored sequence
#   - compare against within-family stored-vs-stored identity distribution
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment ..."
const _EXPERIMENT_ROOT = @__DIR__
include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))

const ANALYSIS_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "data")
const PFAM_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "..", "pfam-families", "data")
mkpath(ANALYSIS_DATA_DIR)

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

# ══════════════════════════════════════════════════════════════════════════════
# Compute number of substitutions between two sequences
# ══════════════════════════════════════════════════════════════════════════════

function count_substitutions(seq1::String, seq2::String)
    L = min(length(seq1), length(seq2))
    subs = 0
    compared = 0
    for i in 1:L
        (seq1[i] == '-' || seq2[i] == '-') && continue
        compared += 1
        seq1[i] != seq2[i] && (subs += 1)
    end
    return subs, compared
end

# ══════════════════════════════════════════════════════════════════════════════
# Main analysis
# ══════════════════════════════════════════════════════════════════════════════

@info "Computing identity distributions ..."

summary_rows = Vector{NamedTuple}()
detail_rows  = Vector{NamedTuple}()

for fam in FAMILIES
    @info "  $(fam.name) ($(fam.id)) ..."
    fam_dir = joinpath(PFAM_DATA_DIR, fam.id)

    # load stored sequences
    stored_data = parse_fasta(joinpath(fam_dir, "stored_sequences.fasta"))
    stored_seqs = [s[2] for s in stored_data]
    K = length(stored_seqs)
    L = length(stored_seqs[1])

    # load generated sequences
    gen_data = parse_fasta(joinpath(fam_dir, "sa_generation_sequences.fasta"))
    gen_seqs = [s[2] for s in gen_data]
    S = length(gen_seqs)

    # ── Generated vs stored: 150 x K identity matrix ──
    max_ids   = Float64[]   # max identity for each generated seq
    min_subs  = Int[]       # min substitutions for each generated seq

    for (i, gseq) in enumerate(gen_seqs)
        best_id = 0.0
        best_subs = typemax(Int)
        for sseq in stored_seqs
            id_val = sequence_identity(gseq, sseq)
            nsubs, _ = count_substitutions(gseq, sseq)
            best_id = max(best_id, id_val)
            best_subs = min(best_subs, nsubs)
        end
        push!(max_ids, best_id)
        push!(min_subs, best_subs)

        push!(detail_rows, (
            Family    = fam.name,
            SeqIdx    = i,
            MaxID     = best_id,
            MinSubs   = best_subs,
            SeqLen    = count(c -> c != '-', gseq),
        ))
    end

    # ── Within-family stored-vs-stored identity ──
    stored_ids = Float64[]
    for i in 1:(K-1)
        for j in (i+1):K
            push!(stored_ids, sequence_identity(stored_seqs[i], stored_seqs[j]))
        end
    end

    # ── Summary statistics ──
    frac_gt95  = count(x -> x > 0.95, max_ids) / S
    frac_gt90  = count(x -> x > 0.90, max_ids) / S
    frac_gt80  = count(x -> x > 0.80, max_ids) / S
    frac_1sub  = count(x -> x <= 1, min_subs) / S
    frac_2sub  = count(x -> x <= 2, min_subs) / S
    frac_5sub  = count(x -> x <= 5, min_subs) / S

    push!(summary_rows, (
        Family           = fam.name,
        K                = K,
        S                = S,
        L                = L,
        MeanMaxID        = mean(max_ids),
        StdMaxID         = std(max_ids),
        MedianMaxID      = median(max_ids),
        FracGt95         = frac_gt95,
        FracGt90         = frac_gt90,
        FracGt80         = frac_gt80,
        Frac1Sub         = frac_1sub,
        Frac2Sub         = frac_2sub,
        Frac5Sub         = frac_5sub,
        MeanStoredID     = mean(stored_ids),
        StdStoredID      = std(stored_ids),
        MedianStoredID   = median(stored_ids),
    ))

    @info "    Max ID: mean=$(round(mean(max_ids), digits=3)), " *
          ">95%: $(round(frac_gt95*100, digits=1))%, " *
          ">90%: $(round(frac_gt90*100, digits=1))%, " *
          "<=2 subs: $(round(frac_2sub*100, digits=1))%"
    @info "    Stored pairwise ID: mean=$(round(mean(stored_ids), digits=3))"
end

# ══════════════════════════════════════════════════════════════════════════════
# Save results
# ══════════════════════════════════════════════════════════════════════════════

summary_df = DataFrame(summary_rows)
detail_df  = DataFrame(detail_rows)

CSV.write(joinpath(ANALYSIS_DATA_DIR, "identity_distribution_summary.csv"), summary_df)
CSV.write(joinpath(ANALYSIS_DATA_DIR, "identity_distribution_detail.csv"), detail_df)

# ══════════════════════════════════════════════════════════════════════════════
# Print LaTeX table
# ══════════════════════════════════════════════════════════════════════════════

println("\n======================================================")
println("IDENTITY DISTRIBUTION SUMMARY")
println("======================================================")
println()
println("% LaTeX table: near-duplicate analysis")
println("\\begin{tabular}{@{}lrrrrrrr@{}}")
println("\\toprule")
println("Family & Mean MaxID & \\%>95\\% & \\%>90\\% & \\%>80\\% & \\%\\leq2sub & Mean Stored ID \\\\")
println("\\midrule")
for r in eachrow(summary_df)
    println("$(r.Family) & $(round(r.MeanMaxID, digits=3)) & " *
            "$(round(r.FracGt95*100, digits=1)) & " *
            "$(round(r.FracGt90*100, digits=1)) & " *
            "$(round(r.FracGt80*100, digits=1)) & " *
            "$(round(r.Frac2Sub*100, digits=1)) & " *
            "$(round(r.MeanStoredID, digits=3)) \\\\")
end
println("\\bottomrule")
println("\\end{tabular}")

@info "Done. Results saved to identity_distribution_*.csv"
