#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Prepare FASTA input files for ColabFold (AlphaFold2) structure prediction.
#
# Exports the same 50 SA-generated + 50 stored sequences per family that were
# used for ESMFold validation, in clean FASTA format suitable for
# `colabfold_batch`.
# ──────────────────────────────────────────────────────────────────────────────

const _EXPERIMENT_ROOT = @__DIR__
const PFAM_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "..", "pfam-families", "data")
const OUTPUT_DIR = joinpath(_EXPERIMENT_ROOT, "colabfold_input")
const MAX_SEQ = 50

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

const AA_ALPHABET = Set("ACDEFGHIKLMNPQRSTVWY")

function parse_fasta_simple(path::String)
    sequences = Tuple{String,String}[]
    name = ""
    seq_parts = String[]
    for line in eachline(path)
        line = strip(line)
        if startswith(line, '>')
            if !isempty(name) && !isempty(seq_parts)
                push!(sequences, (name, join(seq_parts)))
            end
            name = line[2:end]
            seq_parts = String[]
        elseif !isempty(line)
            push!(seq_parts, line)
        end
    end
    if !isempty(name) && !isempty(seq_parts)
        push!(sequences, (name, join(seq_parts)))
    end
    return sequences
end

function clean_sequence(seq::String)
    return String(filter(c -> c in AA_ALPHABET, uppercase(seq)))
end

function safe_name(name::String)
    return replace(name, "/" => "_", "\\" => "_", " " => "_", "." => "_")
end

function main()
    mkpath(OUTPUT_DIR)

    total = 0
    for fam in FAMILIES
        sa_gen_path = joinpath(PFAM_DATA_DIR, fam.id, "sa_generation_sequences.fasta")
        stored_path = joinpath(PFAM_DATA_DIR, fam.id, "stored_sequences.fasta")

        if !isfile(sa_gen_path) || !isfile(stored_path)
            @warn "Missing FASTA for $(fam.id), skipping"
            continue
        end

        sa_gen = parse_fasta_simple(sa_gen_path)
        stored = parse_fasta_simple(stored_path)

        # Take first MAX_SEQ of each (matching ESMFold experiment)
        sa_gen = sa_gen[1:min(MAX_SEQ, length(sa_gen))]
        stored = stored[1:min(MAX_SEQ, length(stored))]

        # Write per-family FASTA for ColabFold
        fam_dir = joinpath(OUTPUT_DIR, fam.id)
        mkpath(fam_dir)

        # SA generated
        open(joinpath(fam_dir, "sa_generation.fasta"), "w") do f
            for (i, (name, seq)) in enumerate(sa_gen)
                clean = clean_sequence(seq)
                sn = safe_name("$(fam.id)_sa_gen_$(i)")
                println(f, ">$sn")
                println(f, clean)
            end
        end

        # Stored
        open(joinpath(fam_dir, "stored.fasta"), "w") do f
            for (i, (name, seq)) in enumerate(stored)
                clean = clean_sequence(seq)
                sn = safe_name("$(fam.id)_stored_$(i)")
                println(f, ">$sn")
                println(f, clean)
            end
        end

        n_gen = length(sa_gen)
        n_sto = length(stored)
        total += n_gen + n_sto
        @info "$(fam.name) ($(fam.id)): $n_gen SA gen + $n_sto stored"
    end

    # Also write a single combined FASTA for convenience
    open(joinpath(OUTPUT_DIR, "all_sequences.fasta"), "w") do f
        for fam in FAMILIES
            fam_dir = joinpath(OUTPUT_DIR, fam.id)
            for fname in ["sa_generation.fasta", "stored.fasta"]
                fpath = joinpath(fam_dir, fname)
                isfile(fpath) || continue
                write(f, read(fpath, String))
            end
        end
    end

    @info "Total: $total sequences written to $OUTPUT_DIR"
    @info "Upload the colabfold_input/ directory to Google Colab and run colabfold_predict.py"
end

main()
