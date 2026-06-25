#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Pfam Family Size Census
#
# Downloads the complete Pfam-A seed and full alignment files (release 36.0)
# and counts the number of sequences per family. Outputs a CSV with columns:
#   accession, family_id, n_seed, n_full
# plus summary statistics used in the paper.
#
# Data source: https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam36.0/
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment …"
const _EXPERIMENT_ROOT = @__DIR__
include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))

const CENSUS_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "data")
const FIG_DIR = joinpath(dirname(_EXPERIMENT_ROOT), "..", "arxiv", "figs")
mkpath(CENSUS_DATA_DIR)
@info "Environment loaded."

# ══════════════════════════════════════════════════════════════════════════════
# Configuration
# ══════════════════════════════════════════════════════════════════════════════

const PFAM_VERSION = "Pfam36.0"
const FTP_BASE = "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/$(PFAM_VERSION)"
const SEED_URL = "$(FTP_BASE)/Pfam-A.seed.gz"
const FULL_URL = "$(FTP_BASE)/Pfam-A.full.gz"

# ══════════════════════════════════════════════════════════════════════════════
# Download helpers
# ══════════════════════════════════════════════════════════════════════════════

"""
    download_if_needed(url, dest) -> String

Download `url` to `dest` if not already present. Returns `dest`.
"""
function download_if_needed(url::String, dest::String)
    if isfile(dest)
        @info "  Already downloaded: $dest"
        return dest
    end
    @info "  Downloading $(basename(dest)) from $(PFAM_VERSION) …"
    @info "  URL: $url"
    @info "  This may take a while for large files."
    Downloads.download(url, dest)
    @info "  Saved to $dest ($(round(filesize(dest) / 1e9, digits=2)) GB)"
    return dest
end

# ══════════════════════════════════════════════════════════════════════════════
# Stockholm stream parser (counts only, no sequence storage)
# ══════════════════════════════════════════════════════════════════════════════

"""
    count_sequences_per_family(gz_path) -> DataFrame

Stream-parse a gzipped multi-entry Stockholm alignment file.
For each family entry (delimited by `//`), extract the accession,
family ID, and number of unique sequences. Returns a DataFrame.
"""
function count_sequences_per_family(gz_path::String)
    accessions = String[]
    family_ids = String[]
    seq_counts = Int[]

    current_acc = ""
    current_id  = ""
    current_seqs = Set{String}()
    n_families = 0

    # stream decompress via gunzip pipe
    @info "  Parsing $(basename(gz_path)) …"
    io = open(`gunzip -c $gz_path`)
    try
        for line in eachline(io)
            if startswith(line, "#=GF AC")
                # e.g., "#=GF AC   PF00001.23" — strip version suffix
                acc_raw = strip(line[8:end])
                current_acc = split(acc_raw, '.')[1]
            elseif startswith(line, "#=GF ID")
                current_id = strip(line[8:end])
            elseif startswith(line, "//")
                # end of entry — save counts
                if !isempty(current_acc)
                    push!(accessions, current_acc)
                    push!(family_ids, current_id)
                    push!(seq_counts, length(current_seqs))
                    n_families += 1
                    if n_families % 2000 == 0
                        @info "    Processed $n_families families …"
                    end
                end
                # reset for next entry
                current_acc = ""
                current_id  = ""
                empty!(current_seqs)
            elseif !startswith(line, "#") && !startswith(line, "//")
                stripped = strip(line)
                if !isempty(stripped)
                    parts = split(stripped)
                    if length(parts) >= 2
                        push!(current_seqs, parts[1])
                    end
                end
            end
        end
    finally
        close(io)
    end

    @info "  Done: $n_families families parsed from $(basename(gz_path))"
    return DataFrame(accession=accessions, family_id=family_ids, n_seq=seq_counts)
end

# ══════════════════════════════════════════════════════════════════════════════
# Summary statistics
# ══════════════════════════════════════════════════════════════════════════════

function print_summary(df::DataFrame, label::String)
    n = df.n_seq
    println("\n", "="^60)
    println("  $label alignment census ($PFAM_VERSION)")
    println("="^60)
    println("  Total families:       $(nrow(df))")
    println("  Median sequences:     $(median(n))")
    println("  Mean sequences:       $(round(mean(n), digits=1))")
    println("  Q25:                  $(quantile(n, 0.25))")
    println("  Q75:                  $(quantile(n, 0.75))")
    println("  Min:                  $(minimum(n))")
    println("  Max:                  $(maximum(n))")
    println()
    for threshold in [10, 20, 50, 100, 500, 1000]
        frac = mean(n .< threshold)
        count = sum(n .< threshold)
        println("  Families with < $(lpad(threshold, 4)) sequences: " *
                "$(lpad(count, 6)) ($(round(100*frac, digits=1))%)")
    end
    println("="^60)
end

# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

function main()
    # download files
    seed_gz = download_if_needed(SEED_URL, joinpath(CENSUS_DATA_DIR, "Pfam-A.seed.gz"))
    full_gz = download_if_needed(FULL_URL, joinpath(CENSUS_DATA_DIR, "Pfam-A.full.gz"))

    # parse seed alignments
    @info "Parsing seed alignments …"
    df_seed = count_sequences_per_family(seed_gz)
    rename!(df_seed, :n_seq => :n_seed)
    print_summary(rename(df_seed, :n_seed => :n_seq), "Seed")

    # parse full alignments
    @info "Parsing full alignments …"
    df_full = count_sequences_per_family(full_gz)
    rename!(df_full, :n_seq => :n_full)
    print_summary(rename(df_full, :n_full => :n_seq), "Full")

    # merge on accession
    df = outerjoin(df_seed, df_full[:, [:accession, :n_full]], on=:accession)
    sort!(df, :n_seed)

    # save
    out_path = joinpath(CENSUS_DATA_DIR, "pfam_family_sizes.csv")
    CSV.write(out_path, df)
    @info "Saved census data to $out_path ($(nrow(df)) families)"

    return df
end

df = main()
