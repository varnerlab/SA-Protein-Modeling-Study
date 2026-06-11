#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# Structure Validation: ESMFold predictions + pLDDT & TM-score analysis
#
# Predicts structures for SA-generated, SA-retrieval, bootstrap, and stored
# sequences using the ESMFold REST API, then evaluates structural plausibility
# via pLDDT confidence and TM-score to a family reference structure.
# ──────────────────────────────────────────────────────────────────────────────

@info "Loading environment …"
const _EXPERIMENT_ROOT = @__DIR__

# Activate the project environment from the parent code/ directory
using Pkg
Pkg.activate(joinpath(_EXPERIMENT_ROOT, ".."))

include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))

using HTTP
using Printf

const PFAM_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "..", "pfam-families", "data")
const STRUCT_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "data")
const RESULTS_DIR = joinpath(_EXPERIMENT_ROOT, "results")
const FIGS_DIR = joinpath(_EXPERIMENT_ROOT, "figs")
const TMALIGN_BIN = joinpath(_EXPERIMENT_ROOT, "bin", "TMalign")

const FAMILIES = [
    (id="PF00076", name="RRM",     pdb="1FXL", chain="A"),
    (id="PF00018", name="SH3",     pdb="1SHG", chain="A"),
    (id="PF00397", name="WW",      pdb="1PIN", chain="A"),
    (id="PF00014", name="Kunitz",  pdb="1BPI", chain="A"),
    (id="PF00096", name="zf-C2H2", pdb="1ZAA", chain="C"),
    (id="PF00595", name="PDZ",     pdb="1BE9", chain="A"),
    (id="PF00069", name="Pkinase", pdb="1ATP", chain="E"),
    (id="PF00711", name="Defensin_beta", pdb="1E4S", chain="A"),
]

const ESMFOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
const RATE_LIMIT_SLEEP = 1.0  # seconds between API calls
const MAX_SEQ_FOR_STRUCTURE = 50  # predict structures for first N per method (file order; no fold-quality selection)

# ══════════════════════════════════════════════════════════════════════════════
# ESMFold API
# ══════════════════════════════════════════════════════════════════════════════

function esmfold_predict(sequence::String; retries=3, timeout=120)
    for attempt in 1:retries
        try
            response = HTTP.post(ESMFOLD_URL, [], sequence;
                                 readtimeout=timeout, connect_timeout=30,
                                 retry=false)
            return String(response.body)
        catch e
            if attempt < retries
                wait_time = 2^attempt
                @warn "  ESMFold attempt $attempt failed: $e — retrying in $(wait_time)s"
                sleep(wait_time)
            else
                @error "  ESMFold failed after $retries attempts: $e"
                return nothing
            end
        end
    end
    return nothing
end

function extract_mean_plddt(pdb_string::String)
    plddts = Float64[]
    for line in split(pdb_string, '\n')
        if startswith(line, "ATOM") && length(line) >= 66
            try
                plddt = parse(Float64, strip(line[61:66]))
                push!(plddts, plddt)
            catch
                continue
            end
        end
    end
    if isempty(plddts)
        return 0.0
    end
    m = mean(plddts)
    # ESMFold API returns pLDDT in [0,1]; convert to [0,100] if needed
    return m < 1.5 ? m * 100.0 : m
end

# ══════════════════════════════════════════════════════════════════════════════
# Reference structure download from RCSB
# ══════════════════════════════════════════════════════════════════════════════

function download_reference_pdb(pdb_id::String, chain::String, outdir::String)
    pdb_file = joinpath(outdir, "$(pdb_id)_$(chain).pdb")
    if isfile(pdb_file)
        @info "  Reference PDB already exists: $pdb_file"
        return pdb_file
    end

    url = "https://files.rcsb.org/download/$(uppercase(pdb_id)).pdb"
    @info "  Downloading reference PDB: $url"
    full_pdb = Downloads.download(url)
    full_content = read(full_pdb, String)

    # Extract the specified chain
    chain_lines = String[]
    for line in split(full_content, '\n')
        if (startswith(line, "ATOM") || startswith(line, "HETATM")) && length(line) >= 22
            line_chain = string(line[22])
            if line_chain == chain && startswith(line, "ATOM")
                push!(chain_lines, line)
            end
        end
    end
    push!(chain_lines, "END")

    mkpath(outdir)
    open(pdb_file, "w") do f
        for line in chain_lines
            println(f, line)
        end
    end
    @info "  Saved reference: $pdb_file ($(length(chain_lines)-1) ATOM records)"
    return pdb_file
end

# ══════════════════════════════════════════════════════════════════════════════
# TM-align
# ══════════════════════════════════════════════════════════════════════════════

function run_tmalign(query_pdb::String, reference_pdb::String)
    if !isfile(TMALIGN_BIN)
        @error "TMalign binary not found at $TMALIGN_BIN"
        return nothing
    end

    try
        output = read(`$TMALIGN_BIN $query_pdb $reference_pdb`, String)
        # Parse TM-score (normalized by reference length)
        for line in split(output, '\n')
            if contains(line, "TM-score") && contains(line, "Chain_2")
                m = match(r"TM-score=\s*([\d.]+)", line)
                if m !== nothing
                    return parse(Float64, m.captures[1])
                end
            end
        end
        # Fallback: take first TM-score line
        for line in split(output, '\n')
            if contains(line, "TM-score")
                m = match(r"TM-score=\s*([\d.]+)", line)
                if m !== nothing
                    return parse(Float64, m.captures[1])
                end
            end
        end
    catch e
        @warn "  TMalign failed: $e"
    end
    return nothing
end

# ══════════════════════════════════════════════════════════════════════════════
# Predict structures for a set of sequences
# ══════════════════════════════════════════════════════════════════════════════

function predict_structures(sequences::Vector{Tuple{String,String}},
                           method_name::String,
                           struct_dir::String;
                           max_seqs=MAX_SEQ_FOR_STRUCTURE)
    mkpath(struct_dir)
    n = min(length(sequences), max_seqs)
    plddts = Float64[]
    pdb_files = String[]

    for i in 1:n
        name, seq = sequences[i]
        # Strip gaps and non-standard characters
        clean_seq = filter(c -> c in AA_ALPHABET, uppercase(seq))
        if isempty(clean_seq)
            @warn "  Skipping empty sequence: $name"
            continue
        end

        # Sanitize filename: replace / and other problematic chars
        safe_name = replace(name, "/" => "_", "\\" => "_", " " => "_")
        pdb_file = joinpath(struct_dir, "$(safe_name).pdb")

        # Check if already predicted (checkpoint — reuse existing PDB)
        if isfile(pdb_file)
            pdb_content = read(pdb_file, String)
            plddt = extract_mean_plddt(pdb_content)
            push!(plddts, plddt)
            push!(pdb_files, pdb_file)
            if i % 10 == 0
                @info "  $method_name: $i/$n (cached, pLDDT=$(round(plddt, digits=1)))"
            end
            continue
        end

        # Call ESMFold API
        pdb_content = esmfold_predict(clean_seq)
        if pdb_content === nothing
            @warn "  Failed to predict structure for $name"
            continue
        end

        # Save PDB
        open(pdb_file, "w") do f
            write(f, pdb_content)
        end

        plddt = extract_mean_plddt(pdb_content)
        push!(plddts, plddt)
        push!(pdb_files, pdb_file)

        if i % 10 == 0
            @info "  $method_name: $i/$n predicted (pLDDT=$(round(plddt, digits=1)))"
        end

        sleep(RATE_LIMIT_SLEEP)
    end

    return plddts, pdb_files
end

# ══════════════════════════════════════════════════════════════════════════════
# Generate bootstrap sequences (resample from stored)
# ══════════════════════════════════════════════════════════════════════════════

function generate_bootstrap_sequences(stored_data, n=MAX_SEQ_FOR_STRUCTURE)
    rng = Random.MersenneTwister(42)
    indices = rand(rng, 1:length(stored_data), n)
    return [("bootstrap_$i", stored_data[idx][2]) for (i, idx) in enumerate(indices)]
end

# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

function run_structure_validation()
    mkpath(RESULTS_DIR)
    mkpath(FIGS_DIR)

    all_rows = Vector{NamedTuple}()

    for fam in FAMILIES
        @info "════════════════════════════════════════════════════════════"
        @info "  $(fam.name) ($(fam.id))"
        @info "════════════════════════════════════════════════════════════"

        fam_struct_dir = joinpath(STRUCT_DATA_DIR, fam.id)
        mkpath(fam_struct_dir)

        # ── Load sequences ──────────────────────────────────────────────
        sa_gen_fasta = joinpath(PFAM_DATA_DIR, fam.id, "sa_generation_sequences.fasta")
        sa_ret_fasta = joinpath(PFAM_DATA_DIR, fam.id, "sa_retrieval_sequences.fasta")
        stored_fasta = joinpath(PFAM_DATA_DIR, fam.id, "stored_sequences.fasta")

        if !isfile(sa_gen_fasta) || !isfile(stored_fasta)
            @warn "  Missing FASTA files for $(fam.id), skipping"
            continue
        end

        sa_gen_data  = parse_fasta(sa_gen_fasta)
        sa_ret_data  = isfile(sa_ret_fasta) ? parse_fasta(sa_ret_fasta) : Tuple{String,String}[]
        stored_data  = parse_fasta(stored_fasta)
        bootstrap_data = generate_bootstrap_sequences(stored_data)

        @info "  Loaded: $(length(sa_gen_data)) SA gen, $(length(sa_ret_data)) SA ret, $(length(stored_data)) stored"

        # ── Download reference PDB ──────────────────────────────────────
        ref_dir = joinpath(fam_struct_dir, "reference")
        ref_pdb = download_reference_pdb(fam.pdb, fam.chain, ref_dir)

        # ── Predict structures ──────────────────────────────────────────
        methods = [
            ("SA generation", sa_gen_data,    joinpath(fam_struct_dir, "structures", "sa_generation")),
            ("SA retrieval",  sa_ret_data,    joinpath(fam_struct_dir, "structures", "sa_retrieval")),
            ("Bootstrap",     bootstrap_data, joinpath(fam_struct_dir, "structures", "bootstrap")),
            ("Stored",        stored_data,    joinpath(fam_struct_dir, "structures", "stored")),
        ]

        for (method_name, seqs, struct_dir) in methods
            if isempty(seqs)
                @warn "  No sequences for $method_name, skipping"
                continue
            end

            @info "  Predicting structures: $method_name ($(min(length(seqs), MAX_SEQ_FOR_STRUCTURE)) sequences) …"
            plddts, pdb_files = predict_structures(
                [(s[1], s[2]) for s in seqs],
                method_name, struct_dir
            )

            if isempty(plddts)
                @warn "  No structures predicted for $method_name"
                continue
            end

            # ── TM-scores ───────────────────────────────────────────────
            @info "  Computing TM-scores for $method_name …"
            tm_scores = Float64[]
            for pdb_file in pdb_files
                tm = run_tmalign(pdb_file, ref_pdb)
                if tm !== nothing
                    push!(tm_scores, tm)
                end
            end

            # ── Summary stats ───────────────────────────────────────────
            n_seq = length(plddts)
            plddt_mean = mean(plddts)
            plddt_se = n_seq > 1 ? std(plddts) / sqrt(n_seq) : 0.0
            frac_plddt70 = mean(p -> p > 70 ? 1.0 : 0.0, plddts)

            tm_mean = isempty(tm_scores) ? NaN : mean(tm_scores)
            tm_se = length(tm_scores) > 1 ? std(tm_scores) / sqrt(length(tm_scores)) : 0.0
            frac_tm05 = isempty(tm_scores) ? NaN : mean(t -> t > 0.5 ? 1.0 : 0.0, tm_scores)

            @info @sprintf("  %s: pLDDT=%.1f±%.1f (%.0f%% >70), TM=%.3f±%.3f (%.0f%% >0.5)",
                method_name, plddt_mean, plddt_se, 100*frac_plddt70,
                tm_mean, tm_se, 100*frac_tm05)

            push!(all_rows, (
                Family      = fam.name,
                PfamID      = fam.id,
                Method      = method_name,
                N           = n_seq,
                pLDDT_mean  = plddt_mean,
                pLDDT_se    = plddt_se,
                pLDDT_gt70  = frac_plddt70,
                TM_mean     = tm_mean,
                TM_se       = tm_se,
                TM_gt05     = frac_tm05,
            ))
        end

        # ── Save per-family CSV ─────────────────────────────────────────
        fam_rows = filter(r -> r.PfamID == fam.id, all_rows)
        if !isempty(fam_rows)
            fam_df = DataFrame(fam_rows)
            CSV.write(joinpath(fam_struct_dir, "structure_metrics.csv"), fam_df)
        end
    end

    # ══════════════════════════════════════════════════════════════════════
    # Aggregate results
    # ══════════════════════════════════════════════════════════════════════

    if isempty(all_rows)
        @error "No results collected!"
        return nothing
    end

    df = DataFrame(all_rows)
    CSV.write(joinpath(RESULTS_DIR, "structure_validation_results.csv"), df)
    @info "Results saved to $(joinpath(RESULTS_DIR, "structure_validation_results.csv"))"

    pretty_table(df[:, [:Family, :Method, :pLDDT_mean, :pLDDT_se, :pLDDT_gt70, :TM_mean, :TM_se, :TM_gt05]])

    # ══════════════════════════════════════════════════════════════════════
    # Generate figures
    # ══════════════════════════════════════════════════════════════════════
    generate_structure_figures(df)

    return df
end

# ══════════════════════════════════════════════════════════════════════════════
# Figures
# ══════════════════════════════════════════════════════════════════════════════

function generate_structure_figures(df::DataFrame)
    families = unique(df.Family)
    methods = ["SA generation", "SA retrieval", "Bootstrap", "Stored"]
    method_colors = [
        RGB(0.122, 0.467, 0.706),  # SA gen - blue
        RGB(1.0, 0.498, 0.055),    # SA ret - orange
        RGB(0.839, 0.153, 0.157),  # Bootstrap - red
        RGB(0.173, 0.627, 0.173),  # Stored - green
    ]

    n_fam = length(families)
    n_meth = length(methods)
    bar_width = 0.18
    offsets = [-1.5, -0.5, 0.5, 1.5] .* bar_width

    # ── Panel 1: pLDDT ──────────────────────────────────────────────────
    p1 = plot(size=(500, 300), ylabel="Mean pLDDT",
              xticks=(1:n_fam, families), xrotation=30,
              legend=:topright, legendfontsize=7,
              ylims=(0, 100), grid=:y,
              left_margin=5Plots.mm, bottom_margin=8Plots.mm)
    hline!(p1, [70], color=:gray, linestyle=:dash, linewidth=0.8, label="pLDDT=70")

    for (j, method) in enumerate(methods)
        xs = Float64[]
        ys = Float64[]
        es = Float64[]
        for (i, fam) in enumerate(families)
            row = filter(r -> r.Family == fam && r.Method == method, df)
            if nrow(row) > 0
                push!(xs, i + offsets[j])
                push!(ys, row.pLDDT_mean[1])
                push!(es, row.pLDDT_se[1])
            end
        end
        bar!(p1, xs, ys, yerror=es, bar_width=bar_width,
             color=method_colors[j], label=method, linewidth=0.5)
    end

    # ── Panel 2: TM-score ───────────────────────────────────────────────
    p2 = plot(size=(500, 300), ylabel="TM-score",
              xticks=(1:n_fam, families), xrotation=30,
              legend=:topright, legendfontsize=7,
              ylims=(0, 1.0), grid=:y,
              left_margin=5Plots.mm, bottom_margin=8Plots.mm)
    hline!(p2, [0.5], color=:gray, linestyle=:dash, linewidth=0.8, label="TM=0.5")

    for (j, method) in enumerate(methods)
        xs = Float64[]
        ys = Float64[]
        es = Float64[]
        for (i, fam) in enumerate(families)
            row = filter(r -> r.Family == fam && r.Method == method, df)
            if nrow(row) > 0 && !isnan(row.TM_mean[1])
                push!(xs, i + offsets[j])
                push!(ys, row.TM_mean[1])
                push!(es, row.TM_se[1])
            end
        end
        bar!(p2, xs, ys, yerror=es, bar_width=bar_width,
             color=method_colors[j], label=method, linewidth=0.5)
    end

    # ── Composite ───────────────────────────────────────────────────────
    composite = plot(p1, p2, layout=(1, 2), size=(1000, 350),
                     title=["(A) ESMFold Confidence (pLDDT)" "(B) Structural Similarity (TM-score)"],
                     titlefontsize=10)

    savefig(composite, joinpath(FIGS_DIR, "structure_validation_composite.pdf"))
    @info "Figure saved: $(joinpath(FIGS_DIR, "structure_validation_composite.pdf"))"

    # Copy to paper figs directory
    paper_figs = joinpath(_EXPERIMENT_ROOT, "..", "..", "paper", "figs")
    if isdir(paper_figs)
        cp(joinpath(FIGS_DIR, "structure_validation_composite.pdf"),
           joinpath(paper_figs, "structure_validation_composite.pdf"), force=true)
        @info "Copied to paper/figs/"
    end

    return composite
end

# ── Run ──────────────────────────────────────────────────────────────────────
df_struct = run_structure_validation()
