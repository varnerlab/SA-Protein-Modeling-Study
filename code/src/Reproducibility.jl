# ──────────────────────────────────────────────────────────────────────────────
# Reproducibility.jl — Utilities for reproducible experiment execution
#
# Provides: explicit RNG factory, run metadata logging, input validation.
# ──────────────────────────────────────────────────────────────────────────────

"""
    make_experiment_rng(seed::Int) -> AbstractRNG

Create an explicit, local random number generator seeded with `seed`.
Uses the same RNG type as Julia's default (Xoshiro on Julia 1.7+) so that
the output stream is identical to `Random.seed!(seed); rand(...)`.

# Example
```julia
rng = make_experiment_rng(42)
x = rand(rng, 10)            # does not affect global RNG
idx = StatsBase.sample(rng, 1:100, 5)
```
"""
function make_experiment_rng(seed::Int)
    rng = copy(Random.default_rng())
    Random.seed!(rng, seed)
    return rng
end

"""
    save_run_metadata(output_dir; seed=nothing, kwargs...)

Write a JSON sidecar file (`run_metadata.json`) recording the seed, Julia version,
git revision, timestamp, and any additional keyword arguments.
"""
function save_run_metadata(output_dir::String; seed=nothing, kwargs...)
    mkpath(output_dir)

    # try to capture git revision
    git_rev = try
        strip(read(`git -C $output_dir rev-parse --short HEAD`, String))
    catch
        "unknown"
    end

    metadata = Dict{String,Any}(
        "seed"           => seed,
        "julia_version"  => string(VERSION),
        "git_revision"   => git_rev,
        "timestamp"      => Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS"),
    )
    for (k, v) in kwargs
        metadata[string(k)] = v
    end

    path = joinpath(output_dir, "run_metadata.json")
    open(path, "w") do io
        # simple JSON serialization (no external dependency)
        println(io, "{")
        entries = collect(pairs(metadata))
        for (i, (k, v)) in enumerate(entries)
            val_str = v isa AbstractString ? "\"$(v)\"" :
                      v === nothing ? "null" : string(v)
            comma = i < length(entries) ? "," : ""
            println(io, "  \"$k\": $val_str$comma")
        end
        println(io, "}")
    end
    @info "  Run metadata saved to $path"
    return path
end

"""
    assert_inputs_exist(files; context="")

Check that every path in `files` exists on disk. Throws an informative error
listing all missing files if any are absent.

# Example
```julia
assert_inputs_exist([
    joinpath(data_dir, "stored_sequences.fasta"),
    joinpath(data_dir, "PF00076_seed.sto"),
]; context="run_pfam_families")
```
"""
function assert_inputs_exist(files; context::String="")
    missing_files = filter(f -> !ispath(f), files)
    if !isempty(missing_files)
        ctx = isempty(context) ? "" : " [$context]"
        msg = "Missing required input files$ctx:\n" *
              join(["  - $f" for f in missing_files], "\n")
        error(msg)
    end
end
