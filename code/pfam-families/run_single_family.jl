#!/usr/bin/env julia
# Run SA experiment for a single family (avoids re-running all families)
# Usage: julia --project=.. run_single_family.jl PF00711

length(ARGS) == 0 && error("Usage: julia --project=.. run_single_family.jl <PFAM_ID>")
const TARGET_ID = ARGS[1]

@info "Loading environment …"
const _EXPERIMENT_ROOT = @__DIR__
include(joinpath(_EXPERIMENT_ROOT, "..", "Include.jl"))

const PFAM_DATA_DIR = joinpath(_EXPERIMENT_ROOT, "data")
const PFAM_FIG_DIR  = joinpath(_EXPERIMENT_ROOT, "figs")
@info "Environment loaded."

# Load full family list from main script (just the definitions)
include_string(Main, """
const ALL_FAMILIES = [
    (id="PF00076", name="RRM",     desc="RNA Recognition Motif"),
    (id="PF00018", name="SH3",     desc="SH3 domain"),
    (id="PF00397", name="WW",      desc="WW domain"),
    (id="PF00014", name="Kunitz",  desc="Kunitz/BPTI domain"),
    (id="PF00096", name="zf-C2H2", desc="Zinc finger C2H2"),
    (id="PF00595", name="PDZ",     desc="PDZ domain"),
    (id="PF00069", name="Pkinase", desc="Protein kinase domain"),
    (id="PF00711", name="Defensin_beta", desc="Beta defensin"),
]
""")

fam = findfirst(f -> f.id == TARGET_ID, ALL_FAMILIES)
isnothing(fam) && error("Family $TARGET_ID not found")
fam = ALL_FAMILIES[fam]

# Parameters (same as main experiment)
const α_step        = 0.01
const n_chains      = 30
const T_per_chain   = 5000
const T_burnin      = 2000
const thin_interval = 100
const samples_per_chain = 5
const S             = n_chains * samples_per_chain
const σ_init        = 0.01
const R_pca         = 0.95

# Include the SA chain functions from run_pfam_families.jl
# (re-define them here to avoid const conflicts)

function run_sa_chains_local(X̂::Matrix{Float64}, β::Float64)
    d, K = size(X̂)
    all_samples = Vector{Vector{Float64}}()
    rng_init = make_experiment_rng(42)
    pattern_indices = StatsBase.sample(rng_init, 1:K, n_chains, replace=(n_chains > K))
    for (c, k) in enumerate(pattern_indices)
        seed_c = 12345 + c
        rng_c = make_experiment_rng(seed_c)
        ξ₀ = X̂[:, k] .+ σ_init .* randn(rng_c, d)
        (_, Ξ) = sample(X̂, ξ₀, T_per_chain; β=β, α=α_step, seed=seed_c)
        chain_pool = Vector{Vector{Float64}}()
        for tᵢ in (T_burnin+1):thin_interval:T_per_chain
            push!(chain_pool, Ξ[tᵢ, :])
        end
        n_avail = length(chain_pool)
        n_avail == 0 && continue
        idxs = round.(Int, range(1, n_avail, length=min(samples_per_chain, n_avail)))
        for idx in idxs
            push!(all_samples, chain_pool[idx])
        end
    end
    return all_samples
end

function run_mala_chains_local(X̂::Matrix{Float64}, β::Float64)
    d, K = size(X̂)
    all_samples = Vector{Vector{Float64}}()
    accept_rates = Float64[]
    rng_init = make_experiment_rng(42)
    pattern_indices = StatsBase.sample(rng_init, 1:K, n_chains, replace=(n_chains > K))
    for (c, k) in enumerate(pattern_indices)
        seed_c = 12345 + c
        rng_c = make_experiment_rng(seed_c)
        ξ₀ = X̂[:, k] .+ σ_init .* randn(rng_c, d)
        (_, Ξ, ar) = mala_sample(X̂, ξ₀, T_per_chain; β=β, α=α_step, seed=seed_c)
        push!(accept_rates, ar)
        chain_pool = Vector{Vector{Float64}}()
        for tᵢ in (T_burnin+1):thin_interval:T_per_chain
            push!(chain_pool, Ξ[tᵢ, :])
        end
        n_avail = length(chain_pool)
        n_avail == 0 && continue
        idxs = round.(Int, range(1, n_avail, length=min(samples_per_chain, n_avail)))
        for idx in idxs
            push!(all_samples, chain_pool[idx])
        end
    end
    return all_samples, mean(accept_rates)
end

function run_simple_baselines_local(X̂::Matrix{Float64}, β::Float64, S::Int)
    d, K = size(X̂)
    σ_noise = sqrt(2 * α_step / β)
    rng_bs = make_experiment_rng(12345)
    bs = [copy(X̂[:, rand(rng_bs, 1:K)]) for _ in 1:S]
    rng_gp = make_experiment_rng(12345)
    gp = [X̂[:, rand(rng_gp, 1:K)] .+ σ_noise .* randn(rng_gp, d) for _ in 1:S]
    dirichlet = Dirichlet(K, 1.0)
    rng_rc = make_experiment_rng(12345)
    rc = [X̂ * rand(rng_rc, dirichlet) for _ in 1:S]
    return Dict("Bootstrap" => bs, "Gaussian perturbation" => gp, "Convex combination" => rc)
end

function evaluate_method_local(samps, seqs, X̂, β, stored_seqs, L)
    novelty_vals  = [sample_novelty(ξ, X̂) for ξ in samps]
    energy_vals   = [hopfield_energy(ξ, X̂, β) for ξ in samps]
    seq_id_vals   = [nearest_sequence_identity(s, stored_seqs) for s in seqs]
    return (
        Novelty   = mean(novelty_vals),
        Diversity = sample_diversity(samps),
        Energy    = mean(energy_vals),
        SeqID     = mean(seq_id_vals),
        ValidAA   = mean(valid_residue_fraction(s) for s in seqs),
        KL_AA     = aa_composition_kl(seqs, stored_seqs),
    )
end

# ── Run the single family ────────────────────────────────────────────────────

@info "Running family: $(fam.name) ($(fam.id))"

fam_data_dir = joinpath(PFAM_DATA_DIR, fam.id)
mkpath(fam_data_dir)

# Step 1: Load alignment
@info "Step 1: Loading alignment …"
sto_file = download_pfam_seed(fam.id; cache_dir=fam_data_dir)
raw_seqs = parse_stockholm(sto_file)
if isempty(raw_seqs)
    @info "  No Stockholm sequences, trying FASTA …"
    raw_seqs = parse_fasta(sto_file)
end
@info "  Parsed $(length(raw_seqs)) sequences"

# Step 2: Clean alignment
@info "Step 2: Cleaning alignment …"
char_mat, seq_names = clean_alignment(raw_seqs)
K, L = size(char_mat)
stored_seqs = [String(char_mat[k, :]) for k in 1:K]
@info "  K=$K, L=$L"

# Step 3: Build memory matrix
@info "Step 3: Building memory matrix …"
X̂, pca_model, L, d_full = build_memory_matrix(char_mat; pratio=R_pca)
d = size(X̂, 1)
@info "  d=$d (from d_full=$d_full)"

# Step 4: MSA statistics
@info "Step 4: Computing MSA statistics …"
col_ent = msa_column_entropy(char_mat)
H_col_mean = mean(col_ent)
K_eff = effective_num_sequences(char_mat; threshold=0.8)
X_onehot = onehot_encode(char_mat)
C_mat = cov(X_onehot'; dims=1)
eigvals_C = sort(eigvals(C_mat), rev=true)
λ1_ratio = eigvals_C[1] / sum(eigvals_C)
@info "  H̄_col = $(round(H_col_mean, digits=3)), K_eff = $(round(K_eff, digits=1)), λ₁/tr(C) = $(round(λ1_ratio, digits=4))"

# Step 5: Phase transition
@info "Step 5: Phase transition analysis …"
pt = find_entropy_inflection(X̂; α=α_step)
β_ret = round(Int, max(20 * pt.β_star, 50))
β_gen = round(Int, max(2 * pt.β_star, 5))
@info "  β* = $(round(pt.β_star, digits=2)), β_ret = $β_ret, β_gen = $β_gen"

# Step 6: Run SA
@info "Step 6: Running SA (retrieval β=$β_ret) …"
sa_ret_samps = run_sa_chains_local(X̂, Float64(β_ret))
@info "  SA retrieval: $(length(sa_ret_samps)) samples"

@info "Step 6b: Running SA (generation β=$β_gen) …"
sa_gen_samps = run_sa_chains_local(X̂, Float64(β_gen))
@info "  SA generation: $(length(sa_gen_samps)) samples"

# Step 7: Run MALA
@info "Step 7: Running MALA (β=$β_ret) …"
mala_samps, mala_ar = run_mala_chains_local(X̂, Float64(β_ret))
@info "  MALA: $(length(mala_samps)) samples, accept rate = $(round(mala_ar, digits=4))"

# Step 8: Simple baselines
@info "Step 8: Running simple baselines …"
baselines = run_simple_baselines_local(X̂, Float64(β_ret), S)

# Step 9: Decode and evaluate
@info "Step 9: Decoding and evaluating …"
all_rows = Vector{NamedTuple}()
methods = Dict{String, Vector{Vector{Float64}}}(
    "SA (retrieval)"  => sa_ret_samps,
    "SA (generation)" => sa_gen_samps,
    "MALA"            => mala_samps,
)
merge!(methods, baselines)

β_eval = Float64(β_ret)
for (method_name, samps) in methods
    seqs = [decode_sample(ξ, pca_model, L) for ξ in samps]
    metrics = evaluate_method_local(samps, seqs, X̂, β_eval, stored_seqs, L)
    push!(all_rows, (
        Family   = fam.name,
        PfamID   = fam.id,
        K        = K,
        L        = L,
        d_PCA    = d,
        Method   = method_name,
        β        = method_name == "SA (generation)" ? β_gen : β_ret,
        MALA_AR  = method_name == "MALA" ? mala_ar : missing,
        metrics...
    ))
end

# Step 10: Save sequences
@info "Step 10: Saving sequences …"
for (label, samps) in [("sa_retrieval", sa_ret_samps), ("sa_generation", sa_gen_samps)]
    seqs = [decode_sample(ξ, pca_model, L) for ξ in samps]
    open(joinpath(fam_data_dir, "$(label)_sequences.fasta"), "w") do io
        for (i, seq) in enumerate(seqs)
            println(io, ">$(label)_$(i)")
            println(io, seq)
        end
    end
end
open(joinpath(fam_data_dir, "stored_sequences.fasta"), "w") do io
    for (i, seq) in enumerate(stored_seqs)
        name = i <= length(seq_names) ? seq_names[i] : "stored_$i"
        println(io, ">$name")
        println(io, seq)
    end
end

# Print results
df_results = DataFrame(all_rows)
println("\n══════════════════════════════════════════════════════")
println("RESULTS: $(fam.name) ($(fam.id))")
println("══════════════════════════════════════════════════════")
pretty_table(df_results[:, [:Method, :Novelty, :Diversity, :SeqID, :KL_AA, :ValidAA]])

# Family summary
println("\nFamily summary: K=$K, L=$L, d=$d, β*=$(round(pt.β_star, digits=2)), MALA AR=$(round(mala_ar, digits=4))")

# Save individual results CSV
CSV.write(joinpath(fam_data_dir, "results.csv"), df_results)

# Also append to master summary if it exists
summary_file = joinpath(PFAM_DATA_DIR, "pfam_family_summary.csv")
if isfile(summary_file)
    df_summary = CSV.read(summary_file, DataFrame)
    # remove existing entry for this family if present
    filter!(r -> r.PfamID != fam.id, df_summary)
    push!(df_summary, (
        Family=fam.name, PfamID=fam.id, K=K, L=L, d_PCA=d, d_full=d_full,
        H_col_mean=H_col_mean, K_eff=K_eff, λ1_ratio=λ1_ratio,
        β_star=pt.β_star, β_star_theory=pt.β_star_theory,
        β_retrieval=β_ret, β_generation=β_gen,
    ); promote=true)
    CSV.write(summary_file, df_summary)
    @info "Updated $summary_file"
end

# Append to master results if it exists
results_file = joinpath(PFAM_DATA_DIR, "pfam_results.csv")
if isfile(results_file)
    df_master = CSV.read(results_file, DataFrame)
    filter!(r -> r.PfamID != fam.id, df_master)
    append!(df_master, df_results; promote=true)
    CSV.write(results_file, df_master)
    @info "Updated $results_file"
end

@info "Done with $(fam.name)!"
