# Code Review: SA-Protein-Modeling-Study Simulation Scripts

Date: 2026-03-16  
Scope: Base code and simulation scripts across sequence/structure/family/scaling/baseline/beta-prediction experiments.

## Severity-ordered findings

1. High — AF2 rankings can be misassociated with sequences in `analyze_af2_results.jl`
- In `analyze_af2_results.jl`, AF2 PDB files are collected by broad filename/rank patterns and matched by rank/order only.
- The sequence-specific lookup helper exists (`find_af2_pdb`) but is not used in the main metrics collection path.
- Impact: metrics for one sequence/family can be assigned to another, silently biasing structure validation summaries.
- Affected file: `code/structure-validation/analyze_af2_results.jl` (collect + ranking file handling logic near `collect_af2_metrics`).

2. High — Empty alignment crash path in `clean_alignment`
- `clean_alignment` assumes `raw_seqs` is non-empty and accesses `seqs[1]` unguarded.
- Impact: empty FASTA/Stockholm inputs or parser edge-cases raise hard errors and can abort runs.
- Affected file: `code/src/Protein.jl` (`clean_alignment` implementation).

3. High — Global RNG seeded at include-time
- `Random.seed!(1234)` is executed during module include in:
  - `code/Include.jl`
  - `code/sequence-experiment/Include-Sequence.jl`
- This mutates global RNG state for every include.
- Impact: reproducibility depends on include order and can be invalidated by unrelated imports.

4. Medium-High — Reproducibility drift from duplicated seed logic
- Multiple scripts define local versions of chain-running and seeding logic instead of using a central deterministic pathway.
- Impact: two experiments can use different random streams even with same declared base seed, and behavior changes if any pre-sampling code changes.
- Affected file set:
  - `code/pfam-families/run_pfam_families.jl`
  - `code/pfam-families/run_single_family.jl`
  - `code/scaling-study/run_scaling_study.jl`
  - `code/multifamily-scaling/run_multifamily_scaling.jl`
  - `code/beta-prediction/run_beta_prediction.jl`
  - `code/beta-prediction/run_beta_loocv.jl`

5. Medium — Weak integrity checks for intermediate artifacts
- Several experiment pipelines consume outputs from prior steps (alignments, stored sequences, AF2/structure outputs) without hash/version assertions.
- Impact: stale or partial artifacts can be silently consumed and produce numerically coherent but incorrect results.
- Affected file set includes:
  - `code/pfam-families/run_pfam_families.jl`
  - `code/baselines/run_hmm_baseline.jl`
  - `code/baselines/run_consensus_baseline.jl`
  - `code/structure-validation/run_structure_validation.jl`
  - `code/structure-validation/analyze_af2_results.jl`

6. Low — Implicit sequence hygiene assumptions
- Downstream pipelines often convert matrix rows to `String` without explicit validation beyond existing cleaning steps.
- Impact: malformed or unexpected symbols may propagate into metrics in edge cases.
- Affected areas: `code/src/Protein.jl` in combination with `baselines`, `scaling`, and `pfam` scripts.

## Recommended fixes (publication-facing)

1. Remove all include-time seeding (`Random.seed!`) from shared include files.
2. Use explicit RNG objects per experiment run and persist per-run seed metadata with outputs.
3. In AF2 analysis, map results by explicit `(category, seq_id)` keys and fail fast if an expected ranked model is missing.
4. Add defensive checks in `clean_alignment` for empty input and emit actionable error messages with source filename.
5. Introduce immutable input manifests (hash + source file timestamp + code version) for each experiment stage.
6. Save per-result run metadata: seed, RNG state/seed, package versions, input artifacts, and script git revision.

## Residual risks / assumptions

- Randomness in some algorithms is intentionally stochastic (e.g., bootstrap/permutation baselines); reproducibility controls should capture intended stochastic protocol, not suppress randomness.
- No runtime validation was executed for this review pass; this is a code-level review only.
