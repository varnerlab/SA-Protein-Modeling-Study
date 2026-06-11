# Simulated Peer Review

**Manuscript:** Training-Free Generation of Protein Sequences from Small Family Alignments via Stochastic Attention
**Author:** J. D. Varner
**Target venue:** PNAS (Research Article)
**Reviewers:** three domain experts in protein sequence generation, constraint-based / coevolutionary modeling, MCMC sampling, and uncertainty quantification.

> Note on scope: every critique below is tied to a specific passage, equation, table, or figure. Where a control or baseline already exists in the manuscript, the review acknowledges it rather than asking for it again. Suggested references are limited to work the reviewers are confident exists and is directly relevant.

---

## Reviewer 1 — Moderate

### Summary
The paper proposes stochastic attention (SA), a training-free generator that treats the modern Hopfield energy over a family alignment as a Boltzmann distribution and samples it with Langevin dynamics, using the closed-form attention residual as the score. The author evaluates SA on eight Pfam families, derives an empirical law predicting the critical inverse temperature from PCA dimensionality, compares against four baselines, and validates with two structure predictors, ESM2, and two DMS datasets. The central idea is elegant and the empirical program is unusually thorough for a method paper; my concerns are mostly about clarity, a few overclaims, and missing methodological detail rather than soundness.

### Strengths
1. **A genuinely useful idea, cleanly motivated.** The observation that the softmax attention residual is the exact score of the Hopfield Boltzmann distribution (Eq. 1–2, lines 98–106) yields a sampler with no trained parameters. The Pfam census (SI Section S20, Table S15) establishes that the small-family regime is the common case (median seed 22 sequences), so the motivation is real and quantified rather than rhetorical.
2. **Breadth of validation.** Eight families spanning a ten-fold size and seven-fold length range, two independent structure predictors (ESMFold and AlphaFold2, Fig. 4 and SI Fig. S3–S4), an independent language model (ESM2-650M, Table S7), and two DMS cross-references (Table S8) together make a far more complete case than is typical.
3. **Honest limitations section.** The Discussion (lines 93–94) candidly flags the fixed-length constraint, the selection bias toward well-structured globular domains, and the absence of wet-lab validation. This transparency is commendable.
4. **The β* prediction is practically valuable.** If the relationship β* ≈ 1.52 + 0.28√d (line 66, Fig. 2) holds up, it removes the temperature sweep and makes the method genuinely automatic from a seed alignment, which is the difference between a demo and a tool.
5. **Thoughtful negative controls.** The consensus-with-noise (Table S14) and column-permuted-alignment (Table S15) controls directly address the "is this just consensus plus noise?" objection at the level of composition and covariation.

### Weaknesses
1. **Gap handling in the encoding is unspecified.** The one-hot encoding (SI Eq. S2, line 760–764) is over the 20-letter amino acid alphabet, but seed alignments retain gaps in sequences with up to 30% gap content after cleaning (Methods, line 153). How is a gap character encoded in a 20-channel one-hot vector, and how is it decoded by argmax (which can only emit one of 20 amino acids)? This affects both the memory matrix and the decoded output length/content. Please state the gap convention explicitly and confirm whether decoded sequences can contain gaps.
2. **Reference structures for TM-align are not identified.** The structural validation (Methods, line 155; Fig. 4–5) compares predicted structures to "an experimentally determined reference," but no PDB accession codes are given for any of the eight families. TM-score depends on the reference chosen. A small SI table listing the PDB ID (and chain/residue range) used per family is needed for reproducibility.
3. **The novelty metric is measured before decoding, the identity metric after.** Novelty is defined as 1 − max_k cos(ξ̂, m_k) in PCA space (Fig. 1 caption), while sequence identity is computed on the argmax-decoded sequence. Because argmax decoding is many-to-one, the continuous-space novelty number (0.40–0.65) does not map transparently onto anything biological. Please report a sequence-space novelty (e.g., 1 − max nearest-neighbor identity) as the primary novelty measure, and keep the PCA-space cosine as a secondary diagnostic.
4. **"Training-free" and "no learned representations" are slightly overstated.** PCA is an unsupervised linear representation learned from the family by SVD (SI Section S17), and β* is fit from the data. The honest and still-strong claim is "no iterative parameter optimization, no backpropagation, no pretraining corpus." I would soften "no learned representations" (lines 56, 86) accordingly; the contrast with deep models survives the correction.
5. **The β_gen = max(⌈2β*⌉, 5) rule is heuristic and its sensitivity is untested.** The factor of two (line 125, 148) is asserted as "just above the entropy inflection" but no sweep over the multiplier is shown. Given that the whole β* story is about principled temperature selection, a short sensitivity analysis over, say, 1.5×–3× would close the loop.

### Questions for Authors
1. How are gaps in retained seed sequences represented in the 20-channel one-hot encoding, and can a decoded sequence contain a gap?
2. How exactly is the "entropy inflection point" located numerically on the 50-point β grid (second-difference? spline?), and how stable is it to grid choice?
3. Which experimentally determined structures (PDB IDs) serve as the TM-align references for each family?
4. Does the distribution over *decoded* sequences inherit any of the guarantees that hold for the continuous Boltzmann distribution, or is decoding a heuristic projection whose effect on the sampled distribution is uncharacterized?

### Requested Experiments / Analyses
1. **Sequence-space diversity and novelty.** Report pairwise Hamming/identity diversity among the 150 generated sequences (and number of distinct clusters) per family, alongside the existing PCA cosine measures, so that "well-separated, not mode-collapsed" (line 60) is supported in sequence space.
2. **β_gen multiplier sensitivity.** Repeat generation at β_gen ∈ {1.5β*, 2β*, 3β*} for two or three representative families and report KL/novelty/identity, to show the 2× choice is not finely tuned.

### Minor Comments
1. The amino acid composition KL for "SA gen" is reported with different standard errors in Table S2 (e.g., RRM 0.060 ± 0.005, bootstrap) and Table S9 (RRM 0.060 ± 0.036, sub-block); the means match but the error bars do not. Reconcile or footnote the two SE conventions.
2. The abstract states "over 90% of matched SA-generated substitutions are experimentally tolerated"; this rests on two families (PDZ, SH3). Please qualify the breadth in the abstract (line 46).
3. Table 1 reports K_eff "effective number of sequences at 80% identity" but the clustering method/tool is not stated.
4. "Critical temperature" and "phase transition" are used for what is, at finite K, a smooth entropy crossover located on a coarse grid; consider "entropy inflection" consistently to avoid implying a true phase transition.
5. Define the significance-bracket convention (`***/**/*`) once and reference it across Figs. 3–4 and Tables S7.

### Recommendation
**Minor Revision.** The method is sound and the evaluation is strong; the issues above are clarifications, one or two overclaims, and two small analyses.

---

## Reviewer 2 — Hard

### Summary
This is a careful and well-instrumented study, and the closed-form-score derivation (lines 98–106) checks out: −β∇E = β(T(ξ) − ξ) with T the softmax attention map, and the ULA update with step α/β reduces exactly to Eq. 2. My concerns are with the statistical framing of the headline comparative claim, the resolution and independence of the β* regression, the diversity diagnostics, and a few reporting inconsistencies. Several of these are load-bearing for the paper's central message and should be addressed before publication.

### Strengths
1. **Correct, transparent mathematics.** The score is exact, the ULA derivation is right, and the regularity analysis (SI Section S19) honestly bounds the Lipschitz constant and dissipativity.
2. **Sampling diagnostics are above the norm.** Integrated autocorrelation times, per-chain ESS, MALA acceptance rates, and a burn-in margin analysis (SI Table S10) are all reported. Few generative-model papers do this.
3. **Two strong ablations.** PCA variance-threshold sensitivity (Table S16) and initialization sensitivity (Table S11) are exactly the right robustness checks, and the random-sphere init result (max ΔKL = 0.0003) is convincing evidence against an initialization artifact.

### Weaknesses
1. **The "30–70% within-family range" band does not match the measured identities of the studied families, yet the central comparative claim rests on it.** The main text repeatedly uses a 30–70% band as the "natural within-family identity range" to argue that HMM and MSA-Transformer outputs "moved outside the family" (lines 60, 68, 70, 86; Figs. 1C, 3C). But SI Table S12 reports the mean within-family pairwise identity of these eight families as **0.217–0.368** (22–37%), and the Discussion itself states the families span "22–37% mean pairwise identity" (line 92). HMM emit produces nearest-neighbor identity of 0.21 for RRM (Table S2-baselines), essentially equal to RRM's own 0.23 mean pairwise identity — so calling that "below the family range" is not supported by these families' statistics. The 30–70% band appears to be a generic literature value, not a quantity derived from the data being analyzed. This weakens the paper's single most-repeated claim. **Fix:** compute, per family, the leave-one-out nearest-neighbor identity distribution among the *stored* sequences and use that empirical distribution (not a fixed 30–70% band) as the reference envelope; re-evaluate which baselines actually fall outside it.
2. **The empirical β* is grid-quantized, and the regression is fit to nested, non-independent data.** The Table 1 β* values (3.23, 3.85, 4.58, 5.45) are consecutive points of the 50-point log grid on [0.1, 500] (ratio ≈1.19 between successive entries), so the "empirical" β* has ~19% multiplicative resolution. The regression target is therefore a step function of the true critical temperature, which inflates apparent R². Moreover, 25 of the 33 fitted points are nested subsamples of a single family (WW), so they are not independent observations. The author does refit on the eight families alone (R² = 0.95, line 66; SI Section S18), which partially addresses independence, but the headline "R² = 0.97, n = 33" overstates the evidence. **Fix:** refine the β grid (or fit on a continuous inflection estimate), report grid-induced uncertainty on β*, and lead with the eight-family fit and its leave-one-family-out result rather than the n = 33 number.
3. **"Mean pairwise cosine distance ≈ 1.0" is not evidence against mode collapse.** In d ≈ 34–186 dimensions, independent unit vectors are near-orthogonal by concentration of measure regardless of whether the sampler has collapsed; a cosine distance near 1.0 (line 60) is therefore expected even for a degenerate sampler and cannot support "distinct, well-separated sequences rather than collapsing to a single mode." **Fix:** demonstrate coverage with a sequence-space metric — pairwise Hamming-distance distribution, cluster counts (e.g., at 80% identity), or per-chain vs. across-chain diversity — not high-dimensional cosine.
4. **ULA-vs-MALA equivalence is argued only from acceptance rate.** A MALA acceptance of 99.6–99.8% at step size α = 0.01 (line 60; Table S10) is largely a statement that the step is small, and high acceptance alone is a weak guarantee of low ULA discretization bias. **Fix:** show that the *downstream observables* (KL, novelty, identity, and ideally the energy histogram) agree between ULA and MALA on two or three families; that is the evidence that justifies using ULA in production.
5. **The key control for the structural claims was not run through the structure pipeline.** The consensus-with-noise control (Table S14) and the column-permuted control (Table S15) are evaluated only on KL, novelty, identity, and MI — not on ESMFold/AlphaFold2 TM-score or pLDDT. But the most striking structural result is that SA's TM-scores *exceed* those of natural sequences (line 76), and the natural reading of that is regression toward consensus, which predictors favor. Without running consensus-with-noise (and permuted) sequences through the same structure predictors, the "structured exploration, not consensus" interpretation is not established at the structural level. **Fix:** fold the consensus and permuted ensembles through ESMFold/AF2 and compare TM/pLDDT against SA.
6. **Reporting inconsistencies in the uncertainty estimates.** Beyond the KL SE mismatch between Tables S2 and S9 (same means, different SEs), the main results Table S2 lists zf-C2H2 SA-gen KL as 0.013 ± 0.028 — a standard error larger than the mean for a non-negative quantity, which signals a poorly behaved bootstrap or a heavy-tailed per-block distribution. Please reconcile the SE conventions across tables and check the zf-C2H2 KL bootstrap.
7. **Multiple-comparison control is described for structure tests but not for the ESM2/DMS analyses.** Methods (line 155) specify Benjamini–Hochberg FDR at q = 0.05 for the Wilcoxon structural comparisons, but Table S7 (ESM2, eight families × tests) and the DMS analysis report raw p-values. With n in the thousands, the DMS enrichment is hugely significant yet modest in magnitude (1.18× over a 77.5% baseline for PDZ); the more informative question is whether the null (random single mutations) is the right comparator, given that SA substitutions are biased toward conserved-tolerant positions. Please state the correction applied across all families and justify the DMS null.

### Questions for Authors
1. What numerical procedure locates the entropy inflection on the 50-point grid, and what is the implied uncertainty on β* from grid spacing?
2. Why do the SA-generation KL standard errors differ between Table S2 and Table S9 for identical means, and which is the intended estimator for the figures?
3. For the DMS enrichment, is the comparator "all possible single mutations" or "single mutations at the positions SA actually mutated"? The two give very different baselines.
4. The convergence guarantees you cite (Durmus–Moulines, Dalalyan) hold for βσ²_max < 2, but SI Section S19 states the operating regime has βσ²_max ≫ 2. On what basis are sample-quality guarantees claimed (Discussion, line 92)?

### Requested Experiments / Analyses
1. **Empirical within-family identity envelope.** Replace the fixed 30–70% band with the per-family leave-one-out nearest-neighbor identity distribution of stored sequences, and re-derive the baseline "inside/outside the family" conclusions.
2. **Consensus and permuted controls through the structure predictors.** Run consensus-with-noise and column-permuted ensembles through ESMFold/AF2; report TM-score and pLDDT against SA and natural sequences.
3. **ULA vs. MALA on observables.** Compare the full metric set (and energy histograms) between ULA and MALA for two or three families.
4. **β grid refinement and β_gen multiplier sweep.** Refine the grid (report β* uncertainty) and sweep the β_gen multiplier to show the headline metrics are not sensitive to the 2× choice.

### Minor Comments
1. State the plmDCA regularization, sweep count, and software used for the Potts baseline in the main text, not only the SI (currently SI Section S9).
2. Table S2 caption: clarify that "SE via 1,000-fold bootstrap" applies only to the KL column.
3. The integrated autocorrelation cutoff ρ(τ) < 0.05 (SI Section S10) is a design choice; a one-line justification would help.
4. Report the number of decoded sequences that were exact duplicates of each other (intra-ensemble), not only duplication against stored patterns (Table S12).

### Recommendation
**Major Revision.** The method is correct and the instrumentation is good, but the central comparative claim depends on an identity band that the families' own statistics contradict, the β* regression is over-credited, and the structural superiority over natural sequences is not yet distinguished from consensus regression. These are addressable with re-analysis of existing data plus two modest new runs.

---

## Reviewer 3 — Very Hard

### Summary
The paper claims a training-free method that matches or exceeds learned generators on small protein families. The engineering is competent and the writing is careful, but the evidence for the strongest claims is entirely computational and, in several places, points toward an unflattering alternative explanation: that SA regresses toward the family consensus, which the validators (structure predictors, ESM2, DMS) systematically reward. The paper does not exclude that explanation, omits the most relevant baseline, and contains an internal contradiction about convergence guarantees. I would not accept it in its current form.

### Strengths
1. **The targeted gap is real.** The small-family regime is genuinely underserved, and the Pfam census (Table S15) quantifies it honestly.
2. **Self-awareness.** The author repeatedly notes that all validators are trained on evolutionary data and that single-mutant tolerance does not establish multi-mutant function (lines 80, 88). I credit this candor even as I argue it is not enough.

### Weaknesses
1. **Every "superiority" signal is consistent with consensus regression, which the validators reward, and the paper does not exclude it.** SA achieves TM-scores *higher than the natural sequences themselves* in six of eight families (line 76), passes a biophysical AMP filter *more often than natural defensins* (51% vs. 42%, line 82), and gets ESM2 pseudo-perplexity *better than natural* in Kunitz (line 80). A method cannot be more "natural" than nature unless the metric rewards something other than naturalness — the parsimonious explanation is that SA produces consensus-proximal sequences, and consensus sequences sit near the mode that predictors, language models, and central-tendency filters all favor. The author offers "Boltzmann averaging over lineage-specific deviations" as an alternative (line 88), but the decisive control — running the consensus-with-noise ensemble (Table S14) through the *same* structure/ESM2/DMS pipeline — was not performed. Until it is, "structural plausibility" reads as "proximity to consensus." **Fix:** put consensus-with-noise and column-permuted ensembles through every validator and show SA separates from them.
2. **The covariation claim is contradicted by the family where covariation matters most.** The paper argues the Hopfield energy "implicitly captures higher-order correlations" (lines 70, 78). Yet for Pkinase (L = 262, the longest family, where long-range contacts dominate), SA's MI correlation is r = 0.05 — *worse* than the position-independent HMM (0.23) — while Potts achieves 0.95 on the identical alignment (Table S6). The "high data-to-position ratio limits MI estimation" defense (line 78) is undercut by Potts succeeding on the same data. The most that can be claimed is that SA captures *some* covariation in *short* families; the higher-order-correlation narrative is not supported where it is testable and hardest.
3. **Internal contradiction on convergence guarantees.** The Discussion states the analytic score "interfaces with the convergence theory of Langevin dynamics, enabling formal guarantees on sample quality" (line 92), but SI Section S19 states plainly that the operating regime (βσ²_max ≫ 2) lies *outside* the log-concave regime where those guarantees hold, in a non-convex landscape with multiple basins. The guarantee claim should be removed or explicitly restricted to a regime the method does not use.
4. **The most relevant baseline is missing.** The paper's framing — deep models "overfit and collapse" at small K (lines 52, 90) — is asserted, not demonstrated. A VAE on the MSA (DeepSequence, Riesselman et al. 2018, already in the bibliography) is the canonical alignment-based deep generative model and the natural thing to show collapsing at K ≈ 22–45. Its absence means the central motivating claim is untested in this paper. Including it (and showing it collapse) would convert an assertion into evidence. Separately, the paper should engage with Russ et al., *Science* 2020 (chorismate mutase designed from a DCA/Potts model and experimentally validated in vivo): it is the strongest precedent that alignment-based generation can yield *functional* proteins, it sets the bar for what "validation" means in this area, and it shows experimental validation of exactly this class of method is feasible.
5. **The novelty over the author's own prior preprint is not delineated.** Reference 1 (Alswaidan & Varner 2026) "introduced SA and validated it on synthetic data, MNIST, and a single protein family" (line 54). This paper must state precisely what is new: scaling to eight families, the β* law, the baseline panel, and the structural/DMS validation, versus what was already established. As written, a reader cannot tell how much of the method is novel here.
6. **The speed claim leans on crippled baselines.** "3–5 orders of magnitude faster than EvoDiff and the MSA Transformer" (lines 74, 46) compares CPU-only runs of methods designed for GPUs. The author hedges ("GPU inference would reduce this substantially"), but the headline number in the abstract and Results still frames a CPU-vs-CPU artifact as a property of the methods. Fixed-length alignment generation also forces variable-length generators into a regime they were not built for, making the comparison partly apples-to-oranges.
7. **DMS validation: single-mutant tolerance is the wrong unit for a 17–25-substitution sequence.** The PDZ/SH3 enrichment (Table S8) scores substitutions *independently*, but SA sequences carry 17–25 simultaneous substitutions (line 62). Epistasis means individually tolerated mutations need not be jointly tolerated. The author concedes this (line 80), but then the DMS result cannot be offered as evidence for the "genuine structural convergence" interpretation in the Discussion (line 88) — it speaks only to marginal position/residue choice.

### Questions for Authors
1. When you run consensus-with-noise sequences through ESMFold/AF2/ESM2/DMS, does SA separate from them on any structural or fitness metric? If not, how is "structural plausibility" more than "proximity to consensus"?
2. Across families, does predicted TM-score (or pLDDT) correlate with each generated sequence's identity to the family consensus? A positive correlation would directly support the consensus-regression explanation.
3. What does a VAE (DeepSequence) trained on these same seed alignments produce at K ≈ 22–45, and does it in fact collapse as claimed?
4. Precisely which results here are new relative to Alswaidan & Varner 2026?
5. Why are SA's TM-scores higher than those of the natural sequences the method is imitating, and why should a reviewer read that as a feature rather than as evidence of metric bias?

### Requested Experiments / Analyses
1. **Run the consensus and permuted controls through every validator** (ESMFold, AlphaFold2, ESM2, and where applicable DMS), not only the composition/MI metrics, and report whether SA is separable from consensus-plus-noise.
2. **Add a VAE/DeepSequence-on-MSA baseline** at the studied family sizes to test, rather than assert, the "deep models collapse at small K" premise.
3. **Consensus-distance vs. validator-score regression.** For all generated sequences, regress predicted TM-score/pLDDT/ESM2 score on identity-to-consensus; report the slope. This is the direct test of the predictor-bias hypothesis the Discussion raises but does not resolve.
4. **Broaden DMS coverage using a standardized benchmark.** The PDZ and SH3 assays used here are part of ProteinGym (Notin et al., NeurIPS 2023); extending the cross-reference to the additional families with DMS data in that benchmark would make the tolerance-enrichment claim less anecdotal and the null better defined.

### Minor Comments
1. The abstract's "over 90% ... experimentally tolerated" generalizes from two families; restrict the wording.
2. "No learned representations" (lines 56, 86) is inaccurate given PCA; see also Reviewer 1.
3. The prior-work and EvoDiff citations are preprints (arXiv/bioRxiv); for a PNAS submission, note peer-review status where it matters to the claims.
4. State software versions (ESMFold, ColabFold/AF2, EvoDiff, ESM2, plmDCA) and random seeds, beyond the GitHub link, for reproducibility.
5. Figure 5 selects the single highest-pLDDT generated sequence per family for display; label it clearly as a best-case visualization (the caption does, but the main-text pointer at line 76 should too).

### Recommendation
**Major Revision** (with the understanding that, absent the consensus-regression controls and the missing baseline, the strongest structural and covariation claims are not yet supported). If the requested controls show SA separating from consensus-plus-noise across validators, my assessment would rise substantially.

---

## Summary of Actionable Items (Consolidated, Deduplicated, Prioritized)

### Tier 1 — Required to support the central claims
1. **Replace the fixed 30–70% "within-family" band with an empirical envelope.** Compute the per-family leave-one-out nearest-neighbor identity distribution of stored sequences and re-derive which baselines fall "outside the family." The current band contradicts the families' own 22–37% mean pairwise identity (R2-W1; R1 also). *Re-analysis of existing data.*
2. **Run the consensus-with-noise and column-permuted controls through the full structure/ESM2/DMS pipeline,** not only composition/MI. This is the decisive test that SA is more than consensus-plus-noise and underpins the TM-score-above-natural result (R2-W5, R3-W1). *New runs on existing ensembles.*
3. **Add a VAE/DeepSequence-on-MSA baseline** at K ≈ 22–45 to test, not assert, that deep models collapse in the small-family regime (R3-W4). *New baseline.*
4. **Regress validator scores on identity-to-consensus** to directly probe the predictor-bias / consensus-regression hypothesis (R3-W2, R3-Q2). *Re-analysis.*
5. **Remove or restrict the "formal convergence guarantees" claim** (Discussion line 92), which contradicts SI Section S19's βσ²_max ≫ 2 statement (R3-W3, R2-Q4). *Text.*

### Tier 2 — Strengthens rigor and fixes overclaims
6. **β* regression:** refine the β grid, report grid-induced uncertainty (current values are consecutive log-grid points, ~19% resolution), and lead with the eight-family fit rather than the WW-dominated n = 33 (R2-W2).
7. **Diversity in sequence space:** report pairwise Hamming/identity diversity and cluster counts; drop high-dimensional cosine distance as evidence against mode collapse (R2-W3, R1-Exp1).
8. **ULA vs. MALA on observables:** show metric and energy-histogram agreement, not just acceptance rate (R2-W4).
9. **β_gen multiplier sensitivity:** sweep 1.5×–3× on representative families (R1-W5, R2-Exp4).
10. **Pkinase covariation:** address directly that SA's MI (r = 0.05) is below HMM while Potts reaches 0.95 on the same data; temper the higher-order-correlation narrative (R3-W2).
11. **Delineate novelty** relative to Alswaidan & Varner 2026 explicitly (R3-W5).
12. **DMS:** state the multiple-comparison correction across families, justify the null (positions actually mutated vs. all positions), and consider extending to additional ProteinGym families; do not use single-mutant tolerance to argue multi-mutant structural convergence (R2-W7, R3-W7, R3-Exp4).

### Tier 3 — Clarity, reproducibility, and wording
13. **Specify gap handling** in the 20-channel one-hot encoding and decoding (R1-W1, R1-Q1).
14. **List PDB reference structures** (IDs, chains, residue ranges) used for TM-align (R1-W2, R1-Q3).
15. **Reconcile the KL standard errors** between Tables S2 and S9, and check the zf-C2H2 KL bootstrap (SE > mean) (R1-Minor1, R2-W6).
16. **Soften "no learned representations"** to "no iterative parameter optimization / no pretraining" (R1-W4, R3-Minor2).
17. **Qualify "over 90% tolerated"** in the abstract as based on two families (R1-Minor2, R3-Minor1).
18. **Reframe the speed comparison:** present CPU-only baseline timings without implying a GPU-vs-CPU artifact is a property of the methods; note the fixed-length constraint on variable-length baselines (R3-W6).
19. **Define the entropy-inflection numerical procedure** and the K_eff clustering method; standardize "entropy inflection" over "phase transition / critical temperature" (R1-Q2, R1-Minor3, R1-Minor4, R2-Q1).
20. **Report software versions and seeds**, and note peer-review status of preprint citations (R3-Minor3, R3-Minor4).

### Consensus across reviewers
All three reviewers regard the core method as sound and the empirical program as unusually thorough. The shared blocking concern is that the paper's strongest claims — structural plausibility, superiority to baselines, and capture of higher-order correlation — are not yet separated from the simpler hypothesis that SA regresses toward the family consensus, which every computational validator rewards. Tier 1 items (especially #1, #2, and #4) are the path from "Major Revision" to acceptance.

**Overall:** Reviewer 1 — Minor Revision; Reviewer 2 — Major Revision; Reviewer 3 — Major Revision.
