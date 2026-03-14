# Reviewer 3

**Manuscript:** Training-Free Generation of Protein Sequences from Small Family Alignments via Stochastic Attention

**Recommendation:** Reject (substantial concerns about novelty, validation, and overclaiming)

**Overall assessment:** The paper proposes using Langevin dynamics on the modern Hopfield energy as a training-free protein sequence generator. The mathematical framework is correct, and the experimental coverage across eight families is thorough by computational standards. However, the contribution is narrower than the paper claims, the validation is entirely circular within computational tools trained on the same evolutionary data being modeled, and several core claims do not survive careful scrutiny. The paper would benefit from a substantially more honest framing of what the method can and cannot do.

---

## Major Concerns

**1. The methodological contribution is incremental.**
The core technical idea, Langevin sampling on the modern Hopfield energy with an analytic score function, was already published by the same author (Alswaidan & Varner 2025, ref [6]). That paper derived the score function, proved the convergence properties, and demonstrated stochastic attention on a single protein family (RRM). The present manuscript extends this to seven additional families and adds computational validation. While the protein-specific analyses are competent, the intellectual contribution of the present work reduces to "we ran the same method on more families and added structure prediction." The theoretical sections (Eqs. 1-6) are reproduced verbatim from the prior work. For a venue like PNAS, which requires significant conceptual advance, this is a concern. The $\beta^*$ prediction (Eq. 6) is the most novel element, but it is a purely empirical correlation fit on 33 data points (25 of which are WW replicates from a single family) with no held-out validation. The R^2=0.97 is in-sample and likely inflated by the autocorrelation among WW replicates. The 8-family-only R^2 drops to 0.95 on 8 points, which is not statistically distinguishable from many simple functions of d.

**2. All validation is computationally circular, and the paper does not adequately acknowledge this.**
The generated sequences come from an alignment of evolutionary homologs. They are validated by ESMFold and AlphaFold2 (trained on evolutionary homologs from PDB + UniRef), scored by ESM2 (trained on UniRef), and cross-referenced against DMS data from well-characterized families. Every validation tool has seen the same evolutionary signal that defines the input alignment. A method that generates consensus-like sequences from a family alignment will, almost by definition, score well on tools trained to recognize evolutionary consensus. This is not a bug in any individual validation, but the cumulative effect is that the paper provides zero evidence that SA-generated sequences have any property that a consensus sequence or a position-frequency-reweighted sample would not also have. The paper acknowledges this for structure prediction ("both predictors are neural models trained on evolutionary data") but then proceeds as if the acknowledgment resolves the concern. It does not. The appropriate control would be to generate sequences from a scrambled alignment (same amino acid frequencies, destroyed covariation) and show that SA produces better-scoring sequences than scrambled-alignment SA. Without such a control, the results are consistent with the null hypothesis that SA is simply regenerating the consensus with noise.

**3. The "higher TM-scores than natural sequences" result is almost certainly an artifact, and the paper's treatment of it is misleading.**
The paper presents two explanations (consensus fold vs. predictor bias) and then says "the consensus-fold explanation is favored." But these are not mutually exclusive, and the predictor-bias explanation is far more parsimonious. Both ESMFold and AlphaFold2 are trained on PDB structures, which are biased toward well-folding, stable, consensus-like proteins. A method that generates consensus-proximal sequences will get higher confidence and TM-scores from these predictors, regardless of actual structural quality. The concordance between two predictors does not resolve this, since both share the same training data bias. The paper should present this result with appropriate skepticism rather than as a positive finding. As written, a reader might conclude that SA-generated sequences fold better than natural ones, which is an extraordinary claim that requires experimental evidence.

**4. The comparison framework is constructed to favor SA, and the framing obscures this.**
The paper evaluates four baselines: profile HMMs, EvoDiff, MSA Transformer, and Potts (the latter on only 2 families). Each comparison is framed to highlight SA's advantage:

- Profile HMMs are penalized for producing sequences with low identity to stored patterns. But low identity is the point of HMMs: they model the family probability distribution, not individual members. A sequence with 30% identity that is within the family distribution is arguably more useful for protein engineering than one with 60% identity that is a consensus echo.
- EvoDiff is penalized for computational cost and for a single poor result on Pkinase (KL=0.20). But EvoDiff is designed for open-ended generation from large MSAs, not for small-family generation. The comparison is like benchmarking a bicycle against a truck on a narrow path and concluding bicycles are superior vehicles.
- The MSA Transformer is used for generation via iterative masked filling, which is not its intended use case. The paper acknowledges none of this.
- The Potts model, which is the only apples-to-apples comparison, achieves comparable results on both families tested. The paper acknowledges this but then argues SA is better because it avoids fitting parameters. This is true but cuts both ways: the Potts model fits parameters that capture real pairwise couplings, while SA captures them only "implicitly" through an uncontrolled mechanism.

The paper should have included ancestral sequence reconstruction (a training-free method that generates biologically meaningful variants) and simple consensus-with-noise generation as baselines. Without these, it is impossible to distinguish SA from a consensus generator with a diversity knob.

**5. The PCA representation is a fundamental bottleneck, and the MI analysis does not adequately address this.**
PCA on one-hot encodings captures linear correlations in binary indicator space. This is not the same as capturing pairwise residue couplings, which are nonlinear relationships between categorical variables. The MI analysis (Table S5) shows Pearson r=0.53-0.92 between stored and generated MI matrices, but this comparison is confounded: any method that preserves marginal amino acid frequencies will show positive MI correlation with the stored alignment, because MI is influenced by marginal frequencies. The appropriate control is to compare SA's MI preservation against a position-independent model that perfectly matches the marginal frequencies (i.e., the null model for MI). The comparison against profile HMMs is helpful but insufficient, since HMMs are known to be poor at preserving pairwise structure. The Pkinase failure (r=0.05) is particularly telling: in the family with the most extreme dimensionality reduction (5240-dimensional one-hot projected to d=34, a 154x compression), the method fails to preserve any pairwise structure. This is the regime the paper claims to address (small families with many positions), and it is exactly where PCA discards the most information.

**6. The defensin biophysical analysis proves less than it appears to.**
The antimicrobial plausibility filter (charge >= +2, hydrophobic ratio in [0.3, 0.7], >= 4 cysteines) is so permissive that 42% of natural defensins fail it. A filter with a 58% false negative rate on known functional sequences is not a meaningful assessment of biological plausibility. The fact that SA achieves 51% pass rate (compared to 42% for stored sequences) does not indicate that SA-generated sequences are more antimicrobial; it indicates that the filter favors consensus-like sequences, which natural defensins (with their lineage-specific deviations) sometimes violate. The entire analysis could be removed without loss.

---

## Minor Concerns

**7. The abstract and significance statement overclaim.**
"Competitive, training-free alternative to learned generative models" is not supported by the data. SA produces consensus-proximal sequences with moderate novelty from small alignments. It does not generate diverse functional variants, design proteins with specific binding properties, or produce sequences with validated biological activity. The significance statement claims SA "generates novel protein sequences that ... fold into the correct three-dimensional structure," which is unsupported in the absence of experimental structure determination.

**8. The scaling study is still effectively single-family.**
The K=20 replication on three additional families (Table S2) uses families that are already small (SH3 K=55, Kunitz K=99, zf-C2H2 K=151). Subsampling these to K=20 removes 64-87% of the data. The resulting metrics (KL < 0.03, novelty ~0.47) are nearly identical across all four families, which is suspicious: it suggests the method converges to a universal output profile at K=20 regardless of family identity, rather than capturing family-specific structure. Are the generated sequences at K=20 actually different across families, or has the method collapsed to a generic "small protein family" mode?

**9. The paper is a single-author study with no experimental collaborator.**
This is not inherently disqualifying, but it explains the absence of any experimental validation, which multiple reviewers have flagged as a critical gap. For a paper claiming relevance to protein engineering and therapeutic peptide design, the absence of even a single experimental test (expression, solubility, CD spectrum) is a significant limitation. The paper should either include experimental data or substantially narrow its claims to the computational domain.

**10. Statistical testing concerns.**
The paper reports Wilcoxon tests and Cohen's d across 8 families and multiple metrics without any correction for multiple comparisons. With ~40 statistical tests, some nominally significant results are expected by chance. The effect sizes for pLDDT comparisons are often small (d < 0.5), making the practical significance questionable even when statistical significance is achieved. The bootstrapped KL standard errors (1000 resamples of 150 sequences) may underestimate true uncertainty because the 150 sequences come from 30 chains with potentially correlated samples.

**11. The β* prediction is not validated out-of-sample.**
The regression β* = 1.57 + 0.28√d is fit on 33 points and evaluated on the same 33 points. The paper presents this as a predictive model, but no cross-validation or held-out test is reported. Leave-one-family-out cross-validation would be straightforward and much more convincing. As it stands, the R^2=0.97 may reflect overfitting to the WW-dominated training set.

**12. The method has no mechanism for functional design.**
The paper's Discussion includes a speculative paragraph about directed generation of binding proteins, but the method as presented has no way to incorporate functional constraints, binding targets, stability requirements, or any objective beyond "sample from the family distribution." This paragraph should be removed or clearly labeled as speculation, as it creates expectations the method cannot meet.

**13. Runtime comparisons are misleading.**
The paper compares SA runtime (seconds) against EvoDiff and MSA Transformer (hours) on CPU. But EvoDiff and the MSA Transformer are designed for GPU execution, where they run orders of magnitude faster. The CPU timing for EvoDiff (14 hours for Pkinase) is not a fair representation of that method's practical runtime. If GPU timings are unavailable, the paper should acknowledge that the CPU comparison inflates SA's speed advantage.

**14. The convex combination baseline is still present despite Reviewer 2 flagging it as a straw man.**
It appears in Table S1 and inflates SA's apparent advantage. KL > 2.0 for a convex combination is expected (averaging one-hot vectors produces a mixture that looks nothing like any individual sequence) and provides no information about SA's quality. This baseline should be removed or relegated to a brief footnote.

---

## Questions for the Authors

1. What happens if you generate sequences from the family consensus plus calibrated Gaussian noise in one-hot space (i.e., a simple consensus-with-noise baseline)? How do the metrics compare to SA?
2. Can you provide leave-one-family-out cross-validation for the β* prediction model?
3. For the K=20 scaling experiment, are the generated sequences actually family-specific, or do all four families produce similar outputs? What is the inter-family cosine similarity of the generated ensembles?
4. What fraction of SA-generated sequences are identical to or within 1-2 substitutions of a stored sequence? The novelty metric (PCA-space cosine) may mask near-duplicates in sequence space.
5. Have you tested SA on a family where the seed alignment is known to contain errors or misaligned sequences? Robustness to alignment quality is critical for real-world use.
6. The Discussion mentions extending SA to variable-length generation via "length-agnostic embeddings." Have you attempted this? If the answer is no, this limitation effectively restricts SA to the subset of protein design problems where indels are irrelevant, which is a small minority.
7. EvoDiff and the MSA Transformer are penalized for generating sequences with low identity to stored patterns (11-59%). But these methods are sampling from a broader distribution that may include functional variants not represented in the seed alignment. Has any analysis been done to determine whether the "drifted" baseline sequences are actually non-functional, or is low identity being used as a proxy for poor quality without justification?

---

## Summary

The mathematical framework is sound, and the experimental evaluation is thorough within its computational bounds. However, the paper suffers from three interconnected problems: (1) the methodological contribution over the prior publication is incremental, (2) all validation is circular within the ecosystem of evolutionary-data-trained tools, and (3) the framing systematically overclaims what the results demonstrate. The method generates consensus-proximal sequences from small alignments, which is useful but modest. Presenting this as "competitive with learned generative models" or relevant to "therapeutic peptide design" is not supported by the evidence. A revised version that honestly scopes the contribution, includes a consensus-with-noise baseline, performs leave-one-family-out validation for β*, and either adds experimental data or explicitly limits claims to the computational domain could be a solid contribution to a computational biology venue. In its current form, the gap between claims and evidence is too large for PNAS.
