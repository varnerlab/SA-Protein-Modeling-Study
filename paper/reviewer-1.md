# PNAS Peer Review

**Manuscript:** Training-Free Generation of Protein Sequences from Small Family Alignments via Stochastic Attention

**Recommendation:** Revise and resubmit (major revision, but publishable with revisions)

**Overall assessment:** This is a creative and well-executed paper that connects modern Hopfield networks to training-free protein sequence generation. The core idea is elegant, the experimental coverage across eight families is thorough, and the structure validation is convincing. However, several significant weaknesses need to be addressed before the paper meets the PNAS bar.

---

## Strengths

1. **Genuinely novel framing.** Treating the attention mechanism as a Boltzmann sampler and deriving an exact, closed-form score function is a clean theoretical contribution. The connection between Hopfield energy, attention, and Langevin dynamics is not just repackaging.

2. **Practical relevance to the long tail of protein space.** The paper correctly identifies a real gap: most families are too small for deep generative models. The method's O(dK) cost and laptop-level compute requirements are compelling.

3. **Structure validation is thorough.** Cross-validation with both ESMFold and AlphaFold2, with proper statistical tests and effect sizes, is above the standard for this type of paper. The resolution of the SH3 anomaly across predictors is well handled.

4. **The beta-star prediction (R^2 = 0.97) is striking** and provides a genuinely useful automation result.

5. **The Kunitz sequence analysis (Fig. 6)** is convincing evidence that conservation patterns emerge naturally from the energy landscape rather than being imposed.

---

## Major Concerns

**1. No experimental validation of any generated sequence.**
This is the most significant weakness. The paper claims to generate "structurally plausible" and even "antimicrobial" sequences, but all validation is computational (ESMFold, AlphaFold2). For a PNAS paper, at least one family should have wet-lab validation: express a handful of generated sequences, confirm they fold (CD spectroscopy), and ideally show some functional activity. The defensin analysis is particularly exposed here, as you claim 51% pass an "antimicrobial plausibility filter" but this filter is just charge + hydrophobicity + cysteine count. Any reviewer in the AMP field will note that these are necessary but far from sufficient conditions for antimicrobial activity. The limitations section acknowledges this ("wet-lab characterization ... is a natural next step") but for PNAS this may be a requirement, not a suggestion. **Mitigation if wet-lab is infeasible:** tone down the functional claims substantially, especially for defensins. Frame structure predictions as "consistent with foldability" rather than implying the sequences would function.

**2. The baselines comparison has an apples-to-oranges element.**
Profile HMMs, EvoDiff, and the MSA Transformer are designed for different tasks than SA. Profile HMMs are alignment models, not generative models in the modern sense. EvoDiff is pretrained on millions of sequences and is designed for open-ended generation, not family-constrained generation. You penalize them for producing sequences "too far from the family," but a user of EvoDiff might *want* that diversity. The framing should acknowledge that "staying close to the family" is a design choice of SA (due to the energy landscape), not an inherent advantage. A fairer comparison would also include:
- **Potts model / DCA-based generation** (Trinquier et al., which you cite but do not benchmark against). This is the most direct competitor: it also uses only the MSA, captures pairwise couplings, and generates without pretraining.
- **CARP or other VAE-based methods** trained on the family alignment alone, to test whether a lightweight learned model can match SA at small K.

**3. The method cannot generate insertions, deletions, or variable-length sequences.**
This is mentioned briefly in the limitations but deserves more prominence. Most real protein engineering tasks involve indels. The fixed-length constraint means SA cannot, for example, generate loop variants, which is a core application of protein design. This significantly narrows the practical applicability claim.

**4. The "higher TM-scores than natural sequences" claim needs more careful interpretation.**
You attribute this to SA generating sequences closer to a "consensus fold." But an alternative explanation is simpler: ESMFold/AF2 may systematically predict more canonical structures for consensus-like sequences because their training data is biased toward well-characterized (and thus consensus-like) proteins. The structure prediction confidence could be conflated with structural accuracy. The paper should discuss this confound explicitly, and if possible, test it by predicting structures of the actual family consensus sequence and comparing its TM-score.

**5. The PCA representation discards higher-order correlations by construction.**
PCA captures only pairwise (linear) correlations. The paper argues this is sufficient, but a DCA/Potts model would capture pairwise epistasis directly, and the claim that SA "faithfully recapitulates the position-specific conservation landscape" (Kunitz analysis) shows single-site conservation, not pairwise coupling preservation. Do the generated sequences preserve known covariation patterns (e.g., between contacting residues)? Without this analysis, the claim of "faithful recapitulation" is limited to marginal statistics.

---

## Minor Concerns

**6. The MALA validation is perfunctory.** 99.6% acceptance is reported, but this does not tell us about mixing. What is the effective sample size? What is the autocorrelation time? How many independent samples does each chain actually produce? You discard 2,000 burn-in and thin every 100, but the justification for these choices is not given.

**7. The eight families are all well-studied, well-structured domains.** Would SA work on intrinsically disordered regions, multi-domain proteins, or families with high structural heterogeneity? Acknowledging this bias in family selection would strengthen the paper.

**8. The defensin analysis in the Discussion feels out of place.** It introduces new quantitative results (cysteine counts, charge distributions, antimicrobial filter pass rates) that are not in the Results section. Either move this analysis to Results with a proper figure/table, or reduce it to a qualitative remark.

**9. Some numerical inconsistencies between sections.** The beta-star regression is reported as both R^2 = 0.969 (results) and R^2 = 0.962 (theory/appendix), and the intercept is 1.52 (results) vs 1.57 (theory) vs 1.565 (appendix). These likely reflect rounding or slightly different datasets, but they will confuse reviewers. Pick one set of numbers and use them consistently.

**10. The introduction is long and reads like a related-work survey.** For PNAS, consider tightening the first two paragraphs (the literature review of all generative models for proteins) and getting to the key insight faster.

**11. Writing style.** The results section is dense and reads like a long list of numbers. Consider whether some of the per-family statistics could be summarized more concisely with references to tables/figures, rather than enumerated inline.

**12. The scaling study uses only WW subsamples.** This means you are testing scaling within a single family's sequence space, not across families. The claim "stable down to K=20" may not generalize to families with very different sequence landscapes.

---

## Questions for the Authors

1. Have you tested whether the generated sequences preserve known pairwise covariation (e.g., mutual information or DCA scores between contacting residue pairs)?
2. What happens when the PCA variance threshold is changed from 95% to 90% or 99%? How sensitive are the results?
3. For the beta-star regression, the 25 WW replicates and 8 families are pooled. The WW points dominate (25/33). Is the R^2 still high using only the 8 families?
4. The multi-chain protocol initializes near stored patterns. What happens with random initialization on the sphere?

---

## Summary

The core contribution (training-free generation via exact Hopfield scores) is sound and the experimental evaluation is above average for a methods paper. The main vulnerability is the absence of any experimental validation, which PNAS reviewers in the protein design space will likely flag as a hard requirement. The baseline comparison should include DCA/Potts-based generation as the most apples-to-apples competitor. The numerical inconsistencies and the defensin analysis placement are easily fixed. With these revisions, this could be a strong PNAS paper.
