# arXiv Pre-Submission Review

## Critical Issues

**1. MI table inconsistency (SI Appendix)** ✅ Resolved

`tab:mi` and `tab:permuted-alignment` both claimed to report Pearson $r$ between stored and SA-generated pairwise MI vectors, but the values differed dramatically (e.g., RRM: 0.53 vs 0.87; Pkinase: 0.05 vs 0.71).

**Root cause:** `run_permuted_alignment.jl` used raw counts (no pseudocount) while `run_covariation_analysis.jl` used pseudocount=1.0. Fixed by adding pseudocount=1.0 to `compute_mi_vector` in `run_permuted_alignment.jl` and rerunning the analysis. Tables are now consistent.

---

**2. Incorrect count in permuted-alignment text** ✅ Resolved

`appendix.tex` stated MI differences "consistently favor the real alignment in seven of eight families", but the old (buggy) table showed only 4 families where SA (real) > SA (permuted), with 4 ties.

**Resolution:** After fixing the pseudocount bug (Issue #1), recomputed data shows 7 of 8 families with real > permuted (only zf-C2H2 is marginal at 0.92 vs 0.93). The "seven of eight" claim is now correct. The magnitude range in the accompanying text was also corrected from "0.00–0.03" to "0.01–0.08".

---

**3. Misleading KL vs. bootstrap comparison** ✅ Resolved

`results.tex:7` said SA KL values are "substantially lower than the bootstrap baseline (KL range: $0.007$--$0.133$)." But in SH3 and zf-C2H2, bootstrap achieves *lower* KL than SA gen (SH3: SA 0.036 vs. BS 0.030; zf-C2H2: SA 0.013 vs. BS 0.007). The claim was misleading.

**Resolution:** Replaced the "substantially lower" claim and the incorrect causal explanation with an accurate description: SA outperforms bootstrap in 6 of 8 families; the two exceptions are noted explicitly, with the explanation that large stored sets reduce bootstrap's finite-sample resampling variance.

---

## Moderate Issues

**4. MALA acceptance rate range mismatch** ✅ Resolved

`results.tex` said "range: 99.6%--99.9%", but `tab:sampling-diag` shows the maximum MALA AR is 0.998 = 99.8%. Corrected to "99.6%--99.8%".

---

## Minor Issues

**5. Abstract KL claim** ✅ Resolved

Abstract said "KL $< 0.06$" but the RRM SA gen value is exactly 0.060. Corrected to "KL $\leq 0.06$".

---

**6. Abstract sequence identity range** ✅ Resolved

Abstract said "SA maintains $52$--$66\%$ identity" but Pkinase SA gen SeqID is 0.515 = 51.5%, below the stated lower bound. Corrected to "$51$--$66\%$".

---

**7. Unreferenced figure** ✅ Resolved

`paper/figs/mi_covariation_barplot.pdf` was unreferenced. The figure visually shows the HMM < SA < Potts ≤ EvoDiff MI ordering across all 8 families simultaneously, directly supporting the paper's narrative. Added to the SI Appendix covariation section as `fig:mi-barplot` with caption and a reference in the body text.

---

**8. Style inconsistency** ✅ Resolved

Using NeurIPS preprint format for arXiv. The Significance Statement will be removed when reformatting for journal submission.

---

**9. Multiple `\label`s on one section** ✅ Resolved

`appendix.tex` had 4 `\label{}` commands on a single `\section{}`. All four cross-references from `theory.tex` pointed to the same section heading, even though they referred to distinct sub-results (concentration-of-measure, τ* Gaussian calculation, stored-pattern analysis, OLS regression). Fixed by retaining `\label{app:beta-star-analysis}` on the `\section{}` heading and inserting `\phantomsection\label{...}` anchors at the appropriate locations in the body text for the other three labels.

---

**10. Possible citation key typo** ✅ Resolved (no change needed)

`theory.tex` cites `ackeleyLearningAlgorithmBoltzmann1985`. The `.bib` file uses the same key, and the author field is correctly spelled `Ackley, David H.`. The extra 'e' is only in the BibTeX key name (cosmetic), not in the rendered bibliography. Citation compiles and displays correctly; no change required.
