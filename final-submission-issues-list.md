# Final Submission Issues

> Status added 2026-06-14 after a reviewer-response editing pass. `[x]` = handled in the
> manuscript/repo; items needing your action are marked **USER ACTION**.

## Critical Scientific Issues

- [x] **Qualify the Pfam family-size premise.** Intro now identifies the seed alignment as "the curated set of representative sequences that alignment-based methods take as input"; abstract and significance use "seed alignment"; SI Section S20 already explains seed-vs-full and reports the 841 full-alignment median. `pnas/main.tex` intro + `pnas/si-appendix.tex` S20.

- [x] **Correct the identity-envelope claims.** Abstract's generic `30--70%` band removed (now "each family's nearest-neighbor identity range"). Results state the envelope holds for seven of eight families and call out Pkinase as the exception (generated 0.515 vs natural max 0.454). Figure captions retain the 30–70% band but label it "approximate cross-family … envelope" and point to per-family SI values (Table nn-band).

- [x] **Structure control — SOFTENED (your choice).** Discussion now states explicitly that we do not read the elevated TM-scores as evidence SA folds better than natural sequences, and that the decisive control (folding the consensus-with-noise ensemble through the same predictors) was not performed here and remains future work. (CN sequences were never folded — only SA/stored/bootstrap were.)

- [x] **Support or narrow the ULA discretization-bias claim — NARROWED.** "Negligible discretization bias" removed in both main text and SI; replaced with the narrower statement that <0.4% of moves are rejected, which bounds rather than quantifies stationary-distribution bias.

- [x] **Qualify and strengthen the DMS result — DONE + COMPUTED.** Ran a position-matched null (and a conservation-matched cross-check). After controlling for SA's preference for tolerant positions, enrichment drops but stays significant: PDZ 91.7% vs 84.3% (1.09×, p<10⁻⁷⁷); SH3 90.4% vs 74.9% (1.21×, p<10⁻¹²¹). Abstract qualified to "two domains"; Results report the position-matched null; new SI table added (`tab:dms-null`). Script: `code/dms-validation/analyze_dms_position_matched_null.py`.

- [x] **Remove "no learned representations."** Replaced at both sites with "no iterative parameter optimization, backpropagation, or external pretraining"; significance now says "without model training."

## Important Claim and Framing Revisions

- [x] **Delineate novelty relative to the prior SA paper.** Added a compact intro paragraph listing the eight-family evaluation, scarce-data scaling, beta prediction, four-baseline panel, covariation, structure (ESMFold/AF2), ESM2, and DMS as new here.

- [x] **Temper the baseline-superiority language.** Results now state the timing figures reflect "the CPU-only, fixed-length evaluation protocol used here rather than a general method-level speed ranking," and note alanine padding/truncation penalize variable-length generators.

- [x] **Keep the Pkinase covariation limitation prominent.** Results already restrict covariation to seven of eight families with Pkinase (r=0.05) as the named exception; added an explicit covariation limitation to the Discussion limitations paragraph.

- [x] **Avoid equating nearest-neighbor identity with function.** "within the family's functional envelope" → "within the observed sequence-similarity envelope."

- [x] **Lead with the conservative beta regression evidence.** The eight-family fit (R²=0.95) + leave-one-family-out is now the primary evidence; the pooled 33-point fit is explicitly demoted to supporting because 25 points are nested WW subsamples.

## Submission and Reproducibility

- [x] **Manuscript length.** Light trim applied; checklist target corrected (the ~4,000-word figure was a published-length approximation, not enforced for a format-neutral initial submission). Main text ~5,750 words.

- [x] **Shorten the abstract and significance statement.** Abstract 240 words (≤250); significance 116 words (≤120).

- [x] **Remove the placeholder DOI.** `\doi{}` emptied; 0 occurrences of the XXXX placeholder in the compiled PDF.

- [ ] **Create a versioned software/data archive. — USER ACTION.** Availability statement rewritten to point to a Zenodo archive with DOI placeholder `[to be assigned on deposit]` and to drop the prior-preprint citation as the data record. Mint the archive and fill the DOI.

- [x] **Document software versions and seeds.** New SI section "Software Versions, Random Seeds, and Reproducibility": Julia 1.12.5 + package versions (Flux 0.16.9, NNlib 0.9.33, MultivariateStats 0.10.4, Distributions 0.25.123, StatsBase 0.34.10, DataFrames 1.8.1, CSV 0.10.16, Plots 1.41.6); model checkpoints (esm2_t33_650M_UR50D, esm_msa1b_t12_100M_UR50S, MSA_OA_DM_RANDSUB, ColabFold alphafold2_ptm); HMMER 3.4; seeds table. **USER ACTION:** fill exact Python dep versions (PyTorch, fair-esm, EvoDiff, ColabFold, NumPy/SciPy, TM-align) in the archived environment spec.

- [ ] **Confirm all promised data are public. — USER ACTION.** Verify seed alignments, generated/baseline FASTAs, structure metrics, DMS mappings, and per-table scripts are pushed to the public repo / archive.

- [x] **Preprint identifiers reconciled.** Confirmed the three IDs are different works: 2603.14717 (this manuscript's preprint), 2603.06875 (prior SA), 2603.20115 (conditioning). Cover letter reworded to say the preprint "is distinct from our earlier methods paper … (arXiv 2603.06875)." **USER ACTION:** confirm each ID resolves on arXiv.

## Final Package Checks

- [x] Rebuilt `main.pdf` (16 pp), `si-appendix.pdf` (27 pp), and `cover-letter.pdf` (2 pp).
- [x] No undefined citations/references, `??` markers, or DOI placeholder in the current build.
- [ ] Proofread main + SI from the rendered PDFs (equations, labels, legends, tables, cross-refs). — USER ACTION (recommended).
- [ ] Verify every edited numerical claim against the final tables. — partially: edited claims are tied to existing SI tables; do a final pass.
- [ ] Confirm portal metadata (title, ORCID, contributions, competing-interest, funding, acknowledgments, availability). — USER ACTION.
- [ ] Verify reviewer/editor suggestions for affiliations, suitability, COI. — USER ACTION.

## Build Status (post-edit)

- Main PDF: 16 pages; builds clean (exit 0).
- SI PDF: 27 pages; builds clean (exit 0); +1 new section (reproducibility) and +1 new table (DMS nulls).
- No unresolved references, undefined citations, or `??` in either build.
