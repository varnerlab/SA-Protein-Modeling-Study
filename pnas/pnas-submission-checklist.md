# PNAS Submission Checklist

## Formatting & Compliance
- [x] Word count: main text ~5,750 words (texcount). PNAS format-neutral initial submission does not enforce a strict count; the earlier "~4,000 word" target was a published-length approximation. Confirm two-column production length only if requested at revision.
- [x] Reference count: 32 cited in the main manuscript (within PNAS's ~50 guideline)
- [x] Abstract 240 words (≤ 250)
- [x] Significance statement 116 words (50–120)
- [x] Confirm 6 main-text figures (PNAS allows up to 6)
- [x] Data availability statement present (Zenodo DOI to be minted; placeholder `[to be assigned on deposit]` in main.tex)
- [x] Code availability statement present
- [x] ORCID for corresponding author (0000-0002-2558-7026 in main.tex)
- [x] Cover letter (drafted; includes editor suggestions and six suggested reviewers; arXiv preprint wording clarified)

## Content Quality
- [ ] Full proofread of main text (from rendered PDF)
- [ ] Full proofread of SI appendix (from rendered PDF)
- [x] Verify all SI cross-references match (re-verified after the reviewer-response edits; both PDFs build with 0 undefined references and no `??` in the output)
- [ ] Check all equations render correctly in PDF
- [ ] Verify figure quality/readability at print size
- [x] Confirm all citations/refs resolve: 0 `??` in both PDFs after the reviewer-response rebuild

## Reviewer-response (final-submission-issues-list.md)
- [x] All in-manuscript text revisions applied (see issues list for per-item status)
- [x] Position-matched + conservation-matched DMS null computed and written into Results + SI (Table S-dms-null)
- [x] Software versions + random seeds documented in SI (new "Software Versions, Random Seeds, and Reproducibility" section)
- [ ] Mint Zenodo archive and replace the DOI placeholder in the availability statement (user)
- [ ] Confirm all promised source/generated data are public on GitHub (user)
- [ ] Fill exact Python dependency versions (PyTorch, fair-esm, EvoDiff, ColabFold, NumPy/SciPy, TM-align) in the archived environment spec (user)
- [ ] Verify each arXiv identifier resolves: 2603.14717 (this preprint), 2603.06875 (prior SA), 2603.20115 (conditioning) (user)
- [ ] Verify portal metadata: title, ORCID, contributions, competing-interest, funding, acknowledgments, data-availability (user)
- [ ] Rebuild cover-letter.pdf after the arXiv wording change
