# Reviewer 2: Second Round Review

The authors have substantially improved the manuscript. The MI analysis, sampling diagnostics, PCA sensitivity study, initialization experiment, and 8-family regression are all welcome additions that address Reviewer 1's concerns comprehensively. The paper is now better substantiated. I have the following remaining concerns, organized by severity.

## Major

**1. The "diversity" metric is defined but never reported.** The experiments section (line 7) defines five metrics including "Diversity" (mean pairwise cosine distance among generated samples), but I cannot find diversity values anywhere in the results, tables, or SI. Either report it or remove it from the metrics list. This is a minor fix but looks careless.

**2. The AF2 cross-validation is incomplete.** The results state "additional families are in progress" for AlphaFold2 (results line 27). This is not acceptable for publication. Either complete the analysis for all eight families or explicitly limit the AF2 claim to the families tested. The phrase "in progress" belongs in a preprint, not a journal submission.

**3. The Pkinase family is a persistent outlier that undermines universality claims.** It appears as an outlier in MI analysis (r = 0.05), has the most extreme data-to-parameter ratio, the flattest PCA sensitivity curve (because d ~ K), and is the family where EvoDiff catastrophically fails. The paper acknowledges some of these individually but never confronts the pattern: is SA actually working well on Pkinase, or is the family simply too data-starved for any method? A candid paragraph acknowledging Pkinase as a stress test (37 sequences encoding 262 positions) rather than a success story would strengthen credibility.

## Minor

**4. The discussion is now quite long (~7 substantial paragraphs plus limitations, future directions, and summary).** The binder design paragraph (paragraph 7) is speculative and unsupported by any data in this paper. For a PNAS article, consider whether this paragraph belongs here or in a separate perspective/commentary. It reads like a grant aim rather than a discussion point arising from the current results.

**5. "Predicted to fold" hedging is inconsistent.** The intro (line 11) says "predicted to fold into the correct three-dimensional architecture," which is appropriately hedged. But the discussion opening (line 3) says "fold into the correct three-dimensional architecture" without the hedge. This should be consistent throughout.

**6. The convex combination baseline is a straw man.** It appears in every comparison and always performs terribly (KL > 0.3 everywhere). Its inclusion inflates the apparent advantage of SA. The more informative comparisons are against HMM, EvoDiff, and Potts. Consider whether this baseline adds value or just takes up space in the figures and tables.

**7. The scaling study is single-family.** This is acknowledged in the limitations, but it weakens the claim that SA is robust at small K. The WW domain at K=20 may not be representative of other families at K=20. Even a quick 2-3 family replication at K=20 (subsampling SH3, Kunitz, or zf-C2H2) would substantially strengthen this point.

**8. Missing wall-clock timing table.** The paper repeatedly claims SA runs "in seconds" vs. hours for EvoDiff, but I don't see a systematic timing comparison anywhere. A small table reporting wall-clock time per 150 sequences (SA, HMM, EvoDiff, MSAT, Potts) would be more convincing than scattered prose claims.

**9. The one-hot dimension is stated as 21L in the discussion (line 13) but the encoding section in the SI describes a 20-letter alphabet giving 20L.** Which is correct? If gap characters are included, say so explicitly. If not, fix the discussion.
