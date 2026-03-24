# TODO: Confirm Pfam seed alignment size statistics

## Claim to verify
The main text (Introduction, line 43) currently states:
> "the median Pfam seed alignment contains fewer than 50 sequences"

This cites Mistry et al. 2021, but that paper does not report seed alignment size distributions. We only have seed alignments for the 8 families used in the study (K = 37, 44, 45, 55, 68, 99, 151, 420), which is insufficient to verify a claim about all ~19,000+ Pfam families.

## What to do
1. Download the full Pfam seed alignment file from the EBI FTP: `https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.seed.gz` (~3GB compressed)
2. Parse Stockholm headers (`#=GF SQ` lines) to extract the number of seed sequences per family
3. Compute: median, fraction with <50, fraction with <100, fraction with <20
4. Update the Introduction claim with the actual numbers

## Current text to update
```latex
This creates a blind spot for small families: the median Pfam seed alignment contains fewer than 50 sequences~\cite{mistryPfamProteinFamilies2021}, and thousands of families have fewer than 20, well below the minimum viable training set for any deep generative architecture.
```

## Also check
- The abstract says "the median Pfam family has fewer than 100 members" — update to match whatever the seed alignment data shows
- The significance statement says similar — check consistency across all three locations
