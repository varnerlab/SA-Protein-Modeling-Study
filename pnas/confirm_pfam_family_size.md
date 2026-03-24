# DONE: Pfam seed alignment size statistics verified

## Results (Pfam 36.0, 20,795 families)
- **Median seed alignment: 22 sequences** (IQR: 7--65)
- 47.5% of families have fewer than 20 seed sequences
- 69.2% have fewer than 50
- 84.3% have fewer than 100
- Mean: 59.5, Min: 1, Max: 3,769

## What was done
1. Downloaded `Pfam-A.seed.gz` (151 MB) and `Pfam-A.full.gz` (19.2 GB) from Pfam 36.0 release
2. Parsed all Stockholm entries to count unique sequences per family
3. Updated all four text locations (abstract, significance, introduction, cover letter)
4. Added SI Appendix Section S20 with figure and table documenting the census
5. Added InterPro citation (Paysan-Lafosse et al. 2023)

## Scripts
- `code/pfam-census/run_pfam_census.jl` -- download + parse + save CSV
- `code/pfam-census/plot_pfam_census.jl` -- generate distribution figure

## Updated text locations
- Abstract (main.tex line ~37): "the median Pfam seed alignment contains only 22 sequences"
- Significance (main.tex line ~23): "a median of only 22 sequences in their Pfam seed alignments"
- Introduction (main.tex line ~43): "the median Pfam seed alignment contains only 22 sequences, and nearly half of all families have fewer than 20"
- Cover letter (cover-letter.tex line ~17): matching language

## Note
The original claim ("fewer than 50 sequences") was correct but conservative. The actual median (22) makes a stronger case for the paper's motivation.
