# Structure Validation Experiment Plan

**Goal:** Demonstrate that SA-generated sequences fold into structures consistent with their family, using ESMFold as an in silico oracle.

**Why this matters:** Reviewers will ask "do these sequences actually fold?" Sequence-level metrics (KL, novelty) don't answer this. Structure prediction provides functional plausibility without wet-lab work.

---

## Data

- **Input:** The 150 SA-generated sequences per family already produced by the pfam-families experiment (stored in `../pfam-families/data/<PF_ID>/sa_generation_sequences.fasta` and `sa_retrieval_sequences.fasta`)
- **Reference structures:** One representative PDB structure per family (download from PDB or use AlphaFold DB)
- **Families:** All 7 Pfam families (RRM, SH3, WW, Kunitz, zf-C2H2, PDZ, Pkinase)

## Method

### Step 1: Structure prediction via ESMFold API
- Use the ESMFold REST API (`https://api.esmatlas.com/foldSequence/v1/pdb/`) — no GPU, no Python dependency, no installation
- Call directly from Julia using `HTTP.jl`: POST the amino acid sequence, receive a PDB file
- Throttle to ~1 request/second to respect rate limits
- Predict structure for:
  - All 150 SA generation sequences per family
  - All 150 SA retrieval sequences per family
  - 150 bootstrap sequences (as baseline)
  - All stored sequences (positive control — these should fold well)
- Save PDB outputs to `data/<PF_ID>/structures/<method>/`
- Total: 7 families × ~500 sequences ≈ 3,500 requests (~1 hour at 1 req/s)
- **Fallback:** If the API is unavailable or rate-limited, run ESMFold locally on CPU. At ~30s/sequence and ~3,500 sequences this is ~29 hours — feasible as a weekend run.

### Step 2: Confidence metric — pLDDT
- Extract per-residue pLDDT from ESMFold PDB output (stored in the B-factor column)
- Compute mean pLDDT per sequence
- Report distribution across generated sequences vs. stored sequences
- Threshold: pLDDT > 70 is "confident", > 90 is "very high confidence"

### Step 3: Structural similarity — TM-score
- Compare each predicted structure to the family's reference structure using TM-align
- TM-score > 0.5 indicates same fold; > 0.7 indicates high structural similarity
- Install/use `TMalign` binary (standard structural bioinformatics tool)
- Report distribution of TM-scores for SA gen vs. SA ret vs. bootstrap vs. stored

### Step 4: Contact map comparison (optional, if time permits)
- Extract Cα contact maps (8Å cutoff) from predicted structures
- Compare to reference contact map using precision/recall
- This is a more fine-grained structural metric than TM-score

## Implementation Details

### ESMFold API from Julia
```julia
using HTTP

function esmfold_predict(sequence::String)
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    response = HTTP.post(url, [], sequence; retry=true, retries=3, readtimeout=120)
    return String(response.body)  # PDB-format string
end

function extract_mean_plddt(pdb_string::String)
    # pLDDT is stored in the B-factor column (columns 61-66) of ATOM records
    plddts = Float64[]
    for line in split(pdb_string, '\n')
        if startswith(line, "ATOM") && length(line) >= 66
            plddt = parse(Float64, strip(line[61:66]))
            push!(plddts, plddt)
        end
    end
    return mean(plddts)
end
```

### TMalign dependency
- Download TMalign binary from https://zhanggroup.org/TM-align/
- Place in `bin/` or install system-wide (`brew install tmalign` if available)
- Call from Julia via `run()` and parse stdout for TM-score

### Reference PDB structures
One representative structure per family (well-resolved, canonical fold):

| Family | Pfam ID | Candidate PDB | Notes |
|--------|---------|---------------|-------|
| RRM | PF00076 | 1FXL or 2CJK | Classic RRM fold |
| SH3 | PF00018 | 1SHG | Src SH3 |
| WW | PF00397 | 1PIN | Pin1 WW domain |
| Kunitz | PF00014 | 1BPI | BPTI |
| zf-C2H2 | PF00096 | 1ZAA | Zif268 |
| PDZ | PF00595 | 1BE9 | PSD-95 PDZ3 |
| Pkinase | PF00069 | 1ATP | PKA catalytic |

Verify these at experiment time; choose chain with best resolution.

### Script structure
```
structure-validation/
├── PLAN.md                          # this file
├── run_structure_validation.jl      # main experiment script (Julia only)
├── data/
│   └── <PF_ID>/
│       ├── structures/
│       │   ├── sa_generation/       # PDB files
│       │   ├── sa_retrieval/
│       │   ├── bootstrap/
│       │   └── stored/
│       ├── plddt_scores.csv         # per-sequence mean pLDDT
│       └── tm_scores.csv           # per-sequence TM-score to reference
├── results/
│   └── structure_validation_results.csv  # aggregated results
└── figs/
    ├── plddt_distribution.pdf
    ├── tmscore_distribution.pdf
    └── structure_validation_composite.pdf
```

### Julia script outline (`run_structure_validation.jl`)
```
1. using HTTP, CSV, DataFrames, Statistics, Plots
2. Include ../Include.jl for parse_fasta, etc.

3. For each family:
   a. Read generated FASTA files from ../pfam-families/data/<PF_ID>/
   b. Strip gap characters from decoded sequences
   c. For each sequence:
      - POST to ESMFold API → PDB string
      - Extract mean pLDDT from B-factor column
      - Save PDB file to data/<PF_ID>/structures/<method>/
      - Sleep 1 second (rate limiting)
   d. Download reference PDB structure from RCSB
   e. Run TMalign on each predicted PDB vs. reference → TM-scores
   f. Save per-family CSV of results

4. Aggregate all families into one DataFrame
5. Compute summary stats: mean ± SE of pLDDT and TM-score per method per family
6. Generate figures
7. Save composite figure for paper
```

### Checkpointing
Since the API calls take ~1 hour and could fail midway:
- Save each PDB file immediately after receiving it
- Before calling the API, check if the PDB file already exists (skip if so)
- Save the pLDDT CSV incrementally after each family completes
- This allows resuming from where we left off if the script is interrupted

## Metrics to Report

| Metric | What it shows | Reported as |
|--------|--------------|-------------|
| Mean pLDDT | Prediction confidence | mean ± SE across sequences |
| Fraction pLDDT > 70 | "Foldable" rate | percentage |
| TM-score to reference | Structural similarity to family | mean ± SE |
| Fraction TM > 0.5 | Same-fold rate | percentage |

## Paper Integration

- **Results section:** New paragraph + figure after the cross-family results
- **Experiments section:** Brief paragraph on ESMFold setup ("We used the ESMFold API to predict structures for generated sequences without additional training or GPU resources.")
- **Figure:** 2-panel composite (pLDDT distributions + TM-score distributions), grouped bar or violin by method × family

## Estimated Compute

- ESMFold API: ~1 request/second, ~3,500 sequences → ~1 hour
- TMalign: seconds per pair → ~10 minutes total
- Total: ~1.5 hours (no GPU needed)
- **Fallback (CPU):** ~29 hours — feasible as a weekend run

## Risks & Mitigations

1. **ESMFold API unavailable or deprecated:** Fall back to local CPU installation (weekend run). Check API availability before starting.
2. **API rate limiting:** Throttle to 1 req/s. If stricter limits encountered, reduce to 20 sequences per method per family (560 total, ~10 min).
3. **Short sequences (WW=31aa, zf-C2H2=23aa):** ESMFold may give low pLDDT even for real sequences. Mitigation: always compare to stored sequence pLDDT as positive control. The comparison is relative (SA gen vs. stored), not absolute.
4. **Gap positions in decoded sequences:** Strip all '-' characters before submitting to ESMFold.
5. **Reference structure selection:** Use well-resolved X-ray structures. If multiple chains, use the one matching the Pfam domain boundaries.
6. **Network interruption during API calls:** Checkpointing (see above) allows seamless resume.
