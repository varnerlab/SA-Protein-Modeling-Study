# Antimicrobial Peptide (AMP) Experiment Plan

**Goal:** Apply SA to a therapeutically relevant, naturally scarce protein family — antimicrobial peptides — and validate generated sequences using biophysical property filters that don't require wet-lab work.

**Why this matters:** AMPs are the ideal use case for training-free generation: small families (10–100 members), short sequences (12–50 residues), and well-characterized physicochemical requirements. This experiment demonstrates practical utility beyond benchmarking.

---

## Data Sources

### Option A: Curated Pfam AMP families (preferred — consistent with rest of paper)
- **Defensins** (PF00323, vertebrate defensin): ~60 seed sequences, L≈35
- **Cecropin** (PF00272): ~20 seed sequences, L≈37
- **Magainin-like** — may not have a clean Pfam entry; check InterPro

### Option B: AMP databases (broader, but requires curation)
- **APD3** (Antimicrobial Peptide Database): https://aps.unmc.edu/
- **DRAMP** (Data Repository of Antimicrobial Peptides): http://dramp.cpu-bioinfor.org/
- Filter for: natural peptides, known activity, length 12–50, non-redundant at 80% identity
- Group by structural class: α-helical, β-sheet, extended, cyclic

### Recommendation
Use **Option A** (Pfam families) for consistency with the rest of the paper. Supplement with 1–2 curated AMP classes from APD3 if Pfam families are too small.

## Method

### Step 1: Data preparation
- Download Pfam seed alignments (same pipeline as pfam-families experiment)
- Clean alignment (remove >50% gap columns, >30% gap sequences)
- Build memory matrix: one-hot → PCA (95% variance) → unit-norm
- Find β* via entropy inflection

### Step 2: Generate sequences with SA
- Same parameters as main experiment: α=0.01, 30 chains × 5000 iterations, 5 samples/chain
- Two regimes: generation (β ≈ 2β*) and retrieval (β ≈ 20β*)
- Same baselines: bootstrap, Gaussian perturbation, convex combination

### Step 3: Biophysical property validation
AMPs have well-known physicochemical requirements. Compute for each generated sequence:

| Property | How to compute | Expected range for AMPs |
|----------|---------------|----------------------|
| **Net charge** at pH 7 | Sum of Lys(+1), Arg(+1), His(+0.5), Asp(−1), Glu(−1) | +2 to +9 (cationic) |
| **Hydrophobic ratio** | Fraction of residues in {A, V, I, L, M, F, W, P} | 0.40–0.60 |
| **Mean hydrophobicity** | Average Kyte-Doolittle score per residue | Slightly positive |
| **Amphipathicity** | Hydrophobic moment (Eisenberg) on α-helical projection | > 0.5 for helical AMPs |
| **Sequence length** | Count non-gap residues | 12–50 |
| **Cysteine count** | Count of Cys residues | 0 for helical, 4–6 for defensins |

These are all computable from sequence alone — no external tools needed.

### Step 4: Comparison to known AMPs
- Compare the biophysical property distributions of generated sequences vs. stored (real) AMPs
- 2D scatter: charge vs. hydrophobic ratio, colored by method
- KL divergence on each biophysical property distribution
- Fraction of generated sequences that pass a "plausible AMP" filter:
  - Net charge ≥ +2
  - Hydrophobic ratio ∈ [0.30, 0.70]
  - Length ∈ [12, 50]
  - Valid residue fraction = 1.0

### Step 5 (optional): AMP classifier prediction
- Use a pre-trained AMP classifier (e.g., AMPlify, iAMPpred) to score generated sequences
- This provides an independent "would this be predicted as antimicrobial?" check
- Requires Python dependency (AMPlify)

## Implementation Details

### Script structure
```
antimicrobial-peptides/
├── PLAN.md                          # this file
├── run_amp_experiment.jl            # main experiment script
├── data/
│   ├── <PF_ID>/                     # per-family data (same structure as pfam-families)
│   │   ├── stored_sequences.fasta
│   │   ├── sa_generation_sequences.fasta
│   │   ├── sa_retrieval_sequences.fasta
│   │   └── biophysical_properties.csv
│   └── amp_results.csv              # aggregated results
└── figs/
    ├── amp_charge_hydrophobicity.pdf     # 2D scatter
    ├── amp_biophysical_distributions.pdf  # violin plots
    └── amp_composite.pdf                  # multi-panel for paper
```

### Julia script outline (`run_amp_experiment.jl`)
1. Define AMP families (Pfam IDs)
2. For each family:
   a. Download + clean alignment (reuse `download_pfam_seed`, `clean_alignment`)
   b. Build memory matrix (reuse `build_memory_matrix`)
   c. Find β* (reuse `find_entropy_inflection`)
   d. Run SA generation + retrieval (reuse `run_sa_chains` pattern from pfam-families)
   e. Run baselines (bootstrap, GP, CC)
   f. Decode all sequences
   g. Compute biophysical properties for all sequences (stored + generated)
3. Aggregate into DataFrame with columns: Family, Method, Sequence, Charge, HydroRatio, MeanHydro, Length, ...
4. Compute pass rates and property KLs
5. Generate figures

### New functions needed (add to Protein.jl or local)
```julia
# Kyte-Doolittle hydrophobicity scale
const KD_HYDRO = Dict('A'=>1.8, 'R'=>-4.5, 'N'=>-3.5, 'D'=>-3.5, 'C'=>2.5,
                       'Q'=>-3.5, 'E'=>-3.5, 'G'=>-0.4, 'H'=>-3.2, 'I'=>4.5,
                       'L'=>3.8, 'K'=>-3.9, 'M'=>1.9, 'F'=>2.8, 'P'=>-1.6,
                       'S'=>-0.8, 'T'=>-0.7, 'W'=>-0.9, 'Y'=>-1.3, 'V'=>4.2)

function net_charge(seq; pH=7.0)  # simplified: Lys/Arg +1, Asp/Glu -1, His +0.5
function hydrophobic_ratio(seq)   # fraction in {A,V,I,L,M,F,W,P}
function mean_hydrophobicity(seq) # mean KD score
```

## Metrics to Report

| Metric | What it shows |
|--------|--------------|
| Charge distribution | Do generated AMPs maintain cationic character? |
| Hydrophobic ratio distribution | Is the hydrophobic content in the functional range? |
| "Plausible AMP" pass rate | Fraction meeting all biophysical criteria |
| Property KL divergence | How well SA preserves each property distribution |
| Standard metrics (KL_AA, Novelty, SeqID) | Consistency with main experiments |

## Paper Integration

- **Results section:** New paragraph: "Applied case study: antimicrobial peptides"
- **Experiments section:** Brief note on AMP family selection and biophysical filters
- **Figure:** Composite showing (A) charge vs. hydrophobicity scatter, (B) pass rate bar chart, (C) property distributions
- **Story:** "SA generates plausible AMP candidates from as few as 20–60 known sequences, preserving the biophysical properties required for antimicrobial activity"

## Risks & Mitigations

1. **Pfam AMP families too small (K < 20):** Supplement with APD3-curated sets, or frame this positively — "SA works even with 20 sequences"
2. **Biophysical filters too loose:** Tighten to match literature values for specific AMP classes (e.g., cathelicidins vs. defensins)
3. **Generated sequences may not be realistic AMPs:** This is expected for a training-free method. Frame as "candidate generation" not "guaranteed activity." The point is that SA preserves biophysical constraints better than baselines.
4. **Alignment quality for short peptides:** Short peptides may have ambiguous alignments. Use seed alignments from Pfam (curated) rather than automated alignment.
