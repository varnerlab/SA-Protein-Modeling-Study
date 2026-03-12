# Learned Baselines Experiment Plan

**Goal:** Compare SA against protein-specific learned generative models on the same families, demonstrating that training-free SA is competitive in the scarce-data regime (K < 200) where learned models struggle.

**Why this matters:** The paper's central claim is "training-free generation from small alignments." Without learned baselines, reviewers will say "this is interesting, but how does it compare to EvoDiff / MSA Transformer / VAE?" This experiment answers that question directly.

---

## Baseline Models

### Tier 1: Must-have (directly comparable, scarce-data relevant)

1. **Profile HMM sampling** (HMMER3)
   - The classical bioinformatics baseline for sequence generation from an MSA
   - `hmmbuild` → `hmmemit --sample`
   - Pure sequence model, no structure, no deep learning
   - Should work at any K but produces limited diversity
   - **Dependency:** HMMER3 (`brew install hmmer` or download from hmmer.org)

2. **EvoDiff** (Microsoft, 2023)
   - Diffusion model for protein sequences, can condition on MSA
   - MSA-conditioned mode (`oa_dm_640M`) is the relevant comparison
   - Pretrained on UniRef50; fine-tuning optional
   - **Dependency:** Python, PyTorch, `evodiff` package
   - **Key question:** How does it perform when the conditioning MSA has only 30–100 sequences?

3. **MSA Transformer** (Meta/FAIR, 2021) — generative sampling
   - Not natively generative, but can be used for masked language model (MLM) sampling
   - Iterative masking + infilling to generate sequences from MSA context
   - **Dependency:** Python, PyTorch, `esm` package

### Tier 2: Nice-to-have (if time permits)

4. **ProtGPT2** (Ferruz et al., 2022)
   - Autoregressive protein language model
   - Unconditional generation (no MSA conditioning) — less directly comparable
   - Can measure: does it generate family-like sequences at all without family context?
   - **Dependency:** Python, HuggingFace `transformers`

5. **Sequence VAE**
   - Train a simple VAE on the MSA (one-hot encoded sequences)
   - Directly tests the "learned model on small data" failure mode
   - Can implement in Julia (Flux.jl) or Python (PyTorch)
   - **Key point:** A VAE trained on K=50 sequences will overfit catastrophically

6. **ArDCA** (Trinquier et al., 2021)
   - Autoregressive model based on Direct Coupling Analysis
   - Specifically designed for protein families
   - Good scarce-data baseline (uses statistical coupling, not deep learning)
   - **Dependency:** Julia package `ArDCA.jl` or Python reimplementation

## Experimental Design

### Head-to-head comparison protocol
For each baseline and each of the 7 Pfam families:
1. Provide the same cleaned MSA as input
2. Generate S=150 sequences (matching SA output count)
3. Decode/process outputs to amino acid sequences
4. Evaluate using the same 3 metrics: KL_AA, Novelty, SeqID

### Scarce-data stress test (using scaling study data)
For the WW domain at K ∈ {20, 50, 100, 200, 400}:
1. Run each baseline on the same 5 random subsets used in the scaling study
2. Track how each method degrades as K decreases
3. This is the money plot: "SA remains stable while learned models collapse at K < 50"

## Implementation Details

### Script structure
```
baselines/
├── PLAN.md                      # this file
├── run_baselines.jl             # main orchestration script (Julia)
├── run_hmm_baseline.jl          # HMMER3 baseline (Julia + shell calls)
├── run_evodiff_baseline.py      # EvoDiff baseline (Python)
├── run_msa_transformer.py       # MSA Transformer MLM sampling (Python)
├── run_vae_baseline.jl          # Simple VAE (Julia/Flux)  [Tier 2]
├── data/
│   └── <PF_ID>/
│       ├── hmm_sequences.fasta
│       ├── evodiff_sequences.fasta
│       ├── msat_sequences.fasta
│       └── vae_sequences.fasta
├── results/
│   ├── baseline_comparison.csv
│   └── scaling_baseline_comparison.csv
└── figs/
    ├── baseline_comparison_composite.pdf
    └── scaling_baseline_comparison.pdf
```

### HMM baseline (`run_hmm_baseline.jl`)
```julia
# For each family:
# 1. Build HMM from the cleaned alignment
run(`hmmbuild --amino data/$pfam_id/family.hmm data/$pfam_id/alignment.sto`)

# 2. Emit 150 sequences
run(`hmmemit -N 150 --seed 42 -o data/$pfam_id/hmm_sequences.fasta data/$pfam_id/family.hmm`)

# 3. Read and evaluate
hmm_seqs = parse_fasta("data/$pfam_id/hmm_sequences.fasta")
# ... compute KL, novelty, seqid
```

### EvoDiff baseline (`run_evodiff_baseline.py`)
```python
# Key: MSA-conditioned generation
from evodiff.pretrained import OA_DM_MSA
model, collater, tokenizer, scheme = OA_DM_MSA()

# Load the family MSA, generate 150 sequences
# Important: EvoDiff expects MSA input in specific format
# Save to FASTA for Julia evaluation
```

### MSA Transformer baseline (`run_msa_transformer.py`)
```python
import esm
model, alphabet = esm.pretrained.esm_msa1b_t12_100M_UR50S()

# MLM sampling strategy:
# 1. Load MSA
# 2. Add a row of <mask> tokens
# 3. Forward pass → sample from softmax at each masked position
# 4. Iterate: randomly re-mask 15% of positions, re-sample (Gibbs-like)
# 5. After N iterations, decode the generated row
```

### Evaluation (in `run_baselines.jl`)
```julia
# For each baseline method and family:
# 1. Read generated FASTA
# 2. Compute KL_AA against stored sequences
# 3. Compute novelty (1 - max cosine similarity in PCA space)
# 4. Compute SeqID (nearest sequence identity)
# 5. Aggregate into comparison table
```

## Expected Results & Story

### The scarce-data argument (this is the key narrative):
- **K > 200:** Learned models (EvoDiff, MSA-T) likely match or beat SA on KL
- **K = 50–100:** Learned models degrade; SA remains stable
- **K < 50:** Learned models fail (overfit / refuse to generate); SA still works
- **HMM:** Works at any K but produces very low novelty (samples from the HMM's probability model, which is close to a profile — similar to our bootstrap baseline)

### Expected ranking at K=50:
1. SA generation: low KL, high novelty (our method)
2. HMM emit: low-moderate KL, near-zero novelty (boring but correct)
3. EvoDiff: moderate-high KL (insufficient conditioning data)
4. MSA Transformer: high KL (MLM sampling from small MSA is noisy)
5. VAE: very high KL (catastrophic overfitting on 50 training sequences)

## Metrics to Report

| Metric | Reported as |
|--------|-------------|
| KL_AA | mean ± SE (bootstrap) |
| Novelty | mean ± SE (per-chain) |
| SeqID | mean ± SE (per-chain) |
| Training time | wall-clock seconds (0 for SA, varies for learned) |
| Generation time | wall-clock seconds per 150 sequences |

## Paper Integration

- **Results section:** New subsection: "Comparison with learned generative models"
- **Experiments section:** Paragraph on baseline setup (model versions, hyperparameters)
- **Figure:** Grouped bar chart or table comparing all methods across families
- **Key sentence:** "At family sizes below ~100 sequences, SA generation matches or exceeds learned models that require orders of magnitude more compute, while providing exact score-function guarantees."

## Risks & Mitigations

1. **EvoDiff MSA-conditioned mode may not work well on small MSAs:** This is actually our point — document the failure mode honestly
2. **Model versions / reproducibility:** Pin exact model checkpoints (e.g., `oa_dm_640M`, `esm_msa1b_t12_100M_UR50S`)
3. **Compute requirements:** EvoDiff and MSA-T need GPU. Plan for ~4 hours GPU time total.
4. **Unfair comparison concerns:** We're comparing a training-free method to pretrained models. Acknowledge that pretrained models leverage UniRef50 (millions of sequences). Our point is about the conditioning/fine-tuning data, not pretraining data.
5. **HMM baseline too strong:** If HMMER3 produces great KL scores (likely), emphasize that it produces near-zero novelty — it's a density model, not a generative model in the creative sense.

## Priority & Dependencies

1. **Do first:** HMM baseline — easiest to implement, strongest classical baseline, Julia-native
2. **Do second:** EvoDiff — the most relevant modern comparison
3. **Do third:** MSA Transformer — complementary to EvoDiff
4. **If time:** VAE (demonstrates overfitting), ProtGPT2 (unconditional baseline)

**Dependencies:** Requires `hmmer` installed, Python environment with `evodiff` and `esm` packages, GPU access for EvoDiff/MSA-T.
