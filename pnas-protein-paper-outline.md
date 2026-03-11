# PNAS Paper Outline: Training-Free Protein Sequence Generation via Stochastic Attention

**Target journal:** PNAS (Direct submission)

**Working title:** Training-Free Generation of Protein Sequences from Small Family Alignments via Stochastic Attention

**One-sentence pitch:** Attention-based memory retrieval, viewed as an energy landscape, provides an exact score function for generating biologically valid protein sequences from as few as 50 family members, without training.

---

## 1. Introduction / Significance Statement

- Protein sequence generation is critical for drug design, enzyme engineering, and understanding evolution
- Existing generative models (EvoDiff, ProtGPT2, ESM-based, VAEs) require large training sets (thousands to millions of sequences)
- Most protein families have tens to hundreds of known members -- the scarce-data regime
- We show that the modern Hopfield energy provides a closed-form score function, turning standard attention into a training-free sampler
- Key question: **what can you generate from 50-200 sequences with no training?**

## 2. Theory

### 2a. From Hopfield Energy to Score Function
- Modern Hopfield energy and its gradient (attention residual)
- Connection to Langevin dynamics: the score is exact, no approximation
- The stochastic attention update (from current paper)

### 2b. Critical Temperature from MSA Statistics
- $\beta^*$ prediction from family-specific properties
- How conservation entropy, effective number of sequences, and similarity spectrum determine $\beta^*$
- Can we predict $\beta^*$ directly from MSA statistics without grid search?
- Contrast with random-pattern scaling $\beta^* \sim \sqrt{d}$

### 2c. Convergence Guarantees
- $W_2$-geometric convergence in convex regime
- MALA acceptance as discretization diagnostic
- What guarantees carry over to the protein setting?

## 3. Experiments

### 3a. Multiple Pfam Families
- Expand beyond RRM to families with varying:
  - Size (K = 20, 50, 100, 200, 500)
  - Conservation level (highly conserved vs. diverse)
  - Sequence length
- Families to consider: RRM, SH3, WW domain, Kunitz, zinc fingers, PDZ, kinase domains

### 3b. Antimicrobial Peptides (AMPs)
- Therapeutically relevant
- Naturally scarce (small families, short sequences)
- Well-characterized validation criteria (charge, hydrophobicity, secondary structure)
- Existing databases (APD, DRAMP) for benchmarking

### 3c. Scaling Study: Quality vs. Family Size
- Fix a large family, subsample K = {20, 50, 100, 200, 500}
- Track generation quality metrics as K shrinks
- Show where learned baselines break down vs. SA

### 3d. Baselines
- **Protein-specific learned models:** DeepSequence, EvoDiff, ProtGPT2, MSA Transformer (generative mode)
- **General generative models:** VAE, DDPM (as in current paper)
- **Simple baselines:** profile HMM sampling, bootstrap, Gaussian perturbation
- Head-to-head on same families, same K

## 4. Validation

### 4a. Sequence-Level
- Amino acid composition (KL divergence)
- Per-position conservation profiles
- Pairwise mutual information (coupling fidelity)
- Higher-order statistics (three-body correlations)
- Pfam HMM E-values and pass rates

### 4b. Structure-Level
- AlphaFold2 or ESMFold structure prediction on generated sequences
- pLDDT scores as confidence metric
- TM-score to known family structures
- Contact map comparison

### 4c. Function-Level (if feasible)
- Molecular dynamics stability simulations on top candidates
- Wet-lab collaboration for AMP candidates (synthesis + activity assays)
- Phylogenetic placement: do generated sequences sit within the family tree?

## 5. Analysis / Discussion

### 5a. When Does Training-Free Beat Learned?
- Characterize the crossover point: at what K do learned models start winning?
- Relate to effective number of parameters vs. family size
- The "nothing to train on" regime

### 5b. What Does $\beta^*$ Tell Us About a Protein Family?
- Biological interpretation of the critical temperature
- Connection to evolutionary constraints and functional selection
- Families with low $\beta^*$ (diverse, fewer constraints) vs. high $\beta^*$ (conserved, tight constraints)

### 5c. Limitations
- Assumes aligned sequences (MSA required)
- Does not model insertions/deletions
- PCA dimensionality reduction may lose information
- No explicit structural or functional constraints

## 6. Conclusion
- Training-free generation is viable for scarce protein families
- The score function is free -- it's the attention map
- Practical tool for protein engineers working with small datasets

---

## Key Differentiators from Current Paper
1. Multiple families (not just RRM)
2. Scaling study (quality vs. K)
3. Structure-level validation (AlphaFold2)
4. Protein-specific baselines (not just VAE/DDPM)
5. $\beta^*$ prediction from MSA statistics (new theory)
6. Biological interpretation of critical temperature

## Open Questions / Risks
- Can we get AlphaFold2 to produce confident structures for generated sequences?
- How much does PCA dimensionality reduction hurt for longer proteins?
- Will the method scale to families with L > 200 residues?
- Reviewer concern: "this is just sampling from a Boltzmann distribution" -- need to emphasize the connection to attention and the training-free aspect

## Timeline Estimate
- TODO: discuss tonight
