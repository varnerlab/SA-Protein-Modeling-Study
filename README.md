# Training-Free Generation of Protein Sequences from Small Family Alignments via Stochastic Attention

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20835924.svg)](https://doi.org/10.5281/zenodo.20835924)

This repository contains the code, data, and manuscript source for the paper *Training-Free Generation of Protein Sequences from Small Family Alignments via Stochastic Attention*. The method treats the energy function of a modern Hopfield network, constructed directly from a seed alignment, as a Boltzmann density and samples from it using the Unadjusted Langevin Algorithm (ULA). The resulting stochastic attention sampler generates novel, compositionally faithful protein sequences that are predicted to fold into the correct three-dimensional architecture, all without training, external data, or GPU resources. The complete pipeline, from a seed alignment to decoded amino acid sequences, is described in Algorithm S1 of the SI Appendix.

## Repository organization

The repository has two top-level directories: `code/` for all computational work (Julia and Python) and `paper/` for the LaTeX manuscript. The eight Pfam seed alignments used throughout the study are pre-cached under `code/pfam-families/data/` and do not need to be downloaded. Each family subdirectory (e.g., `PF00076/`) contains the Stockholm-format seed alignment (`PF00076_seed.sto`) along with the FASTA files produced by the core experiment (SA-generated sequences, SA-retrieval sequences, and stored sequences). Baseline-generated sequences are stored separately under `code/baselines/data/`, organized by family, with one FASTA file per method (HMM emit, EvoDiff, MSA Transformer). Structure validation results, including predicted PDB files and TM-score metrics, live under `code/structure-validation/data/`.

The core library lives in four source files under `code/src/`. `Compute.jl` implements the stochastic attention sampler itself: the `sample` function takes a memory matrix, an initial state, and a number of iterations, and returns the full Langevin trajectory. `Protein.jl` provides the protein data pipeline from end to end, including Stockholm and FASTA parsing, alignment cleaning, one-hot encoding, PCA dimensionality reduction with unit-norm projection, phase-transition detection via entropy inflection, and inverse PCA with argmax decoding. `Utilities.jl` supplies diagnostic functions (cosine similarity, Hopfield energy, attention entropy), and `Data.jl` generates synthetic memory matrices for unit testing. Every Julia script in the repository begins by calling `Include.jl`, which activates the project environment from `Project.toml`, installs any missing dependencies (Plots, StatsPlots, Distributions, NNlib, JLD2, DataFrames, CSV, HTTP, HypothesisTests, MultivariateStats, and others), and loads these four modules.

The experiments are organized into subdirectories under `code/`, each containing a self-contained script that can be run from the command line. The central experiment is `pfam-families/run_pfam_families.jl`, which builds the memory matrix for each of the eight families, detects the critical temperature, runs 30-chain Langevin sampling, decodes the generated sequences, and computes all evaluation metrics (KL divergence, novelty, diversity, sequence identity, valid residue fraction). Supporting scripts in the same directory handle the beta-defensin biophysical analysis (`analyze_defensin_biophysics.jl`), the main results table with standard errors (`compute_table_with_se.jl`), and figure regeneration from cached results (`regenerate_figures.jl`).

The scaling and temperature analyses build on this foundation. `scaling-study/run_scaling_study.jl` subsamples the WW domain from K=420 down to K=20 with five independent replicates at each size, and `multifamily-scaling/run_multifamily_scaling.jl` replicates the K=20 condition on SH3, Kunitz, and zf-C2H2 to confirm that scaling robustness is not family-specific. The critical temperature prediction is handled by three scripts in `beta-prediction/`: `run_beta_prediction.jl` fits the regression model and performs bootstrap uncertainty quantification, `analytical_beta_star.jl` carries out the full analytical investigation (Gaussian mean-field predictions, stored-pattern probe analysis, cross-similarity decomposition), and `investigate_bifurcation.jl` examines the Jacobian bifurcation condition at the critical temperature.

Four baselines are compared against stochastic attention. `baselines/run_hmm_baseline.jl` builds profile HMMs and emits sequences via HMMER3. `baselines/run_evodiff_baseline.py` generates sequences with the MSA-conditioned EvoDiff diffusion model. `baselines/run_msa_transformer_baseline.py` performs iterative masked language model sampling with the MSA Transformer. `baselines/run_potts_baseline.py` fits a Potts model by pseudo-likelihood maximization and generates sequences by Gibbs sampling (implemented from scratch using only NumPy and SciPy, with no external DCA packages). After all baselines have been run, `baselines/analyze_baselines_stats.jl` computes per-chain metrics, runs Wilcoxon rank-sum tests, and generates the composite comparison figure.

Structure validation proceeds in three stages. First, `structure-validation/run_structure_validation.jl` predicts structures for generated and stored sequences using the ESMFold REST API, downloads reference PDB structures (1FXL, 1SHG, 1PIN, 1BPI, 1ZAA, 1BE9, 1ATP, 1E4S) from the RCSB, and computes TM-scores using TM-align (a pre-compiled macOS binary is included at `structure-validation/bin/TMalign`; other platforms should download the binary from the [TM-align website](https://zhanggroup.org/TM-align/)). Second, for AlphaFold2 cross-validation, `prepare_colabfold_input.jl` prepares FASTA input files, and the Google Colab notebook `ColabFold_AF2_Predictions.ipynb` runs the actual predictions on a GPU runtime. The predicted PDB files must then be downloaded from Google Drive and placed under `structure-validation/colabfold_output/`. Third, `analyze_af2_results.jl` and `analyze_structure_results.jl` combine the ESMFold and AlphaFold2 results, run statistical comparisons, and generate the structure validation figures. Publication-quality structure renderings are produced by `render_structure_figure.py` and `render_sequence_figure.py`, which require PyMOL.

The remaining analyses each occupy their own directory. `covariation/run_covariation_analysis.jl` computes pairwise mutual information matrices for stored and generated sequences and measures MI preservation across methods; `covariation/plot_covariation.jl` generates the corresponding figures. `esm-validation/score_sequences.py` computes ESM2-650M pseudo-perplexity scores (the model weights, approximately 2.5 GB, are downloaded automatically on first run). `dms-validation/analyze_dms.py` cross-references generated substitutions against deep mutational scanning data from ProteinGym; the DMS datasets covering four families (DLG4/PDZ, GRB2/SH3, PABP/RRM, WW domain) are pre-cached in `dms-validation/data/`. `timing/run_timing_benchmark.jl` measures wall-clock generation times for SA and HMM emit. `sampling-diagnostics/run_sampling_diagnostics.jl` computes autocorrelation, effective sample size, and burn-in convergence diagnostics. `random-init/run_random_init.jl` compares stored-pattern and random-sphere initialization. `pca-sensitivity/run_pca_sensitivity.jl` sweeps the PCA variance retention threshold from 50% to 99%.

The manuscript source is in `paper/`. The main file `Paper_v1.tex` inputs section files from `paper/sections/` (introduction, results, discussion, theory, experiments, appendix), and references are in `References_v1.bib`. Running `./Build.sh` from the `paper/` directory executes the standard four-pass LaTeX compilation cycle and produces `Paper_v1.pdf`.

## Getting started

The Julia experiments were developed and tested with [Julia](https://julialang.org/downloads/) 1.11. The simplest installation route on all platforms is [juliaup](https://github.com/JuliaLang/juliaup) (`curl -fsSL https://install.julialang.org | sh` on macOS/Linux, or `winget install Julia` on Windows). On the first run, `Include.jl` will download and precompile all Julia dependencies automatically. The notebooks and scripts can be opened directly in Visual Studio Code with the Julia and Jupyter extensions, which is the workflow used during development.

Several baseline and validation scripts are written in Python (3.9+) and depend on specific packages. Each script prints its own install command if a required package is missing, but the complete list is:

| Script(s) | Install command |
|---|---|
| EvoDiff baseline | `pip install evodiff torch "setuptools<81"` |
| MSA Transformer baseline, ESM2 scoring | `pip install fair-esm torch` |
| Potts baseline, DMS validation | `pip install numpy scipy` |

Note that the correct package name for the Meta ESM models is `fair-esm`, not `esm`. The EvoDiff package requires `setuptools<81` due to a packaging constraint. Both the MSA Transformer (~330 MB) and ESM2-650M (~2.5 GB) model weights are downloaded automatically on first use.

Two external command-line tools are also required for specific experiments: [HMMER3](http://hmmer.org/) for the HMM baseline (`brew install hmmer` on macOS) and [TM-align](https://zhanggroup.org/TM-align/) for structure validation (a macOS binary is included; download for other platforms). AlphaFold2 cross-validation runs on [Google Colab](https://github.com/sokrypton/ColabFold) via the included notebook and requires a GPU runtime. [PyMOL](https://pymol.org/) is optional and needed only for structure gallery rendering. Building the paper PDF requires a LaTeX distribution with `pdflatex` and `bibtex` ([TeX Live](https://tug.org/texlive/) or [MiKTeX](https://miktex.org/)).

To reproduce the experiments, clone the repository and run the scripts from their respective directories:

```bash
git clone https://github.com/varnerlab/SA-Protein-Modeling-Study.git
cd SA-Protein-Modeling-Study/code

# Core experiment: generate SA sequences and compute metrics for all eight families
cd pfam-families && julia run_pfam_families.jl

# Baselines (run after the core experiment produces stored sequence files)
cd ../baselines
julia run_hmm_baseline.jl
python run_evodiff_baseline.py
python run_msa_transformer_baseline.py
python run_potts_baseline.py
julia analyze_baselines_stats.jl

# Structure validation (requires internet for ESMFold API and RCSB PDB)
cd ../structure-validation && julia run_structure_validation.jl
# For AF2: julia prepare_colabfold_input.jl, then run the Colab notebook,
# download PDBs to colabfold_output/, then julia analyze_af2_results.jl

# Paper
cd ../../paper && ./Build.sh
```

The experiments can be run independently in any order, with the exception that the baselines and structure validation require the SA-generated and stored sequence FASTA files produced by the core experiment. Each script writes its output (CSVs, figures, FASTA files) to `data/`, `results/`, or `figs/` within its own directory, and copies figures to `paper/figs/` when the paper directory is present. All random seeds are fixed throughout the codebase, so running the scripts from scratch will produce identical results regardless of platform.

## Protein families

The eight Pfam families span a ten-fold range of family sizes and a seven-fold range of sequence lengths:

| Family | Pfam ID | K (sequences) | L (alignment length) | d (PCA dim) |
|---|---|---|---|---|
| RRM | PF00076 | 68 | 71 | 59 |
| SH3 | PF00018 | 55 | 48 | 46 |
| WW | PF00397 | 420 | 31 | 186 |
| Kunitz | PF00014 | 99 | 53 | 80 |
| zf-C2H2 | PF00096 | 151 | 23 | 106 |
| PDZ | PF00595 | 44 | 83 | 37 |
| Pkinase | PF00069 | 37 | 262 | 34 |
| Defensin_beta | PF00711 | 45 | 36 | 37 |

## License

This project is released under the [MIT License](LICENSE).
