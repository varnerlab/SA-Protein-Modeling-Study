"""
    nearest_cosine_similarity(ξ::Vector{Float64}, X::Matrix{Float64}) -> Float64

Compute the maximum cosine similarity between the state vector `ξ` and the columns of the memory matrix `X`.
Returns the cosine similarity to the nearest stored pattern:

``\\text{sim}(\\boldsymbol{\\xi}, \\mathbf{X}) = \\max_{k} \\frac{\\mathbf{m}_k^\\top \\boldsymbol{\\xi}}{\\|\\mathbf{m}_k\\| \\, \\|\\boldsymbol{\\xi}\\|}``

### Arguments
- `ξ::Vector{Float64}`: State vector of dimension `d`.
- `X::Matrix{Float64}`: Memory matrix of size `d × K`, where each column is a stored pattern.

### Returns
- `Float64`: Maximum cosine similarity in `[-1, 1]`. Returns `0.0` if `ξ` is near-zero.
"""
function nearest_cosine_similarity(ξ::Vector{Float64}, X::Matrix{Float64})::Float64
    nξ = norm(ξ)
    nξ < 1e-12 && return 0.0  # guard against zero vector
    max_sim = -Inf
    for k in 1:size(X, 2)
        mk = @view X[:, k]
        sim_k = dot(mk, ξ) / (norm(mk) * nξ)
        max_sim = max(max_sim, sim_k)
    end
    return max_sim
end

"""
    hopfield_energy(ξ::Vector{Float64}, X::Matrix{Float64}, β::Float64) -> Float64

Compute the modern Hopfield energy at state `ξ`:

``E(\\boldsymbol{\\xi}) = \\tfrac{1}{2}\\|\\boldsymbol{\\xi}\\|^2 - \\tfrac{1}{\\beta}\\log\\sum_{k=1}^{K}\\exp(\\beta\\,\\mathbf{m}_k^\\top\\boldsymbol{\\xi})``

Uses the log-sum-exp trick for numerical stability.

### Arguments
- `ξ::Vector{Float64}`: State vector of dimension `d`.
- `X::Matrix{Float64}`: Memory matrix of size `d × K`.
- `β::Float64`: Inverse temperature (must be positive).

### Returns
- `Float64`: The energy value.
"""
function hopfield_energy(ξ::Vector{Float64}, X::Matrix{Float64}, β::Float64)::Float64
    logits = β .* (X' * ξ)
    logits_max = maximum(logits)
    lse = logits_max + log(sum(exp.(logits .- logits_max)))  # log-sum-exp (stable)
    return 0.5 * dot(ξ, ξ) - lse / β
end

"""
    attention_entropy(ξ::Vector{Float64}, X::Matrix{Float64}, β::Float64) -> Float64

Compute the Shannon entropy of the attention weights ``\\mathbf{p} = \\operatorname{softmax}(\\beta\\,\\mathbf{X}^\\top\\boldsymbol{\\xi})``:

``H(\\mathbf{p}) = -\\sum_{k=1}^{K} p_k \\log p_k``

The entropy serves as an order parameter: ``H \\approx \\log K`` at small ``\\beta`` (uniform/disordered)
and ``H \\to 0`` at large ``\\beta`` (concentrated/ordered).

### Arguments
- `ξ::Vector{Float64}`: State vector of dimension `d`.
- `X::Matrix{Float64}`: Memory matrix of size `d × K`.
- `β::Float64`: Inverse temperature (must be positive).

### Returns
- `Float64`: Shannon entropy in nats.
"""
function attention_entropy(ξ::Vector{Float64}, X::Matrix{Float64}, β::Float64)::Float64
    p = NNlib.softmax(β .* (X' * ξ))
    H = 0.0
    for pk in p
        if pk > 1e-30  # guard against log(0)
            H -= pk * log(pk)
        end
    end
    return H
end

"""
    sample_novelty(ξ::Vector{Float64}, X::Matrix{Float64}) -> Float64

Compute the novelty of a generated sample `ξ` with respect to the memory matrix `X`.
Novelty measures how far the sample is from the nearest stored pattern:

``\\text{Novelty}(\\boldsymbol{\\xi}) = 1 - \\max_{k} \\cos(\\hat{\\boldsymbol{\\xi}},\\, \\mathbf{m}_k)``

where ``\\hat{\\boldsymbol{\\xi}}`` is the L2-normalized sample. A novelty of `0` indicates 
that the sample is an exact copy of a stored pattern; higher values indicate greater departure 
from any single memory.

### Arguments
- `ξ::Vector{Float64}`: Generated sample vector of dimension `d`.
- `X::Matrix{Float64}`: Memory matrix of size `d × K`, where each column is a stored pattern.

### Returns
- `Float64`: Novelty score in `[0, 2]`. Typical values are in `[0, 1]`.
"""
function sample_novelty(ξ::Vector{Float64}, X::Matrix{Float64})::Float64
    return 1.0 - nearest_cosine_similarity(ξ, X)
end

"""
    sample_diversity(samples::Vector{Vector{Float64}}) -> Float64

Compute the diversity of a collection of generated samples, defined as the mean pairwise 
cosine distance:

``\\text{Diversity} = \\frac{2}{S(S-1)} \\sum_{i<j} \\left(1 - \\cos(\\boldsymbol{\\xi}_i,\\, \\boldsymbol{\\xi}_j)\\right)``

where ``\\cos(\\mathbf{a}, \\mathbf{b}) = \\mathbf{a}^\\top \\mathbf{b} / (\\|\\mathbf{a}\\| \\, \\|\\mathbf{b}\\|)`` and ``S`` is the number of samples.
A diversity of `0` means all samples are identical; higher values indicate greater spread across 
the sample collection.

### Arguments
- `samples::Vector{Vector{Float64}}`: Collection of `S` generated sample vectors, each of dimension `d`.

### Returns
- `Float64`: Mean pairwise cosine distance. Returns `0.0` if fewer than two samples are provided.
"""
function sample_diversity(samples::Vector{Vector{Float64}})::Float64
    S = length(samples)
    S < 2 && return 0.0
    
    total = 0.0
    count = 0
    for i in 1:(S-1)
        nᵢ = norm(samples[i])
        nᵢ < 1e-12 && continue
        for j in (i+1):S
            nⱼ = norm(samples[j])
            nⱼ < 1e-12 && continue
            cos_ij = dot(samples[i], samples[j]) / (nᵢ * nⱼ)
            total += (1.0 - cos_ij)
            count += 1
        end
    end
    return count > 0 ? total / count : 0.0
end

"""
    sample_quality(samples::Vector{Vector{Float64}}, X::Matrix{Float64}, β::Float64) -> Float64

Compute the quality of a collection of generated samples, defined as the mean Hopfield energy:

``\\bar{E} = \\frac{1}{S}\\sum_{i=1}^{S} E(\\boldsymbol{\\xi}_i)``

where ``E(\\boldsymbol{\\xi}) = \\tfrac{1}{2}\\|\\boldsymbol{\\xi}\\|^2 - \\tfrac{1}{\\beta}\\log\\sum_k \\exp(\\beta\\,\\mathbf{m}_k^\\top \\boldsymbol{\\xi})``.
Lower energy indicates that the samples lie closer to the memory manifold, producing more 
recognizable outputs.

### Arguments
- `samples::Vector{Vector{Float64}}`: Collection of `S` generated sample vectors, each of dimension `d`.
- `X::Matrix{Float64}`: Memory matrix of size `d × K`.
- `β::Float64`: Inverse temperature (must be positive).

### Returns
- `Float64`: Mean Hopfield energy across the sample collection.
"""
function sample_quality(samples::Vector{Vector{Float64}}, X::Matrix{Float64}, β::Float64)::Float64
    S = length(samples)
    S == 0 && return 0.0
    return sum(hopfield_energy(s, X, β) for s in samples) / S
end
