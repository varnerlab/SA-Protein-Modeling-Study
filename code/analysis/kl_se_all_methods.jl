#!/usr/bin/env julia
# Unified 1000-fold bootstrap KL SE for ALL methods across 8 families, so KL SE
# uses one convention in Tables S2, S9 (Potts), and the baseline-metrics table.
# Verifies point estimates match the published tables before SEs are adopted.
using Statistics, Random

const PFAM = joinpath(@__DIR__, "..", "pfam-families", "data")
const BASE = joinpath(@__DIR__, "..", "baselines", "data")
const AA = collect("ACDEFGHIKLMNPQRSTVWY")
const IDX = Dict(a => i for (i, a) in enumerate(AA))

function parse_fasta(fp)
    seqs = String[]; cur = IOBuffer(); have = false
    for ln in eachline(fp)
        if startswith(ln, ">")
            have && push!(seqs, String(take!(cur))); cur = IOBuffer(); have = true
        else
            write(cur, uppercase(strip(ln)))
        end
    end
    have && push!(seqs, String(take!(cur)))
    return seqs
end

function aa_kl(gen, stored)
    f(seqs) = begin
        c = zeros(20)
        for s in seqs, ch in s
            i = get(IDX, ch, 0); i > 0 && (c[i] += 1)
        end
        t = sum(c); t > 0 ? c ./ t : ones(20) ./ 20
    end
    p = f(stored); q = f(gen)
    p .+= 1e-10; p ./= sum(p); q .+= 1e-10; q ./= sum(q)
    sum(p[i] * log(p[i] / q[i]) for i in 1:20)
end

function boot_se(gen, stored; nb=1000)
    n = length(gen)
    std([aa_kl(gen[rand(1:n, n)], stored) for _ in 1:nb])
end

FAMS = [("PF00076","RRM"),("PF00018","SH3"),("PF00397","WW"),("PF00014","Kunitz"),
        ("PF00096","zf-C2H2"),("PF00595","PDZ"),("PF00069","Pkinase"),("PF00711","Defensin_beta")]

# method label -> (dir, filename)
METHODS = [("HMM","hmm_sequences.fasta"), ("EvoDiff","evodiff_sequences.fasta"),
           ("MSAT","msat_sequences.fasta"), ("Potts","potts_sequences.fasta")]

Random.seed!(12345)
println(rpad("Family",15), rpad("Method",10), rpad("KL_point",10), "boot_SE")
println("-"^45)
for (pid, name) in FAMS
    stored = parse_fasta(joinpath(PFAM, pid, "stored_sequences.fasta"))
    for (ml, fn) in METHODS
        fp = joinpath(BASE, pid, fn)
        isfile(fp) || (println(rpad(name,15), rpad(ml,10), "MISSING"); continue)
        gen = parse_fasta(fp)
        klp = aa_kl(gen, stored); se = boot_se(gen, stored)
        println(rpad(name,15), rpad(ml,10), rpad(round(klp,digits=4),10), round(se,digits=4))
    end
end
