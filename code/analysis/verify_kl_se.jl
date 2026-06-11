#!/usr/bin/env julia
# Recompute the 1000-fold bootstrap KL SE for SA-generation across all 8 families
# (matching compute_table_with_se.jl) to confirm Table S2 values and the
# zf-C2H2 SE>mean case, and to provide one canonical KL-SE convention.
using Statistics, Random

const PFAM = joinpath(@__DIR__, "..", "pfam-families", "data")
const AA = collect("ACDEFGHIKLMNPQRSTVWY")
const IDX = Dict(a => i for (i, a) in enumerate(AA))

function parse_fasta(fp)
    seqs = String[]; cur = IOBuffer(); have = false
    for ln in eachline(fp)
        if startswith(ln, ">")
            have && push!(seqs, String(take!(cur)))
            cur = IOBuffer(); have = true
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

FAMS = [("PF00076","RRM"),("PF00018","SH3"),("PF00397","WW"),("PF00014","Kunitz"),
        ("PF00096","zf-C2H2"),("PF00595","PDZ"),("PF00069","Pkinase"),("PF00711","Defensin_beta")]

Random.seed!(12345)
println("Family          KL_point   boot_SE   SE/mean")
println("-"^48)
for (pid, name) in FAMS
    stored = parse_fasta(joinpath(PFAM, pid, "stored_sequences.fasta"))
    gen = parse_fasta(joinpath(PFAM, pid, "sa_generation_sequences.fasta"))
    klp = aa_kl(gen, stored)
    n = length(gen)
    boots = [aa_kl(gen[rand(1:n, n)], stored) for _ in 1:1000]
    se = std(boots)
    println(rpad(name, 15), " ", rpad(round(klp, digits=4), 9), " ",
            rpad(round(se, digits=4), 9), " ", round(se/klp, digits=2))
end
