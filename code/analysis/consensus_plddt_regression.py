#!/usr/bin/env python3
"""
Consensus-regression test (issue #2, feasible part; no new structure prediction).

For each family, parse the EXISTING ESMFold per-sequence PDBs of SA-generated
sequences, compute per-sequence mean pLDDT from CA B-factors, join with each
sequence's identity to the family consensus (review_consensus_identity.csv),
and compute the within-family Spearman correlation between consensus-identity
and pLDDT. A weak/no correlation argues against the hypothesis that the
structure predictor's confidence is merely a function of consensus proximity.
"""
import os, re, glob, csv, math
import numpy as np

def spearmanr(x, y):
    """Spearman rho via Pearson-on-ranks; two-sided p from t-approx (df=n-2)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    n = len(x)
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    r = np.corrcoef(rx, ry)[0, 1]
    if n <= 2 or abs(r) >= 1.0:
        return r, float("nan")
    t = r * math.sqrt((n - 2) / (1 - r * r))
    # two-sided p via normal approx to Student-t (adequate for n~50; approximate)
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(t) / math.sqrt(2))))
    return r, p

CODE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SV   = os.path.join(CODE, "structure-validation", "data")
CONS = os.path.join(CODE, "analysis", "data", "review_consensus_identity.csv")

FAMILIES = [("PF00076","RRM"),("PF00018","SH3"),("PF00397","WW"),
            ("PF00014","Kunitz"),("PF00096","zf-C2H2"),("PF00595","PDZ"),
            ("PF00069","Pkinase"),("PF00711","Defensin_beta")]

# consensus identity: {(family, seqidx) -> consID}
cons = {}
with open(CONS) as f:
    for row in csv.DictReader(f):
        cons[(row["Family"], int(row["SeqIdx"]))] = float(row["ConsensusID"])

def mean_plddt(pdb):
    vals = []
    with open(pdb) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                vals.append(float(line[60:66]))
    if not vals:
        return None
    m = np.mean(vals)
    return m * 100.0 if m <= 1.0 else m   # B-factor stored on 0-1 scale here

print(f"{'Family':14s} {'n':>3s} {'Spearman r':>11s} {'p':>9s}   {'pLDDT range':>14s}  {'consID range':>14s}")
print("-"*72)
allz_c, allz_p = [], []
for pid, name in FAMILIES:
    d = os.path.join(SV, pid, "structures", "sa_generation")
    plddt, cid = [], []
    for pdb in glob.glob(os.path.join(d, "sa_generation_*.pdb")):
        m = re.search(r"sa_generation_(\d+)\.pdb", os.path.basename(pdb))
        idx = int(m.group(1))
        key = (name, idx)
        if key not in cons:
            continue
        p = mean_plddt(pdb)
        if p is None:
            continue
        plddt.append(p); cid.append(cons[key])
    plddt = np.array(plddt); cid = np.array(cid)
    if len(plddt) < 5:
        print(f"{name:14s} {len(plddt):3d}  (insufficient)")
        continue
    r, pval = spearmanr(cid, plddt)
    print(f"{name:14s} {len(plddt):3d} {r:11.3f} {pval:9.3f}   "
          f"{plddt.min():5.1f}-{plddt.max():5.1f}  {cid.min():.3f}-{cid.max():.3f}")
    # standardized within-family for a pooled estimate
    allz_c.append((cid - cid.mean())/ (cid.std()+1e-9))
    allz_p.append((plddt - plddt.mean())/(plddt.std()+1e-9))

zc = np.concatenate(allz_c); zp = np.concatenate(allz_p)
r, pval = spearmanr(zc, zp)
print("-"*72)
print(f"Pooled within-family-standardized Spearman r = {r:.3f} (p = {pval:.3g}, n = {len(zc)})")
