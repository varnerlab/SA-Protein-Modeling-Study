#!/usr/bin/env python3
"""
Consensus-regression test, TM-score prong (issue #2). For each family, run the
bundled TM-align on each EXISTING SA-generated predicted structure vs the family
reference, parse the reference-length-normalized TM-score (the paper's TM), join
with each sequence's identity to the consensus, and compute the within-family
Spearman correlation. Tests whether higher TM tracks consensus proximity.
"""
import os, re, glob, csv, math, subprocess
import numpy as np

CODE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SV   = os.path.join(CODE, "structure-validation", "data")
TM   = os.path.join(CODE, "structure-validation", "bin", "TMalign")
CONS = os.path.join(CODE, "analysis", "data", "review_consensus_identity.csv")

FAMS = [("PF00076","RRM","1FXL_A"),("PF00018","SH3","1SHG_A"),("PF00397","WW","1PIN_A"),
        ("PF00014","Kunitz","1BPI_A"),("PF00096","zf-C2H2","1ZAA_C"),("PF00595","PDZ","1BE9_A"),
        ("PF00069","Pkinase","1ATP_E"),("PF00711","Defensin_beta","1E4S_A")]

cons = {}
with open(CONS) as f:
    for row in csv.DictReader(f):
        cons[(row["Family"], int(row["SeqIdx"]))] = float(row["ConsensusID"])

def spearmanr(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float); n = len(x)
    rx = np.argsort(np.argsort(x)).astype(float); ry = np.argsort(np.argsort(y)).astype(float)
    r = np.corrcoef(rx, ry)[0, 1]
    if n <= 2 or abs(r) >= 1.0: return r, float("nan")
    t = r * math.sqrt((n - 2) / (1 - r * r))
    return r, 2 * (1 - 0.5 * (1 + math.erf(abs(t) / math.sqrt(2))))

# parse TM normalized by Chain_2 (reference)
pat = re.compile(r"TM-score=\s*([0-9.]+).*Chain_2")

print(f"{'Family':14s} {'n':>3s} {'rho(TM,consID)':>15s} {'p':>8s}   {'TM range':>13s}  {'TMmean':>7s}")
print("-"*70)
allz_c, allz_t = [], []
for pid, name, ref in FAMS:
    refpdb = os.path.join(SV, pid, "reference", ref + ".pdb")
    d = os.path.join(SV, pid, "structures", "sa_generation")
    tms, cid = [], []
    for pdb in glob.glob(os.path.join(d, "sa_generation_*.pdb")):
        idx = int(re.search(r"_(\d+)\.pdb", pdb).group(1))
        if (name, idx) not in cons: continue
        out = subprocess.run([TM, pdb, refpdb], capture_output=True, text=True).stdout
        m = pat.search(out)
        if not m: continue
        tms.append(float(m.group(1))); cid.append(cons[(name, idx)])
    tms = np.array(tms); cid = np.array(cid)
    if len(tms) < 5:
        print(f"{name:14s} {len(tms):3d}  (insufficient)"); continue
    r, p = spearmanr(cid, tms)
    print(f"{name:14s} {len(tms):3d} {r:15.3f} {p:8.3f}   {tms.min():.3f}-{tms.max():.3f}  {tms.mean():7.3f}")
    allz_c.append((cid-cid.mean())/(cid.std()+1e-9)); allz_t.append((tms-tms.mean())/(tms.std()+1e-9))

zc = np.concatenate(allz_c); zt = np.concatenate(allz_t)
r, p = spearmanr(zc, zt)
print("-"*70)
print(f"Pooled within-family-standardized Spearman rho(TM, consID) = {r:.3f} (p = {p:.3g}, n = {len(zc)})")
