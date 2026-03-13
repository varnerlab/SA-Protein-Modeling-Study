#!/usr/bin/env python3
"""
Render a sequence-level analysis figure for the Kunitz (PF00014) family.

Panel layout:
  Top:    Sequence alignment -- reference, consensus, top 5 SA-generated
          Color-coded by per-position conservation in the stored MSA.
  Bottom: Per-position entropy comparison (stored MSA vs SA-generated)

Output: figs/sequence_analysis_kunitz.pdf (also copied to paper/figs/)
"""

import os, json, math
import numpy as np
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# ── Paths ─────────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PFAM_DATA = os.path.join(BASE_DIR, "..", "pfam-families", "data", "PF00014")
AF2_DIR = os.path.join(BASE_DIR, "colabfold_output", "PF00014", "sa_generation")
FIG_DIR = os.path.join(BASE_DIR, "figs")
PAPER_FIG_DIR = os.path.join(BASE_DIR, "..", "..", "paper", "figs")

# ── Parse FASTA ───────────────────────────────────────────────────────────
def parse_fasta(path):
    seqs = []
    header, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    seqs.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
    if header:
        seqs.append((header, "".join(seq)))
    return seqs


# ── Load data ─────────────────────────────────────────────────────────────
stored = parse_fasta(os.path.join(PFAM_DATA, "stored_sequences.fasta"))
generated = parse_fasta(os.path.join(PFAM_DATA, "sa_generation_sequences.fasta"))
stored_raw = [s[1] for s in stored]
gen_raw = [s[1] for s in generated]
L = max(len(s) for s in stored_raw)

# Find top 5 by AF2 pLDDT
plddt_scores = {}
for f in os.listdir(AF2_DIR):
    if "scores_rank_001" in f and f.endswith(".json"):
        seq_id = f.split("_scores_")[0]  # e.g. PF00014_sa_gen_48
        idx = int(seq_id.split("_")[-1])
        with open(os.path.join(AF2_DIR, f)) as fh:
            d = json.load(fh)
            plddt_scores[idx] = np.mean(d["plddt"])

top5 = sorted(plddt_scores, key=plddt_scores.get, reverse=True)[:5]
top5_seqs = {}
for idx in top5:
    header = f"sa_generation_{idx}"
    for g in generated:
        if g[0] == header:
            top5_seqs[idx] = g[1]
            break

# Reference 1BPI
aa3to1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
ref_pdb = os.path.join(BASE_DIR, "data", "PF00014", "reference", "1BPI_A.pdb")
ref_seq = []
seen = set()
with open(ref_pdb) as f:
    for line in f:
        if line.startswith("ATOM"):
            resname = line[17:20].strip()
            resnum = line[22:27].strip()
            if resnum not in seen and resname in aa3to1:
                seen.add(resnum)
                ref_seq.append(aa3to1[resname])
ref_seq_str = "".join(ref_seq)

# ── Per-position conservation ─────────────────────────────────────────────
conservation = []
consensus = []
stored_entropy = []
gen_entropy = []

for pos in range(L):
    # Stored MSA
    col = [s[pos] for s in stored_raw if pos < len(s) and s[pos] not in ".~-"]
    counts = Counter(col)
    n = len(col)
    if n > 0:
        top_aa, top_cnt = counts.most_common(1)[0]
        conservation.append(top_cnt / n)
        consensus.append(top_aa)
        H = -sum((c / n) * math.log2(c / n) for c in counts.values())
        stored_entropy.append(H)
    else:
        conservation.append(0)
        consensus.append("-")
        stored_entropy.append(0)

    # Generated
    col_g = [s[pos] for s in gen_raw if pos < len(s)]
    counts_g = Counter(col_g)
    n_g = len(col_g)
    if n_g > 0:
        H_g = -sum((c / n_g) * math.log2(c / n_g) for c in counts_g.values())
        gen_entropy.append(H_g)
    else:
        gen_entropy.append(0)

consensus_str = "".join(consensus)
cys_positions = [i for i, aa in enumerate(consensus) if aa == "C"]

# ── Align 1BPI to the consensus (simple: 1BPI is 58 aa, alignment is 53) ──
# The stored MSA is already aligned; 1BPI may have extra residues.
# Find 1BPI in stored sequences
ref_in_msa = None
for s in stored:
    if "BPT" in s[0].upper() or "1BPI" in s[0].upper():
        ref_in_msa = s[1]
        break
# If not found, use the ITIH4_HUMAN or just use consensus as reference
if ref_in_msa is None:
    # Use the stored sequence most similar to 1BPI PDB sequence
    # For now, just note we'll show consensus as the reference line
    ref_in_msa = None

# ══════════════════════════════════════════════════════════════════════════
# FIGURE
# ══════════════════════════════════════════════════════════════════════════

fig = plt.figure(figsize=(14, 7.5))

# Use gridspec: top 65% for alignment, bottom 35% for entropy
gs = fig.add_gridspec(2, 1, height_ratios=[1.8, 1], hspace=0.35,
                      left=0.08, right=0.95, top=0.94, bottom=0.08)

# ── Color scheme ──────────────────────────────────────────────────────────
# Conservation-based coloring for each cell
def get_cell_color(cons_frac, is_cys=False):
    """Return background color based on conservation level."""
    if is_cys:
        return "#FFD700"  # gold for cysteines
    if cons_frac >= 0.90:
        return "#c62828"  # dark red - highly conserved
    elif cons_frac >= 0.70:
        return "#ef6c00"  # orange - conserved
    elif cons_frac >= 0.50:
        return "#fff9c4"  # light yellow - moderate
    else:
        return "#e3f2fd"  # light blue - variable


def get_text_color(cons_frac, is_cys=False):
    if is_cys:
        return "black"
    if cons_frac >= 0.70:
        return "white"
    return "black"


# ── Panel A: Sequence alignment ──────────────────────────────────────────
ax_aln = fig.add_subplot(gs[0])
ax_aln.set_xlim(-0.5, L - 0.5)

# Rows: conservation bar, consensus, blank, top 5 SA sequences
row_labels = ["Conservation", "Consensus"]
row_labels += [f"SA gen {idx} (pLDDT={plddt_scores[idx]:.0f})" for idx in top5]
n_rows = len(row_labels)

ax_aln.set_ylim(-0.5, n_rows - 0.5)
ax_aln.invert_yaxis()

# Draw cells
cell_h = 0.85
for row_idx, label in enumerate(row_labels):
    if row_idx == 0:
        # Conservation bar - colored blocks, no text
        for pos in range(L):
            is_cys = consensus[pos] == "C"
            color = get_cell_color(conservation[pos], is_cys)
            rect = plt.Rectangle((pos - 0.45, row_idx - cell_h / 2),
                                 0.9, cell_h, facecolor=color,
                                 edgecolor="white", linewidth=0.3)
            ax_aln.add_patch(rect)
            # Show conservation percentage for highly conserved
            if conservation[pos] >= 0.90 or is_cys:
                ax_aln.text(pos, row_idx, f"{conservation[pos]:.0%}",
                           ha="center", va="center", fontsize=5,
                           color=get_text_color(conservation[pos], is_cys),
                           fontfamily="monospace", fontweight="bold")
    elif row_idx == 1:
        # Consensus sequence
        for pos in range(L):
            is_cys = consensus[pos] == "C"
            color = get_cell_color(conservation[pos], is_cys)
            rect = plt.Rectangle((pos - 0.45, row_idx - cell_h / 2),
                                 0.9, cell_h, facecolor=color,
                                 edgecolor="white", linewidth=0.3)
            ax_aln.add_patch(rect)
            ax_aln.text(pos, row_idx, consensus[pos],
                       ha="center", va="center", fontsize=7.5,
                       color=get_text_color(conservation[pos], is_cys),
                       fontfamily="monospace", fontweight="bold")
    else:
        # SA-generated sequence
        sa_idx = top5[row_idx - 2]
        seq = top5_seqs[sa_idx]
        for pos in range(min(L, len(seq))):
            aa = seq[pos]
            cons_aa = consensus[pos]
            is_cys = cons_aa == "C"

            if aa == cons_aa:
                # Match: gray background, dot
                rect = plt.Rectangle((pos - 0.45, row_idx - cell_h / 2),
                                     0.9, cell_h, facecolor="#f5f5f5",
                                     edgecolor="white", linewidth=0.3)
                ax_aln.add_patch(rect)
                ax_aln.text(pos, row_idx, "\u00B7",
                           ha="center", va="center", fontsize=9,
                           color="#999999", fontfamily="monospace")
            else:
                # Substitution: colored by conservation of that position
                if is_cys:
                    bg = "#FFD700"
                elif conservation[pos] >= 0.90:
                    bg = "#ffcdd2"  # light red - substitution at conserved site (rare!)
                elif conservation[pos] >= 0.70:
                    bg = "#ffe0b2"  # light orange
                else:
                    bg = "#bbdefb"  # blue - substitution at variable site
                rect = plt.Rectangle((pos - 0.45, row_idx - cell_h / 2),
                                     0.9, cell_h, facecolor=bg,
                                     edgecolor="white", linewidth=0.3)
                ax_aln.add_patch(rect)
                ax_aln.text(pos, row_idx, aa,
                           ha="center", va="center", fontsize=7.5,
                           color="black", fontfamily="monospace",
                           fontweight="bold")

# Position numbers along top
for pos in range(L):
    if (pos + 1) % 5 == 0:
        ax_aln.text(pos, -0.45, str(pos + 1), ha="center", va="bottom",
                   fontsize=6, color="gray")

# Row labels
for row_idx, label in enumerate(row_labels):
    ax_aln.text(-1, row_idx, label, ha="right", va="center", fontsize=7.5,
               fontfamily="sans-serif",
               fontweight="bold" if row_idx <= 1 else "normal")

# Cysteine disulfide bond annotations
# Kunitz disulfide bonds: C5-C55 -> in alignment: pos 2-52, 11-48 (approximate), 27-35 (approximate)
# Mark with brackets below
disulfide_pairs = [(1, 51), (10, 47), (26, 34)]  # 0-indexed
for p1, p2 in disulfide_pairs:
    y_bot = n_rows - 0.2
    ax_aln.annotate("", xy=(p1, y_bot), xytext=(p2, y_bot),
                    arrowprops=dict(arrowstyle="-", color="#B8860B",
                                   lw=1.5, connectionstyle="arc3,rad=-0.15"))
    mid = (p1 + p2) / 2
    ax_aln.text(mid, y_bot + 0.25, "S-S", ha="center", va="top",
               fontsize=5.5, color="#B8860B", fontweight="bold")

ax_aln.set_yticks([])
ax_aln.set_xticks([])
ax_aln.spines["top"].set_visible(False)
ax_aln.spines["right"].set_visible(False)
ax_aln.spines["bottom"].set_visible(False)
ax_aln.spines["left"].set_visible(False)
ax_aln.set_title("(A)  Kunitz domain: consensus vs. top SA-generated sequences",
                 fontsize=11, fontweight="bold", loc="left", pad=12)

# Legend
legend_elements = [
    mpatches.Patch(facecolor="#c62828", label="Highly conserved (>90%)"),
    mpatches.Patch(facecolor="#ef6c00", label="Conserved (70-90%)"),
    mpatches.Patch(facecolor="#fff9c4", label="Moderate (50-70%)"),
    mpatches.Patch(facecolor="#e3f2fd", label="Variable (<50%)"),
    mpatches.Patch(facecolor="#FFD700", label="Cysteine (disulfide)"),
    mpatches.Patch(facecolor="#bbdefb", label="SA substitution (variable site)"),
]
ax_aln.legend(handles=legend_elements, loc="upper right", fontsize=6,
             ncol=2, framealpha=0.9, bbox_to_anchor=(1.0, 1.15))

# ── Panel B: Per-position entropy ────────────────────────────────────────
ax_ent = fig.add_subplot(gs[1])

positions = np.arange(L)
bar_w = 0.35
ax_ent.bar(positions - bar_w / 2, stored_entropy, bar_w,
          color="#78909c", alpha=0.8, label="Stored MSA", edgecolor="white", linewidth=0.3)
ax_ent.bar(positions + bar_w / 2, gen_entropy, bar_w,
          color="#1565c0", alpha=0.8, label="SA generation", edgecolor="white", linewidth=0.3)

# Highlight cysteine positions
for cp in cys_positions:
    ax_ent.axvspan(cp - 0.5, cp + 0.5, alpha=0.15, color="#FFD700", zorder=0)

# Correlation annotation
r = np.corrcoef(stored_entropy, gen_entropy)[0, 1]
ax_ent.text(0.98, 0.95, f"r = {r:.3f}", transform=ax_ent.transAxes,
           ha="right", va="top", fontsize=9,
           bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

ax_ent.set_xlabel("Position", fontsize=9)
ax_ent.set_ylabel("Shannon entropy (bits)", fontsize=9)
ax_ent.set_title("(B)  Per-position entropy: stored MSA vs. SA-generated sequences",
                fontsize=11, fontweight="bold", loc="left")
ax_ent.set_xlim(-1, L)
ax_ent.set_ylim(0, max(max(stored_entropy), max(gen_entropy)) * 1.15)
ax_ent.legend(fontsize=8, loc="upper left")
ax_ent.tick_params(labelsize=8)

# Position ticks every 5
ax_ent.set_xticks([i for i in range(L) if (i + 1) % 5 == 0])
ax_ent.set_xticklabels([str(i + 1) for i in range(L) if (i + 1) % 5 == 0])
ax_ent.spines["top"].set_visible(False)
ax_ent.spines["right"].set_visible(False)

# ── Save ──────────────────────────────────────────────────────────────────
os.makedirs(FIG_DIR, exist_ok=True)
out_path = os.path.join(FIG_DIR, "sequence_analysis_kunitz.pdf")
fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
print(f"Saved: {out_path}")

# Copy to paper
os.makedirs(PAPER_FIG_DIR, exist_ok=True)
import shutil
paper_path = os.path.join(PAPER_FIG_DIR, "sequence_analysis_kunitz.pdf")
shutil.copy2(out_path, paper_path)
print(f"Copied to: {paper_path}")

# Print summary stats
print(f"\nTop 5 SA-generated sequences (by AF2 pLDDT):")
for idx in top5:
    seq = top5_seqs[idx]
    n_match_hc = sum(1 for i in range(L) if conservation[i] >= 0.9
                     and i < len(seq) and seq[i] == consensus[i])
    n_hc = sum(1 for c in conservation if c >= 0.9)
    n_cys = sum(1 for i in cys_positions if i < len(seq) and seq[i] == "C")
    n_sub = sum(1 for i in range(min(L, len(seq))) if seq[i] != consensus[i])
    print(f"  SA gen {idx}: pLDDT={plddt_scores[idx]:.1f}, "
          f"conserved={n_match_hc}/{n_hc}, Cys={n_cys}/{len(cys_positions)}, "
          f"substitutions={n_sub}/{L}")

print(f"\nEntropy correlation: r = {r:.3f}")
print(f"Cysteine positions: {[i+1 for i in cys_positions]}")
