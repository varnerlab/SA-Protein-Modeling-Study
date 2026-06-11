#!/usr/bin/env python3
"""
Re-compose the multi-panel structure figure (v2) from the already-rendered
PyMOL panel PNGs. No PyMOL required.

Differences from v1:
  1) Larger protein diagrams  -> each PNG is auto-cropped to its non-white
     content so the structure fills the panel, and the per-panel footprint
     is enlarged.
  2) Larger score annotations -> annotation fontsize bumped.
  3) Less white space between structures -> tighter wspace/hspace and the
     auto-crop removes PyMOL's internal padding.
  4) Larger structure labels  -> title fontsize bumped.

Usage: python render_structure_figure_v2.py
"""

import os
import numpy as np

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(BASE_DIR, "figs")
PAPER_FIG_DIR = os.path.join(BASE_DIR, "..", "..", "paper", "figs")

# (family_id, label, "pLDDT=.., RMSD=.. A") -- annotation values taken
# verbatim from the v1 figure so the science is unchanged.
PANELS = [
    ("PF00076", "RRM",           "pLDDT=91, RMSD=0.6 Å"),
    ("PF00018", "SH3",           "pLDDT=96, RMSD=0.7 Å"),
    ("PF00397", "WW",            "pLDDT=91, RMSD=0.9 Å"),
    ("PF00014", "Kunitz",        "pLDDT=97, RMSD=0.6 Å"),
    ("PF00096", "zf-C2H2",       "pLDDT=96, RMSD=0.7 Å"),
    ("PF00595", "PDZ",           "pLDDT=94, RMSD=0.4 Å"),
    ("PF00069", "Pkinase",       "pLDDT=92, RMSD=1.5 Å"),
    ("PF00711", "Defensin_beta", "pLDDT=93, RMSD=1.8 Å"),
]


def autocrop(img, pad_frac=0.02, thresh=0.985):
    """Crop a white-background RGB(A) image to its content bounding box,
    leaving a small uniform margin (fraction of the cropped size)."""
    rgb = img[..., :3] if img.shape[-1] == 4 else img
    # non-white = any channel darker than threshold
    mask = np.any(rgb < thresh, axis=-1)
    if not mask.any():
        return img
    rows = np.any(mask, axis=1)
    cols = np.any(mask, axis=0)
    r0, r1 = np.where(rows)[0][[0, -1]]
    c0, c1 = np.where(cols)[0][[0, -1]]
    h, w = r1 - r0 + 1, c1 - c0 + 1
    pr, pc = int(h * pad_frac), int(w * pad_frac)
    r0 = max(0, r0 - pr); r1 = min(img.shape[0] - 1, r1 + pr)
    c0 = max(0, c0 - pc); c1 = min(img.shape[1] - 1, c1 + pc)
    return img[r0:r1 + 1, c0:c1 + 1]


def compose(output_path):
    from matplotlib import pyplot as plt
    from matplotlib import image as mpimg
    import matplotlib.gridspec as gridspec
    from matplotlib.colors import Normalize, LinearSegmentedColormap
    from matplotlib.cm import ScalarMappable

    panels = [(fid, lab, ann) for (fid, lab, ann) in PANELS
              if os.path.isfile(os.path.join(FIG_DIR, f"structure_panel_{fid}.png"))]
    n = len(panels)
    if n == 0:
        print("No panel PNGs found.")
        return

    nrows, ncols = (2, 4) if n > 4 else (1, n)

    # Bigger per-panel footprint than v1 (was 3.2 x 3.5).
    fig = plt.figure(figsize=(ncols * 3.7, nrows * 3.9 + 0.4))
    gs = gridspec.GridSpec(
        nrows, ncols,
        wspace=0.01, hspace=0.16,           # tighter than v1 (0.08 / 0.25)
        left=0.01, right=0.99, top=0.94, bottom=0.11,
    )

    for i, (fid, label, annot) in enumerate(panels):
        ax = fig.add_subplot(gs[i // ncols, i % ncols])
        img = mpimg.imread(os.path.join(FIG_DIR, f"structure_panel_{fid}.png"))
        ax.imshow(autocrop(img))
        ax.set_title(label, fontsize=16, fontweight="bold", pad=5)   # v1: 11
        ax.text(0.5, -0.015, annot, ha="center", va="top",
                transform=ax.transAxes, fontsize=12, color="0.3",    # v1: 8
                style="italic")
        ax.axis("off")

    for i in range(n, nrows * ncols):
        ax = fig.add_subplot(gs[i // ncols, i % ncols])
        ax.axis("off")

    # pLDDT colorbar at the bottom (matches PyMOL spectrum)
    cbar_ax = fig.add_axes([0.30, 0.035, 0.40, 0.022])
    cmap = LinearSegmentedColormap.from_list(
        "plddt",
        ["#ff0000", "#ff8800", "#ffff00", "#00cc00", "#00cccc", "#0000ff"])
    sm = ScalarMappable(norm=Normalize(vmin=50, vmax=100), cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
    cbar.set_label("pLDDT", fontsize=12)
    cbar.ax.tick_params(labelsize=10)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {output_path}")


def main():
    out = os.path.join(FIG_DIR, "structure_gallery_v2.pdf")
    compose(out)
    import shutil
    paper = os.path.join(PAPER_FIG_DIR, "structure_gallery_v2.pdf")
    shutil.copy2(out, paper)
    print(f"Copied to: {paper}")


if __name__ == "__main__":
    main()
