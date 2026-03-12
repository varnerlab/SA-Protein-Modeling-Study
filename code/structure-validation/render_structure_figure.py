#!/usr/bin/env python3
"""
Render a multi-panel structure figure for the SA protein paper.
Each panel shows the best SA-generated structure (colored by pLDDT)
superimposed on the aligned region of the reference crystal structure (gray).

Only families with completed AF2 predictions are included.

Usage: conda run python render_structure_figure.py
"""

import json, glob, os, sys
import numpy as np

# ── Configuration ──────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
COLABFOLD_DIR = os.path.join(BASE_DIR, "colabfold_output")
DATA_DIR = os.path.join(BASE_DIR, "data")
FIG_DIR = os.path.join(BASE_DIR, "figs")
PAPER_FIG_DIR = os.path.join(BASE_DIR, "..", "..", "paper", "figs")

FAMILIES = [
    ("PF00076", "RRM",     "1FXL_A"),
    ("PF00018", "SH3",     "1SHG_A"),
    ("PF00397", "WW",      "1PIN_A"),
    ("PF00014", "Kunitz",  "1BPI_A"),
    ("PF00096", "zf-C2H2", "1ZAA_C"),
    ("PF00595", "PDZ",     "1BE9_A"),
    ("PF00069", "Pkinase", "1ATP_E"),
]


def find_best_structure(family_id):
    """Find the SA-generated structure with the highest mean pLDDT."""
    score_dir = os.path.join(COLABFOLD_DIR, family_id, "sa_generation")
    if not os.path.isdir(score_dir):
        return None, None, None

    best_plddt = 0
    best_base = None
    best_scores = None
    for f in glob.glob(os.path.join(score_dir, "*scores*.json")):
        with open(f) as fh:
            d = json.load(fh)
            mean_plddt = np.mean(d["plddt"])
            if mean_plddt > best_plddt:
                best_plddt = mean_plddt
                best_scores = d
                best_base = os.path.basename(f).replace(
                    "_scores_rank_001_alphafold2_ptm_model_1_seed_000.json",
                    "_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb",
                )

    if best_base is None:
        return None, None, None

    pdb_path = os.path.join(score_dir, best_base)
    return pdb_path, best_plddt, best_scores


def render_panel(cmd, sa_pdb, ref_pdb, family_name, plddt, panel_png,
                 width=900, height=750):
    """Render one panel: SA structure (pLDDT-colored) with aligned reference (gray)."""
    cmd.reinitialize()

    # Load structures
    cmd.load(ref_pdb, "reference")
    cmd.load(sa_pdb, "generated")

    # Remove solvent and non-protein
    cmd.remove("solvent")
    cmd.remove("organic")
    cmd.remove("inorganic")

    # Align generated onto reference; get alignment object
    aln = cmd.align("generated", "reference", object="aln_obj")
    rmsd = aln[0]
    n_aligned = aln[1]

    # Remove unaligned residues from reference to avoid the
    # multi-domain size mismatch problem. We keep only reference
    # residues that participate in the alignment.
    cmd.select("ref_aligned", "reference and aln_obj")
    cmd.select("ref_unaligned", "reference and not ref_aligned")
    cmd.remove("ref_unaligned")

    # Delete the alignment object so it doesn't render as dashes
    cmd.delete("aln_obj")
    cmd.deselect()

    # Style: reference in light gray cartoon, slightly transparent
    cmd.show("cartoon", "reference")
    cmd.color("gray80", "reference")
    cmd.set("cartoon_transparency", 0.45, "reference")

    # Style: generated colored by b-factor (pLDDT stored in B-factor column)
    cmd.show("cartoon", "generated")
    cmd.spectrum("b", "red_orange_yellow_green_cyan_blue",
                 "generated", minimum=50, maximum=100)

    # Hide everything else
    cmd.hide("lines")
    cmd.hide("sticks")
    cmd.hide("nonbonded")

    # Aesthetics
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.set("cartoon_oval_length", 1.2)
    cmd.set("cartoon_oval_width", 0.25)
    cmd.set("cartoon_loop_radius", 0.25)
    cmd.set("cartoon_tube_radius", 0.4)
    cmd.set("ray_opaque_background", 1)
    cmd.bg_color("white")
    cmd.set("antialias", 2)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_shadows", 1)
    cmd.set("ray_shadow_decay_factor", 0.1)
    cmd.set("depth_cue", 1)
    cmd.set("fog_start", 0.6)
    cmd.set("spec_reflect", 0.25)
    cmd.set("spec_power", 200)
    cmd.set("ambient", 0.3)
    cmd.set("direct", 0.7)
    cmd.set("light_count", 2)

    # Orient on the generated structure and zoom tight
    cmd.orient("generated")
    cmd.zoom("generated", buffer=2)

    # Render
    cmd.ray(width, height)
    cmd.png(panel_png, dpi=300)

    return rmsd, n_aligned


def compose_figure(panel_pngs, labels, annotations, output_path):
    """Compose individual panels into a multi-panel figure using matplotlib."""
    from matplotlib import pyplot as plt
    from matplotlib import image as mpimg
    import matplotlib.gridspec as gridspec
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    import matplotlib.colors as mcolors

    n = len(panel_pngs)
    if n == 0:
        print("No panels to compose.")
        return

    if n <= 2:
        nrows, ncols = 1, n
    elif n <= 4:
        nrows, ncols = 1, n
    elif n <= 6:
        nrows, ncols = 2, 3
    else:
        nrows, ncols = 2, 4

    fig = plt.figure(figsize=(ncols * 3.2, nrows * 3.5 + 0.4))
    # Leave room at bottom for colorbar
    gs = gridspec.GridSpec(nrows, ncols, wspace=0.08, hspace=0.25,
                           left=0.02, right=0.98, top=0.93, bottom=0.10)

    for i, (png, label, annot) in enumerate(zip(panel_pngs, labels, annotations)):
        row = i // ncols
        col = i % ncols
        ax = fig.add_subplot(gs[row, col])

        if os.path.isfile(png):
            img = mpimg.imread(png)
            ax.imshow(img)
        ax.set_title(label, fontsize=11, fontweight="bold", pad=4)
        # Add annotation (pLDDT, RMSD) below the image
        ax.text(0.5, -0.02, annot, ha="center", va="top",
                transform=ax.transAxes, fontsize=8, color="0.3",
                style="italic")
        ax.axis("off")

    # Hide unused panels
    for i in range(n, nrows * ncols):
        row = i // ncols
        col = i % ncols
        ax = fig.add_subplot(gs[row, col])
        ax.axis("off")

    # Add pLDDT colorbar at the bottom
    cbar_ax = fig.add_axes([0.25, 0.03, 0.50, 0.02])
    # Match the PyMOL spectrum: red-orange-yellow-green-cyan-blue
    cmap_colors = ["#ff0000", "#ff8800", "#ffff00", "#00cc00", "#00cccc", "#0000ff"]
    cmap = mcolors.LinearSegmentedColormap.from_list("plddt", cmap_colors)
    norm = Normalize(vmin=50, vmax=100)
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
    cbar.set_label("pLDDT", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved composite figure: {output_path}")
    plt.close(fig)


def main():
    # Initialize PyMOL in headless mode
    import pymol
    pymol.finish_launching(["pymol", "-cq"])
    from pymol import cmd

    os.makedirs(FIG_DIR, exist_ok=True)
    os.makedirs(PAPER_FIG_DIR, exist_ok=True)

    panel_pngs = []
    labels = []
    annotations = []

    for family_id, family_name, ref_name in FAMILIES:
        ref_pdb = os.path.join(DATA_DIR, family_id, "reference", f"{ref_name}.pdb")
        sa_pdb, plddt, scores = find_best_structure(family_id)

        if sa_pdb is None or not os.path.isfile(sa_pdb):
            print(f"  {family_name} ({family_id}): no AF2 predictions yet, skipping")
            continue

        if not os.path.isfile(ref_pdb):
            print(f"  {family_name} ({family_id}): reference PDB not found, skipping")
            continue

        panel_png = os.path.join(FIG_DIR, f"structure_panel_{family_id}.png")
        print(f"  {family_name} ({family_id}): rendering (pLDDT={plddt:.1f})...")
        rmsd, n_aligned = render_panel(
            cmd, sa_pdb, ref_pdb, family_name, plddt, panel_png
        )
        ptm = scores.get("ptm", 0) if scores else 0
        print(f"    RMSD={rmsd:.2f} A, pTM={ptm:.3f}, {n_aligned} aligned atoms")

        labels.append(family_name)
        annotations.append(f"pLDDT={plddt:.0f}, RMSD={rmsd:.1f} A")
        panel_pngs.append(panel_png)

    if not panel_pngs:
        print("No families with AF2 data. Nothing to render.")
        return

    # Compose the multi-panel figure
    output_path = os.path.join(FIG_DIR, "structure_gallery.pdf")
    compose_figure(panel_pngs, labels, annotations, output_path)

    # Copy to paper figs
    import shutil
    paper_path = os.path.join(PAPER_FIG_DIR, "structure_gallery.pdf")
    shutil.copy2(output_path, paper_path)
    print(f"Copied to: {paper_path}")

    cmd.quit()


if __name__ == "__main__":
    main()
