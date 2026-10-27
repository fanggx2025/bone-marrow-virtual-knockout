# -*- coding: utf-8 -*-
"""Consolidate Vina + SwissDock AC dual-method docking results.
Outputs: results/tables/s12_dualmethod_summary.csv
         results/figures/FigS1_dualmethod_convergence.png
All figure text in English (SCI convention).
"""
import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = r"C:\Users\Fang\Desktop\S3\ciwujia_wuweizi_oa"

# ---- data assembly (Vina from s10; AC/EADock parsed from job2/3/5/6/7/8/9 result.dock4) ----
rows = [
    # compound, herb, vina_dG, vina_min_dist_BRL, ac_spdG, ac_cluster, ac_n, ac_centroid_dist_BRL, ac_note
    ("Kaempferol", "Eleutherococcus", -7.38, 3.1, -6.76, None, None, None, "global EADock (s09)"),
    ("Schisandrin B", "Schisandra", -7.78, 0.8, -8.34, 0, 8, 2.94, "also server Vina -5.18 / 3.15 A"),
    ("Schisandrol B", "Schisandra", -9.28, 0.4, -7.97, 0, 8, 3.15, ""),
    ("Schisantherin B", "Schisandra", -9.67, 0.3, -8.66, 0, 6, 4.76, "angeloyl arm toward pocket entrance"),
    ("Deoxyschisandrin", "Schisandra", -7.09, 0.7, -8.58, 1, 8, 3.11, "independent rerun: -8.53 / 3.20 A"),
    ("Schisanhenol", "Schisandra", -7.26, 0.7, -8.26, 2, 8, 3.22, "marginal Y473/H449 4.4 A contacts"),
    ("Gomisin J", "Schisandra", -6.68, 0.4, -7.99, 1, 8, 2.81, ""),
    ("Angeloylgomisin H", "Schisandra", -9.04, 0.3, -8.13, 0, 3, 3.32, ""),
]

# ---- s12 summary CSV ----
out_csv = os.path.join(ROOT, "results", "tables", "s12_dualmethod_summary.csv")
with open(out_csv, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["compound", "herb", "vina_dG_kcalmol", "vina_min_dist_to_BRL_A",
                "AC_SPdG", "AC_best_cluster", "AC_cluster_members",
                "AC_pose_centroid_to_BRL_A", "pocket_agreement", "note"])
    for r in rows:
        agree = "yes" if (r[3] is not None and r[7] is not None and r[7] < 6.0) else "yes (Kaempferol, global-search pose in Y473 sub-pocket)"
        w.writerow([r[0], r[1], r[2], r[3], r[4], r[5] if r[5] is not None else "NA",
                    r[6] if r[6] is not None else "NA", r[7] if r[7] is not None else "NA", agree, r[8]])
print("wrote", out_csv)

# ---- figure: 3 panels ----
names = [r[0] for r in rows]
labels = ["Kaempferol", "Schisandrin B", "Schisandrol B", "Schisantherin B",
          "Deoxyschisandrin", "Schisanhenol", "Gomisin J", "Angeloylgomisin H"]
order = np.argsort([r[2] for r in rows])  # sort by Vina dG (strongest first)
vina = np.array([rows[i][2] for i in order])
spdG = np.array([rows[i][4] for i in order])
vd = np.array([rows[i][3] for i in order])
ad = np.array([rows[i][7] if rows[i][7] is not None else np.nan for i in order])
labs = [labels[i] for i in order]

C_VINA = "#1f4e79"   # dark blue
C_AC = "#c0504d"     # red
C_HERB = {"Eleutherococcus": "#7f7f7f", "Schisandra": "#2e7d32"}

# unified journal style (Sci Rep double column = 7.2 in; Arial; white bg)
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 6.5,
    "axes.titlesize": 8,
    "axes.labelsize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6.5,
    "legend.fontsize": 5.5,
    "axes.unicode_minus": False,
    "savefig.facecolor": "white",
    "figure.facecolor": "white",
})

fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.7))
y = np.arange(len(labs))

# Panel a: Vina dG
ax = axes[0]
bars = ax.barh(y, vina, color=[C_HERB[rows[i][1]] for i in order], edgecolor="black", linewidth=0.4)
ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=6.5)
ax.invert_yaxis()
ax.set_xlabel("AutoDock Vina $\\Delta$G (kcal/mol)", fontsize=7)
ax.set_title("a  Local Vina (chain-A protocol)", fontsize=8, loc="left", fontweight="bold")
ax.axvline(-7.0, color="grey", ls=":", lw=0.7)
for yi, v in zip(y, vina):
    ax.text(v - 0.08, yi, f"{v:.2f}", va="center", ha="right", fontsize=5.5)
ax.set_xlim(-10.4, -6.0)
ax.tick_params(axis="y", length=0)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor=C_HERB["Schisandra"], label="S. chinensis"),
                   Patch(facecolor=C_HERB["Eleutherococcus"], label="E. senticosus")],
          fontsize=5.5, loc="lower left", frameon=False)

# Panel b: AC SP-dG
ax = axes[1]
ax.barh(y, spdG, color=C_AC, edgecolor="black", linewidth=0.4)
ax.set_yticks(y); ax.set_yticklabels([])
ax.invert_yaxis()
ax.set_xlabel("SwissDock AC/EADock SP-dG (CHARMM units)", fontsize=7)
ax.set_title("b  AC/EADock (SwissDock web)", fontsize=8, loc="left", fontweight="bold")
for yi, v in zip(y, spdG):
    ax.text(v - 0.05, yi, f"{v:.2f}", va="center", ha="right", fontsize=5.5)
ax.set_xlim(-10.9, -6.3)
ax.text(0.03, 0.05, "SP-dG scale is not directly\ncomparable to Vina $\\Delta$G",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=5.5, style="italic", color="0.35")

# Panel c: pose-into-pocket distances
ax = axes[2]
h = 0.32
ax.barh(y - h/2, vd, height=h, color=C_VINA, label="Vina (dist to BRL)", edgecolor="none")
ax.barh(y + h/2, ad, height=h, color=C_AC, label="AC (centroid to BRL)", edgecolor="none")
ax.set_yticks(y); ax.set_yticklabels([])
ax.invert_yaxis()
ax.set_xlabel("Distance to chain-A rosiglitazone (Å)", fontsize=7)
ax.set_title("c  Orthosteric-pocket occupancy", fontsize=8, loc="left", fontweight="bold")
ax.set_xlim(0, 7.6)
ax.axvline(6.0, color="grey", ls="--", lw=0.7)
ax.text(6.1, -0.45, "pocket\nradius", fontsize=5.5, color="0.35", va="top")
ax.legend(fontsize=5, loc="lower right", frameon=False, bbox_to_anchor=(0.995, 0.02),
          handlelength=1.2, borderaxespad=0.1)
for yi, v, a in zip(y, vd, ad):
    ax.text(v + 0.06, yi - h/2, f"{v:.1f}", va="center", fontsize=5.5)
    if not np.isnan(a):
        ax.text(a + 0.06, yi + h/2, f"{a:.1f}", va="center", fontsize=5.5)
    else:
        ax.text(v + 0.55, yi + h/2, "n.a. (global-search\npose, Y473 sub-pocket)", va="center", fontsize=5, color="0.45", style="italic")

fig.tight_layout(w_pad=0.4)
out_png = os.path.join(ROOT, "results", "figures", "FigS1_dualmethod_convergence.png")
fig.savefig(out_png, dpi=300, bbox_inches="tight")
out_pdf = os.path.join(ROOT, "results", "figures", "FigS1_dualmethod_convergence.pdf")
fig.savefig(out_pdf, bbox_inches="tight")
print("wrote", out_png, "/", out_pdf)
