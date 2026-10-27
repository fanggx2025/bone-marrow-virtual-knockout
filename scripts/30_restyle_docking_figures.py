#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""30_restyle_docking_figures.py — regenerate ALL docking figures in unified
Scientific Reports style (Arial, white bg, 300 dpi, English only, unified palette).

Regenerates (from SAVED results — no re-docking):
  1. Fig_S_wuweizi_PPARG_convergence.png  (Vina bar chart, herb colours = FigS1 scheme)
  2. Fig_S_<Compound>.png x8              (3D pocket maps from *_out.pdbqt + receptor)
  3. Fig4_docking.png                     (Kaempferol SwissDock EADock 3D map)

Unified palette:
  Schisandra chinensis = #2e7d32 (green) | E. senticosus = #7f7f7f (grey)
  receptor grey #cccccc | rosiglitazone (BRL) cyan #17becf
  docked ligand red #d62728 | key-residue hit blue #1f77b4 / miss grey #888888
"""
import os, csv, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJ = r"C:/Users/Fang/Desktop/S3/ciwujia_wuweizi_oa"
OUT  = os.path.join(PROJ, "docking_files/wuweizi_vina")
FIGR = os.path.join(PROJ, "results/figures")
FULL = os.path.join(PROJ, "docking_files/proteins/PPARG_2PRG.pdb")
SRC_K = r"C:/Users/Fang/Downloads/swissdock_results"   # Kaempferol SwissDock EADock results
S10   = os.path.join(PROJ, "results/tables/s10_wuweizi_vina.csv")

# ---- unified style ----
C_SCHI, C_ELE = "#2e7d32", "#7f7f7f"
C_REC, C_BRL, C_LIG = "#cccccc", "#17becf", "#d62728"
C_HIT, C_MISS = "#1f77b4", "#888888"
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "legend.fontsize": 6,
    "axes.unicode_minus": False,
    "savefig.facecolor": "white",
    "figure.facecolor": "white",
})
HERB_C = {"Schisandra": C_SCHI, "Eleutherococcus": C_ELE}
KEY = [("ARG", 288), ("TYR", 473), ("SER", 289), ("HIS", 323), ("HIS", 449), ("GLN", 286)]


def parse_pdb_atoms(path):
    atoms = []
    for line in open(path):
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        atom = line[12:16].strip()
        atoms.append(dict(atom=atom, resn=line[17:20].strip(), ch=line[21:22],
                          resi=int(line[22:26]),
                          x=float(line[30:38]), y=float(line[38:46]), z=float(line[46:54]),
                          el=line[76:78].strip() or atom[0]))
    return atoms


def ligand_atoms_from_pdbqt(path):
    """Read MODEL 1 (best pose) from a vina out.pdbqt."""
    atoms, model, capture = [], None, False
    for line in open(path):
        if line.startswith("MODEL"):
            model = int(line.split()[1]); capture = (model == 1)
        elif line.startswith("ENDMDL"):
            if model == 1:
                break
        elif (line.startswith("ATOM") or line.startswith("HETATM")) and capture:
            atom = line[12:16].strip()
            atoms.append(dict(atom=atom, resn=line[17:20].strip(), resi=int(line[22:26]),
                              x=float(line[30:38]), y=float(line[38:46]), z=float(line[46:54]),
                              el=line[76:78].strip() or atom[0]))
    return atoms


def brl_chain_A_atoms():
    a = []
    for l in open(FULL):
        if l.startswith("HETATM") and l[17:20].strip() == "BRL" and l[21] == "A":
            a.append(np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
    return np.array(a)


def elem_from_name(name):
    m = re.match(r"[A-Za-z*]+", name.strip().lstrip("0123456789"))
    return (m.group(0)[0].upper() if m else "C")


# ---------- load s10 ----------
rows = list(csv.DictReader(open(S10, encoding="utf-8")))
info = {r["compound"]: r for r in rows}

# =====================================================================
# 1. Fig_S_wuweizi_PPARG_convergence.png — Vina bar chart (FigS1a scheme)
# =====================================================================
recs = sorted(rows, key=lambda r: float(r["vina_dG_kcalmol"]))
names = [r["compound"] for r in recs]
dgs   = [float(r["vina_dG_kcalmol"]) for r in recs]
cols  = [HERB_C["Schisandra"] if r["herb"].startswith("Schisandra") else HERB_C["Eleutherococcus"]
         for r in recs]
dbrl  = [float(r["min_dist_to_BRL_A"]) for r in recs]

fig, ax = plt.subplots(figsize=(3.35, 2.5))
y = np.arange(len(names))
ax.barh(y, dgs, color=cols, edgecolor="black", linewidth=0.4)
for yi, v in zip(y, dgs):
    ax.text(v - 0.08, yi, f"{v:.2f}", va="center", ha="right", fontsize=5.5)
ax.set_yticks(y); ax.set_yticklabels(names)
ax.invert_yaxis()
ax.set_xlabel("AutoDock Vina $\\Delta$G (kcal/mol)")
ax.set_xlim(min(dgs) - 1.1, -5.8)
ax.axvline(-7.0, color="grey", ls=":", lw=0.7)
import matplotlib.patches as mpatches
ax.legend(handles=[mpatches.Patch(color=C_SCHI, label="S. chinensis"),
                   mpatches.Patch(color=C_ELE, label="E. senticosus")],
          loc="lower left", frameon=False)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
fig.tight_layout()
p1 = os.path.join(FIGR, "Fig_S_wuweizi_PPARG_convergence.png")
fig.savefig(p1, dpi=300); plt.close(fig)

# vector PDF of the same horizontal bar chart (matplotlib 3D cannot be true vector)
fig2, ax2 = plt.subplots(figsize=(6.0, 3.2))
y = np.arange(len(names))
ax2.barh(y, dgs, color=cols, edgecolor="black", linewidth=0.4)
for yi, v in zip(y, dgs):
    ax2.text(v - 0.08, yi, f"{v:.2f}", va="center", ha="right", fontsize=5.5)
ax2.set_yticks(y); ax2.set_yticklabels(names)
ax2.invert_yaxis()
ax2.set_xlabel("AutoDock Vina $\\Delta$G (kcal/mol)")
ax2.set_xlim(min(dgs) - 1.1, -5.8)
ax2.axvline(-7.0, color="grey", ls=":", lw=0.7)
import matplotlib.patches as mpatches
ax2.legend(handles=[mpatches.Patch(color=C_SCHI, label="S. chinensis"),
                    mpatches.Patch(color=C_ELE, label="E. senticosus")],
           loc="lower left", frameon=False)
for sp in ("top", "right"): ax2.spines[sp].set_visible(False)
fig2.tight_layout()
p1p = os.path.join(FIGR, "Fig_S_wuweizi_PPARG_convergence.pdf")
fig2.savefig(p1p); plt.close(fig2)
print("wrote", p1, "/", p1p)

# =====================================================================
# 2. 8x Fig_S_<Compound>.png — 3D pocket maps (Vina best pose vs chain-A BRL)
# =====================================================================
rec_pdbqt = os.path.join(OUT, "receptor", "PPARG_2PRG_apo_chainA.pdbqt")
rec = parse_pdb_atoms(rec_pdbqt)
rxyz = np.array([[a["x"], a["y"], a["z"]] for a in rec])
BRL = brl_chain_A_atoms()
recby = {}
for a in rec:
    recby.setdefault((a["resn"], a["resi"]), []).append((a["x"], a["y"], a["z"]))

for r in rows:
    name = r["compound"]
    safe = name.replace(" ", "_")
    pq = os.path.join(OUT, "results", safe + "_out.pdbqt")
    if not os.path.exists(pq):
        print("  !! missing", pq); continue
    lig = ligand_atoms_from_pdbqt(pq)
    if not lig:
        print("  !! empty pose", name); continue
    lxyz = np.array([[a["x"], a["y"], a["z"]] for a in lig])
    herb = r["herb"]
    herb_s = "S. chinensis" if herb.startswith("Schisandra") else "E. senticosus"

    fig = plt.figure(figsize=(3.35, 2.9))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(rxyz[::3, 0], rxyz[::3, 1], rxyz[::3, 2], s=0.8, c=C_REC, alpha=0.30,
               label="PPARG LBD (chain A, apo)")
    ax.scatter(BRL[:, 0], BRL[:, 1], BRL[:, 2], s=14, c=C_BRL, marker="^",
               depthshade=False, label="Rosiglitazone (BRL)")
    ax.scatter(lxyz[:, 0], lxyz[:, 1], lxyz[:, 2], s=16, c=C_LIG, marker="o",
               depthshade=False, label="docked pose")
    keypts = []
    for (rn, ri) in KEY:
        pts = recby.get((rn, ri))
        if not pts:
            continue
        c = np.array(pts).mean(axis=0)
        keypts.append(c)
        d = min(np.linalg.norm(lxyz - np.array(p), axis=1).min() for p in pts)
        col = C_HIT if d <= 4.5 else C_MISS
        ax.scatter([c[0]], [c[1]], [c[2]], s=26, c=col, marker="*", depthshade=False)
        ax.text(c[0], c[1], c[2], " %s%d" % (rn, ri), color=col, fontsize=5.5)
    allp = np.vstack([lxyz, BRL] + ([np.array(keypts)] if keypts else []))
    cen = allp.mean(axis=0); rad = float(np.abs(allp - cen).max()) + 6
    ax.set_xlim(cen[0] - rad, cen[0] + rad)
    ax.set_ylim(cen[1] - rad, cen[1] + rad)
    ax.set_zlim(cen[2] - rad, cen[2] + rad)
    ax.set_xlabel("X (Å)", labelpad=-6); ax.set_ylabel("Y (Å)", labelpad=-6)
    ax.set_zlabel("Z (Å)", labelpad=-6)
    ax.tick_params(pad=-2)
    ax.set_title("%s (%s)\n$\\Delta$G = %.2f kcal/mol, %.1f Å to BRL"
                 % (name, herb_s, float(r["vina_dG_kcalmol"]), float(r["min_dist_to_BRL_A"])),
                 fontsize=7)
    ax.legend(loc="upper left", fontsize=5, frameon=False)
    fig.tight_layout()
    pf = os.path.join(FIGR, "Fig_S_%s.png" % safe)
    fig.savefig(pf, dpi=300)
    fig.savefig(pf[:-4] + ".pdf")
    plt.close(fig)
    print("wrote", pf, "+ pdf")

# =====================================================================
# 3. Fig4_docking.png — Kaempferol SwissDock EADock (best pose 3D map)
# =====================================================================
def parse_pdb_generic(path, het_only=False, exclude_res=None):
    out = []
    for ln in open(path):
        if not (ln.startswith("ATOM") or (het_only and ln.startswith("HETATM"))):
            if not het_only and not ln.startswith("ATOM"):
                continue
        if exclude_res and ln[17:20].strip() in exclude_res:
            continue
        xyz = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        out.append(dict(atm=ln[12:16].strip(), res=ln[17:20].strip(), resi=int(ln[22:26]),
                        xyz=xyz, el=ln[76:78].strip() or elem_from_name(ln[12:16])))
    return out


prot = parse_pdb_generic(os.path.join(SRC_K, "files/receptor.pdb"), exclude_res={"BRL", "HOH"})
lines = open(os.path.join(SRC_K, "result.dock4")).read().splitlines()
lig, in_block = [], False
for ln in lines:
    if ln.startswith("REMARK CLUSTER_NUM"):
        in_block = (int(ln.split(":")[1]) == 0); continue
    if in_block and ln.startswith("REMARK CLUSTER_MEMBER"):
        in_block = (int(ln.split(":")[1]) == 1); continue
    if in_block and ln.startswith("ATOM"):
        lig.append(dict(atm=ln[12:16].strip(),
                        xyz=np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]),
                        el=elem_from_name(ln[12:16])))
print("kaempferol EADock best-pose atoms:", len(lig))

pr = np.array([p["xyz"] for p in prot])
lg = np.array([a["xyz"] for a in lig])
mind = min(np.linalg.norm(a["xyz"] - p["xyz"]) for a in lig for p in prot)

fig = plt.figure(figsize=(3.35, 2.9))
ax = fig.add_subplot(111, projection="3d")
ax.scatter(pr[:, 0], pr[:, 1], pr[:, 2], s=1.2, c=C_REC, alpha=0.35, label="PPARG LBD (receptor)")
ax.scatter(lg[:, 0], lg[:, 1], lg[:, 2], s=14, c=C_LIG, marker="o", depthshade=False,
           label="Kaempferol (EADock best)")
KEY4 = ["ARG288", "SER289", "TYR473", "HIS449", "MET364", "LEU453", "ILE341", "CYS285", "PHE282"]
for key in KEY4:
    rn, ri = re.match(r"([A-Z]+)(\d+)", key).groups(); ri = int(ri)
    sub = [p for p in prot if p["res"] == rn and p["resi"] == ri and p["atm"] == "CA"] or \
          [p for p in prot if p["res"] == rn and p["resi"] == ri]
    if not sub:
        continue
    dd = min(np.linalg.norm(a["xyz"] - p["xyz"]) for a in lig for p in sub)
    c = sub[0]["xyz"]
    hit = dd <= 4.0
    ax.scatter([c[0]], [c[1]], [c[2]], s=34, c=C_HIT if hit else C_MISS,
               marker="*" if hit else "o", depthshade=False, edgecolors="k", linewidths=0.3)
    ax.text(c[0], c[1], c[2], " %s" % key, color=C_HIT if hit else C_MISS, fontsize=5.5)
pad = 13
ax.set_xlim(lg[:, 0].min() - pad, lg[:, 0].max() + pad)
ax.set_ylim(lg[:, 1].min() - pad, lg[:, 1].max() + pad)
ax.set_zlim(lg[:, 2].min() - pad, lg[:, 2].max() + pad)
ax.set_xlabel("X (Å)", labelpad=-6); ax.set_ylabel("Y (Å)", labelpad=-6)
ax.set_zlabel("Z (Å)", labelpad=-6)
ax.tick_params(pad=-2)
ax.set_title("Kaempferol (E. senticosus) — SwissDock EADock\n"
             "$\\Delta$G = −6.76 kcal/mol, min lig–prot = %.2f Å" % mind, fontsize=7)
ax.legend(loc="upper left", fontsize=5.5, frameon=False)
fig.tight_layout()
fp = os.path.join(FIGR, "Fig4_docking.png")
fig.savefig(fp, dpi=300)
fig.savefig(fp[:-4] + ".pdf")
plt.close(fig)
print("wrote", fp, "+ pdf")
print("ALL DONE")
