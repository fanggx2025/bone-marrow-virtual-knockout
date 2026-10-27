#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
26_wuweizi_docking.py  (v2 — chain-A corrected protocol)
=========================================================
Dock the 7 Schisandra chinensis (五味子) lignans — already assigned PPARG+AKT1
in the literature-target layer (s05b) — into the PPARG ligand-binding domain
(PDB 2PRG) with AutoDock Vina 1.2.7, together with Kaempferol (刺五加), to
quantify cross-herb convergence on PPARG at the direct-docking level.

v2 CORRECTION (important):
  v1 averaged the rosiglitazone (BRL) centroids of BOTH chains of the 2PRG
  homodimer (pockets 40.3 A apart) -> box center sat 20.2 A from EACH pocket,
  i.e. inside no orthosteric pocket at all. All v1 poses were between-monomer
  artifacts (archived under archive_2chain_midbox/). v2:
    - receptor rebuilt from chain A ATOM records only
    - box centered on chain A co-crystallized rosiglitazone (BRL) centroid
    - Kaempferol re-docked under the SAME protocol (no hardcoded reference)
    - BRL overlay in figures = chain A BRL only
    - residue parsing no longer merges identical residue numbers across chains

Pipeline:
  1. Read authoritative ChEMBL SMILES for the 7 lignans from s05_compounds.csv
  2. Ligand 3D PDBQT via NCI Cactus 3D SDF -> openbabel
  3. Receptor PDBQT from chain-A-only apo PDB (BRL + waters + chain B/C stripped)
  4. Vina 1.2.7, 20 A box centered on chain-A BRL centroid, exhaustiveness 8, seed 42
  5. Parse best pose, contacts to canonical chain-A pocket residues, distance to BRL
  6. Per-compound figure + cross-herb convergence figure + s10_wuweizi_vina.csv

Run with the managed python venv (openbabel-wheel + numpy installed).
"""
import os, csv, subprocess, sys, glob, re, ssl, shutil
import urllib.request, urllib.parse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa

# ---- paths ----
PROJ = r"C:/Users/Fang/Desktop/S3/ciwujia_wuweizi_oa"
S05  = os.path.join(PROJ, "results/tables/s05_compounds.csv")
FULL = os.path.join(PROJ, "docking_files/proteins/PPARG_2PRG.pdb")   # full 2PRG (has BRL in chains A+B)
APO_A= os.path.join(PROJ, "docking_files/PPARG_2PRG_apo_chainA.pdb") # chain-A-only apo (generated here)
VINA = os.path.join(PROJ, "docking_files/bin/vina_1.2.7_win.exe")
OUT  = os.path.join(PROJ, "docking_files/wuweizi_vina")
for d in ("ligands","receptor","results","figures","complexes"):
    os.makedirs(os.path.join(OUT,d), exist_ok=True)
os.makedirs(os.path.join(PROJ,"results/figures"), exist_ok=True)
os.makedirs(os.path.join(PROJ,"results/tables"), exist_ok=True)

# openbabel data dir (silences MMFF/fragment warnings)
os.environ["BABEL_DATADIR"] = r"C:/Users/Fang/.workbuddy/binaries/python/envs/default/Lib/site-packages/openbabel/bin/data"

# canonical PPARG pocket residues (2PRG numbering, chain A)
POCKET = [("ARG",288),("SER",289),("HIS",323),("HIS",449),("TYR",473),
          ("GLN",286),("CYS",285),("PHE",282),("ILE",281),("GLY",284),
          ("LEU",330),("TYR",327),("ILE",326),("VAL",339),("ILE",341),
          ("LEU",340),("MET",348),("MET",364),("PHE",363),("LEU",453),("LEU",469)]
KEY = [("ARG",288),("TYR",473),("SER",289),("HIS",323),("HIS",449),("GLN",286)]
KEY_S = ["%s%d"%k for k in KEY]   # string keys matching key_hit dict ("ARG288", ...)

# ---------- helpers ----------
def parse_pdb_atoms(path):
    atoms=[]
    for line in open(path):
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        rec=line[0:6].strip()
        atom=line[12:16].strip()
        resn=line[17:20].strip()
        ch=line[21:22]
        resi=int(line[22:26])
        x=float(line[30:38]); y=float(line[38:46]); z=float(line[46:54])
        el=line[76:78].strip() or atom[0]
        atoms.append(dict(rec=rec,atom=atom,resn=resn,ch=ch,resi=resi,x=x,y=y,z=z,el=el))
    return atoms

def ligand_atoms_from_pdbqt(path):
    """Read MODEL 1 (best) from vina pdbqt output."""
    atoms=[]; model=None; capture=False
    for line in open(path):
        if line.startswith("MODEL"):
            model=int(line.split()[1]); capture=(model==1)
        elif line.startswith("ENDMDL"):
            if model==1: break
        elif (line.startswith("ATOM") or line.startswith("HETATM")) and capture:
            atom=line[12:16].strip()
            resn=line[17:20].strip(); resi=int(line[22:26])
            x=float(line[30:38]); y=float(line[38:46]); z=float(line[46:54])
            el=line[76:78].strip() or atom[0]
            atoms.append(dict(atom=atom,resn=resn,resi=resi,x=x,y=y,z=z,el=el))
    return atoms

def brl_chain_A_atoms():
    a=[]
    for l in open(FULL):
        if l.startswith("HETATM") and l[17:20].strip()=="BRL" and l[21]=="A":
            a.append(np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])]))
    return np.array(a)

CTX = ssl.create_default_context(); CTX.check_hostname=False; CTX.verify_mode=ssl.CERT_NONE
def get_3d_sdf(query):
    """Fetch a real 3D SDF from NCI Cactus (robust for natural products; the
    local openbabel OBBuilder fails on complex lignan ring systems)."""
    u="https://cactus.nci.nih.gov/chemical/structure/%s/file?format=sdf&get3d=True"%urllib.parse.quote(query)
    req=urllib.request.Request(u, headers={"User-Agent":"Mozilla/5.0"})
    return urllib.request.urlopen(req, timeout=60, context=CTX).read()

# ---------- 0. chain-A BRL centroid = docking box center ----------
BRL_ATOMS = brl_chain_A_atoms()
CENTER = BRL_ATOMS.mean(axis=0)
print("chain-A BRL atoms: %d ; centroid (box center): %s"%(len(BRL_ATOMS), CENTER.round(2)))

# ---------- 1. compound list ----------
# 7 lignans from s05_compounds.csv (authoritative ChEMBL SMILES) + Kaempferol.
# Kaempferol SMILES is not in s05 (it is the Eleutherococcus reference compound);
# fetch by NAME via Cactus, with fallback to the pre-built ligand pdbqt.
LIGNANS = ["Schisandrin B","Schisantherin B","Deoxyschisandrin",
           "Schisanhenol","Schisandrol B","Gomisin J","Angeloylgomisin H"]
LIGS = []  # (name, query-for-cactus, cn_name, chembl_id, herb)
smiles_map={}; cn_map={}
with open(S05, encoding="utf-8") as f:
    for r in csv.DictReader(f):
        if r["query"] in LIGNANS:
            smiles_map[r["query"]] = r["smiles"]
            cn_map[r["query"]]=(r.get("cn_name",""), r.get("chembl_id",""))
missing=[n for n in LIGNANS if n not in smiles_map]
if missing:
    print("FATAL: missing SMILES for %s"%missing); sys.exit(1)
LIGS=[(n, smiles_map[n], cn_map[n][0], cn_map[n][1], "Schisandra") for n in LIGNANS]
LIGS.insert(0, ("Kaempferol","kaempferol","Kaempferol","CHEMBL150","Eleutherococcus"))
print("Compounds to dock: %d (1 Eleutherococcus + 7 Schisandra)"%len(LIGS))

# ---------- 2. chain-A-only apo receptor PDBQT ----------
from openbabel import pybel
if not os.path.exists(APO_A):
    with open(APO_A,"w") as fh:
        for l in open(FULL):
            if l.startswith("ATOM") and l[21]=="A":   # protein chain A only; drops BRL/HOH/chainB/C
                fh.write(l)
        fh.write("END\n")
    print("wrote chain-A apo PDB:", APO_A)
rec_pdbqt = os.path.join(OUT,"receptor","PPARG_2PRG_apo_chainA.pdbqt")
if not os.path.exists(rec_pdbqt):
    rpm = pybel.readfile("pdb", APO_A).__next__()
    rpm.addh()
    rpm.write("pdbqt", rec_pdbqt, overwrite=True)
    # Vina rigid receptor must NOT contain ROOT/ENDROOT/BRANCH/TORSDOF torsion-tree tags
    lines = open(rec_pdbqt).read().splitlines()
    keep = [l for l in lines if not any(t in l for t in ("ROOT","ENDROOT","BRANCH","TORSDOF"))]
    open(rec_pdbqt, "w").write("\n".join(keep) + "\n")
    print("wrote receptor", rec_pdbqt, "(chain A only, torsion-tree stripped)")

# ---------- 3+4+5. per-compound docking ----------
# docking box: 20 A cube centered on chain-A co-crystallized rosiglitazone.
# The chain-B pocket is 40.3 A away, so no cross-monomer leakage is possible.
BOX = 20.0

SUM=[]  # (name, cn, cid, best, n_contacts, n_key, d2brl, keymsg, herb)
for (name, query, cn, cid, herb) in LIGS:
    safe = name.replace(" ","_")
    lig_pdbqt = os.path.join(OUT,"ligands", safe+".pdbqt")
    sdf_path  = os.path.join(OUT,"ligands", safe+".sdf")
    # build ligand 3D pdbqt (NCI Cactus 3D SDF -> openbabel), with Kaempferol fallback
    if not os.path.exists(lig_pdbqt):
        if not os.path.exists(sdf_path):
            try:
                open(sdf_path,"wb").write(get_3d_sdf(query))
            except Exception as e:
                alt=os.path.join(PROJ,"docking_files/ligands", safe+".pdbqt")
                if os.path.exists(alt):
                    shutil.copy(alt, lig_pdbqt)
                    print("  %s: cactus SDF failed (%s); reusing pre-built pdbqt"%(name,str(e)[:40]))
                else:
                    print("  %s: FAILED to build ligand (%s) — skipping"%(name,str(e)[:60]))
                    SUM.append((name,cn,cid,None,0,0,float('nan'),"",herb)); continue
        if not os.path.exists(lig_pdbqt):
            mol = pybel.readfile("sdf", sdf_path).__next__()
            mol.addh()
            mol.write("pdbqt", lig_pdbqt, overwrite=True)
    out_pdbqt = os.path.join(OUT,"results", safe+"_out.pdbqt")
    log_txt   = os.path.join(OUT,"results", safe+"_log.txt")
    # reuse cached docking when present (Vina is deterministic: same box, seed 42,
    # exhaustiveness 8 -> identical poses on re-run, so re-docking is wasteful)
    if os.path.exists(out_pdbqt) and os.path.exists(log_txt) and os.path.getsize(log_txt)>0:
        with open(log_txt) as fh: res_stdout = fh.read()
        print("  (cached) %s — reusing Vina output"%name, flush=True)
    else:
        print("  docking %s ..."%name, flush=True)
        cmd=[VINA,"--receptor",rec_pdbqt,"--ligand",lig_pdbqt,
             "--center_x","%.3f"%CENTER[0],"--center_y","%.3f"%CENTER[1],"--center_z","%.3f"%CENTER[2],
             "--size_x","%.1f"%BOX,"--size_y","%.1f"%BOX,"--size_z","%.1f"%BOX,
             "--num_modes","9","--exhaustiveness","8","--seed","42",
             "--out",out_pdbqt]
        try:
            res=subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=240)
        except subprocess.TimeoutExpired:
            print("  %-18s TIMEOUT (>240s) — skipping"%name); SUM.append((name,cn,cid,None,0,0,float('nan'),"",herb)); continue
        except subprocess.CalledProcessError as e:
            print("  %-18s VINA ERROR — skipping (%s)"%(name,str(e)[:60])); SUM.append((name,cn,cid,None,0,0,float('nan'),"",herb)); continue
        res_stdout = res.stdout
        with open(log_txt,"w") as fh: fh.write(res_stdout)
    # parse best score from stdout (Vina 1.2.7 results table)
    best=None
    for line in res_stdout.splitlines():
        m=re.match(r"\s*1\s+([-\d.]+)", line)
        if m:
            best=float(m.group(1)); break
    if best is None or not os.path.exists(out_pdbqt):
        print("  %-18s docking produced NO poses; skipping contacts/figure"%name)
        SUM.append((name,cn,cid,None,0,0,float('nan'),"",herb)); continue
    lig = ligand_atoms_from_pdbqt(out_pdbqt)
    rec = parse_pdb_atoms(rec_pdbqt)
    # distance from best pose to nearest chain-A rosiglitazone (BRL) atom
    # -> is the pose actually inside the orthosteric pocket?
    if len(BRL_ATOMS) and len(lig):
        d2brl = float(min(np.linalg.norm(BRL_ATOMS - np.array([a['x'],a['y'],a['z']]), axis=1).min()
                          for a in lig))
    else:
        d2brl = float('nan')
    # contacts (chain-A residues only — no cross-chain residue-number merging)
    ligxyz=np.array([[a['x'],a['y'],a['z']] for a in lig])
    recby={}
    for a in rec:
        recby.setdefault((a['resn'],a['resi']),[]).append([a['x'],a['y'],a['z'],a['el']])
    contacts=[]
    key_hit={}
    for (rn,ri) in POCKET:
        pts=[p for p in recby.get((rn,ri),[])]
        if not pts: continue
        mind=min(np.linalg.norm(ligxyz-np.array(p[:3]),axis=1).min() for p in pts)
        key_hit["%s%d"%(rn,ri)]=mind
        if mind<=4.5:
            contacts.append((rn,ri,round(mind,2)))
    # write complex pdb for pymol
    cplx=os.path.join(OUT,"complexes",safe+"_complex.pdb")
    with open(cplx,"w") as fh:
        fh.write(open(rec_pdbqt).read())
        for i,a in enumerate(lig,1):
            fh.write("HETATM%5d %-4s LIG A   1    %8.3f%8.3f%8.3f  1.00  0.00           %s\n"%(
                i, a['atom'][:4], a['x'],a['y'],a['z'], a['el']))
        fh.write("END\n")
    # contact csv
    ccsv=os.path.join(OUT,"results",safe+"_contacts.csv")
    with open(ccsv,"w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["residue","distance_A"])
        for rn,ri,d in sorted(contacts,key=lambda x:x[2]):
            w.writerow(["%s%d"%(rn,ri), d])
    nkey=sum(1 for k in KEY_S if key_hit.get(k) and key_hit[k]<=4.5)
    SUM.append((name,cn,cid,best,len(contacts),nkey,d2brl,
                ",".join("%s=%.1f"%(k,key_hit[k]) for k in KEY_S if key_hit.get(k)), herb))
    print("  %-18s dG=%.2f  contacts=%d  key_hits=%d  dist2BRL=%.1f"%(name,best,len(contacts),nkey,d2brl))

# ---------- 6a. combined convergence figure (all compounds on ONE Vina scale) ----------
HERB_EN = {"Eleutherococcus":"Eleutherococcus\n(Kaempferol)","Schisandra":"Schisandra\n(lignans)"}
combo = [(s[0], s[3] if s[3] is not None else float('nan'), s[8], s[6]) for s in SUM]
labels=["%s\n%s"%(c[0], HERB_EN.get(c[2],c[2])) for c in combo]
dgs=[c[1] for c in combo]
cols=["#2ca02c" if c[2]=="Eleutherococcus" else "#d62728" for c in combo]
fig=plt.figure(figsize=(11,5.5))
ax=fig.add_subplot(111)
bars=ax.bar(range(len(combo)), dgs, color=cols)
ax.axhline(-6.76, ls=":", c="#2ca02c", lw=1.0, label="Kaempferol SwissDock -6.76 (independent ref)")
ax.set_xticks(range(len(combo))); ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("Vina binding affinity (kcal/mol, more negative = stronger)")
ax.set_title("Cross-herb convergence on PPARG (PDB 2PRG chain A LBD)\n"
             "Vina 1.2.7, 20 A box centered on co-crystallized rosiglitazone (chain A)\n"
             "hatched = pose > 8 A from rosiglitazone (alternative LBD cavity, not the orthosteric pocket)")
for i,b in enumerate(bars):
    if dgs[i]==dgs[i]:  # not nan
        off = (combo[i][3]==combo[i][3]) and (combo[i][3] > 8.0)
        if off:
            b.set_hatch("//")
            b.set_edgecolor("black")
        lbl = "%.2f%s" % (dgs[i], " (off-pocket)" if off else "")
        ax.text(b.get_x()+b.get_width()/2, -0.30, lbl,
                ha="center", va="top", fontsize=7, color="white")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color="#2ca02c",label="Eleutherococcus senticosus"),
                   Patch(color="#d62728",label="Schisandra chinensis")], fontsize=8, loc="lower right")
plt.tight_layout()
cf=os.path.join(OUT,"figures","Fig_S_wuweizi_PPARG_convergence.png")
plt.savefig(cf, dpi=150); plt.close()
print("wrote", cf)

# ---------- 6b. per-compound 3D pocket figures (with chain-A Rosiglitazone overlay) ----------
for (name,cn,cid,best,nc,nk,d2brl,keymsg,herb) in SUM:
    if best is None or (isinstance(best,float) and best!=best):
        print("  skip figure for %s (no valid pose)"%name); continue
    safe=name.replace(" ","_")
    lig=ligand_atoms_from_pdbqt(os.path.join(OUT,"results",safe+"_out.pdbqt"))
    if not lig:
        print("  skip figure for %s (empty pose)"%name); continue
    rec=parse_pdb_atoms(rec_pdbqt)
    lxyz=np.array([[a['x'],a['y'],a['z']] for a in lig])
    rxyz=np.array([[a['x'],a['y'],a['z']] for a in rec])
    bxyz=BRL_ATOMS
    fig=plt.figure(figsize=(7.5,6))
    ax=fig.add_subplot(111,projection="3d")
    ax.scatter(rxyz[::3,0],rxyz[::3,1],rxyz[::3,2],s=1,c="#cccccc",alpha=0.30,label="PPARG LBD (chain A, apo)")
    # co-crystal Rosiglitazone (BRL) chain A — the TRUE orthosteric reference
    ax.scatter(bxyz[:,0],bxyz[:,1],bxyz[:,2],s=22,c="#17becf",marker="^",depthshade=False,label="Rosiglitazone (BRL, chain A)")
    # docked compound
    ax.scatter(lxyz[:,0],lxyz[:,1],lxyz[:,2],s=34,c="#d62728",marker="o",depthshade=False,label="%s (docked)"%name)
    # key residues
    recby={}
    for a in rec: recby.setdefault((a['resn'],a['resi']),[]).append((a['x'],a['y'],a['z']))
    keypts=[]
    for (rn,ri) in KEY:
        pts=recby.get((rn,ri))
        if not pts: continue
        c=np.array(pts).mean(axis=0); keypts.append(c)
        d=min(np.linalg.norm(lxyz-np.array(p),axis=1).min() for p in pts)
        col="#1f77b4" if d<=4.5 else "#bbbbbb"
        ax.scatter([c[0]],[c[1]],[c[2]],s=55,c=col,marker="*",depthshade=False)
        ax.text(c[0],c[1],c[2]," %s%d"%(rn,ri),color=col,fontsize=7)
    # centroids + honest verdict
    lc=lxyz.mean(axis=0); bc=bxyz.mean(axis=0); dcb=float(np.linalg.norm(lc-bc))
    verdict = "in rosiglitazone pocket" if dcb<8 else "LBD cavity, NOT rosiglitazone sub-pocket"
    ax.set_title("%s (%s) -> PPARG LBD chain A\nΔG=%.2f kcal/mol; pose-to-BRL = %.1f A -> %s"%(
        name, herb, best, dcb, verdict), fontsize=9)
    ax.set_xlabel("X (A)"); ax.set_ylabel("Y (A)"); ax.set_zlabel("Z (A)")
    # axis limits encompass ligand + BRL + key residues so spatial relationship is honest
    allp = np.vstack([lxyz, bxyz] + ([np.array(keypts)] if keypts else []))
    cen = allp.mean(axis=0); rad = float(np.abs(allp-cen).max()) + 6
    ax.set_xlim(cen[0]-rad,cen[0]+rad); ax.set_ylim(cen[1]-rad,cen[1]+rad); ax.set_zlim(cen[2]-rad,cen[2]+rad)
    ax.legend(fontsize=7, loc="upper left")
    plt.tight_layout()
    pf=os.path.join(OUT,"figures","Fig_S_%s.png"%safe)
    plt.savefig(pf, dpi=150); plt.close()
    print("wrote", pf)

# ---------- 6c. summary csv ----------
scsv=os.path.join(PROJ,"results/tables","s10_wuweizi_vina.csv")
with open(scsv,"w",newline="") as fh:
    w=csv.writer(fh)
    w.writerow(["compound","cn_name","chembl_id","herb","vina_dG_kcalmol",
                "n_pocket_contacts","n_key_residue_hits","min_dist_to_BRL_A","key_residue_distances"])
    for s in SUM:
        d2 = ("%.1f"%s[6]) if (isinstance(s[6],float) and s[6]==s[6]) else "NA"
        w.writerow([s[0],s[1],s[2],s[8],s[3],s[4],s[5],d2,s[7]])
print("wrote", scsv)

# ---------- print summary ----------
print("\n================ VINA SUMMARY (chain-A protocol v2) ================")
print("%-18s %-14s %8s %6s %6s %10s"%("compound","chembl","dG","cont","key","dist2BRL"))
for s in SUM:
    print("%-18s %-14s %8s %6d %6d %10s"%(s[0],s[2],("%.2f"%s[3] if s[3] is not None else "NA"),s[4],s[5],
          ("%.1f"%s[6] if isinstance(s[6],float) and s[6]==s[6] else "NA")))
print("Box center (chain-A BRL):", CENTER.round(2), " size %.1f A; Vina 1.2.7; seed 42; exhaustiveness 8"%BOX)
print("Receptor:", rec_pdbqt)
