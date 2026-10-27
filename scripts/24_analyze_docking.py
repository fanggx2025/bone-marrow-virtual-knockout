#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""分析 CB-Dock2 输出的 PPARG_Kaempferol 复合物 PDB（修正残基解析）"""
import csv, os, sys, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

SRC = sys.argv[1] if len(sys.argv) > 1 else \
    "C:/Users/Fang/Downloads/PPARG_2PRG_clean_Kaempferol_CHEMBL150_3D_fitdock-out_1.-5.9.complex.pdb"
OUTDIR = "C:/Users/Fang/Desktop/S3/ciwujia_wuweizi_oa"
TAB = os.path.join(OUTDIR, "results/tables"); FIG = os.path.join(OUTDIR, "results/figures")
os.makedirs(TAB, exist_ok=True); os.makedirs(FIG, exist_ok=True)

def parse(lines):
    out = []
    for ln in lines:
        xyz = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        el = ln[76:78].strip() or ln[12:16].strip()[0]
        out.append(dict(atm=ln[12:16].strip(), res=ln[17:20].strip(),
                        resi=int(ln[22:26]), xyz=xyz, el=el,
                        chain=ln[21:22].strip()))
    return out

raw = open(SRC).read().splitlines()
prot = parse([l for l in raw if l.startswith("ATOM")])
lig  = parse([l for l in raw if l.startswith("HETATM")])
print(f"蛋白原子 {len(prot)} | 配体原子 {len(lig)}")

# key residues 字典：(resname, resnum) -> label
KEY = {("ARG",288):"R288",("TYR",473):"Y473",("HIS",323):"H323",
       ("HIS",466):"H466",("SER",289):"S289",("SER",342):"S342",
       ("LYS",264):"K264",("GLN",286):"Q286",("MET",348):"M348",
       ("LEU",333):"L333",("PHE",360):"F360",("ILE",341):"I341",
       ("GLU",343):"E343",("LEU",330):"L330",("VAL",339):"V339",
       ("MET",364):"M364",("GLY",344):"G344",("SER",245):"S245"}
HYDRO = {"ALA","VAL","LEU","ILE","PHE","MET","TRP","TYR","PRO","CYS"}

polar=[]; hydro=[]; key_hit={}
for la in lig:
    for pa in prot:
        d = np.linalg.norm(la['xyz']-pa['xyz'])
        if pa['el'] in ("N","O"):
            if d <= 3.5:
                kind = "polar_contact"
                if d <= 3.3:
                    kind = "H-bond" if (la['el']=="O" and pa['el'] in ("N","O")) else "polar_contact"
                polar.append((la['atm'], f"{pa['res']}{pa['resi']}", round(d,2), kind))
                lab = KEY.get((pa['res'],pa['resi']))
                if lab: key_hit.setdefault(lab,[]).append((la['atm'],round(d,2),kind))
        elif la['el']=="C" and pa['el']=="C" and pa['res'] in HYDRO:
            if d <= 4.0:
                hydro.append((la['atm'], f"{pa['res']}{pa['resi']}", round(d,2)))

# 配体最近蛋白原子（整体贴合度）
mind = min(np.linalg.norm(la['xyz']-pa['xyz']) for la in lig for pa in prot
           if pa['el'] in ("C","N","O"))
print(f"配体-蛋白最小重原子距离 = {mind:.2f} Å")

# 写 CSV
csvp = os.path.join(TAB,"s09_docking_contacts.csv")
with open(csvp,"w",newline="",encoding="utf-8") as f:
    w=csv.writer(f); w.writerow(["ligand_atom","protein_residue","distance_A","interaction"])
    for r in sorted(polar,key=lambda x:x[2]):
        w.writerow(r)
print(f"极性接触 {len(polar)} 条 (H-bond {sum(1 for x in polar if x[3]=='H-bond')}), 疏水 {len(hydro)} 条")
print("CSV ->", csvp)

print("\n===== 经典 PPARG 关键残基命中 =====")
for lab in ["R288","Y473","H323","H466","S289","S342","K264","Q286","M348","L333","F360","E343"]:
    if lab in key_hit:
        best=min(key_hit[lab],key=lambda t:t[1])
        print(f"  {lab:5s} <- {best[0]:4s} d={best[1]:.2f}A ({best[2]})  共{len(key_hit[lab])}接触")
    else:
        print(f"  {lab:5s} 未直接命中")

# ========== Fig4 ==========
O_pos = {"oxo4":(0.0,-1.55),"oh3":(0.6,0.95),"oh5":(-2.7,0.95),"oh7":(-2.7,-0.95),"oh4p":(3.0,0.0)}
def hexes(cx,r=1.0):
    return [(cx+r*np.cos(np.radians(a)),r*np.sin(np.radians(a))) for a in [0,60,120,180,240,300]]
skel=hexes(-2.0)+hexes(0.0)+hexes(2.0)
edges=[(i,(i+1)%6) for i in range(6)]+[(6+i,(6+i+1)%6) for i in range(6)]+[(12+i,(12+i+1)%6) for i in range(6)]

fig=plt.figure(figsize=(11,5.2))
ax1=fig.add_axes([0.04,0.08,0.46,0.84]); ax1.set_aspect("equal"); ax1.axis("off")
ax2=fig.add_axes([0.56,0.12,0.40,0.80])
for (i,j) in edges:
    x1,y1=skel[i]; x2,y2=skel[j]
    ax1.plot([x1,x2],[y1,y2],color="#444",lw=1.6,zorder=1)
for k,(x,y) in O_pos.items():
    ax1.scatter([x],[y],s=70,color="#C0392B",zorder=3,edgecolor="white")
    ax1.text(x,y+0.18,{"oxo4":"O(4-oxo)","oh3":"OH(3)","oh5":"OH(5)","oh7":"OH(7)","oh4p":"OH(4')"}[k],
             ha="center",va="bottom",fontsize=7.5,color="#C0392B")
ax1.text(-2.0,1.35,"A",ha="center",fontsize=10,fontweight="bold")
ax1.text(0.0,1.35,"C",ha="center",fontsize=10,fontweight="bold")
ax1.text(2.0,1.35,"B",ha="center",fontsize=10,fontweight="bold")
# 实际接触残基汇总（无论是否关键）
ctxt = "Observed contacts: " + (", ".join(sorted({r[1] for r in polar})) or "none")
ax1.text(-3.4,-2.6,ctxt,fontsize=7.5,color="#555" if polar else "#C0392B")

# 标注命中关键残基（按官能团就近放）
annot={"oxo4":["R288","S289","H466"],"oh3":["K264","Q286"],"oh5":["M348","L333"],"oh7":["F360","G344"],"oh4p":["Y473","H323"]}
for k,(x,y) in O_pos.items():
    labs=[l for l in annot[k] if l in key_hit]
    if labs:
        txt=", ".join(labs)
        ax1.annotate(txt,xy=(x,y),xytext=(x+0.8,y+(0.5 if y>=0 else -0.6)),
                     fontsize=7,color="#1A5276",arrowprops=dict(arrowstyle="->",color="#1A5276",lw=0.8))
ONP = all(k in key_hit for k in ["R288","Y473"]) or "Y473" in key_hit or "R288" in key_hit
verdict = "IN canonical PPARγ pocket (R288/Y473 contacted)" if ("Y473" in key_hit or "R288" in key_hit) else \
          "OFF-POCKET: missing R288/Y473/H323/H466 contacts (NOT usable for Fig4)"
ax1.set_title("A) Kaempferol in PPARγ LBD\n%s" % verdict, fontsize=9.0, loc="left")

agg={}
for r in polar:
    if r[1] not in agg or r[2]<agg[r[1]][0]: agg[r[1]]=(r[2],r[3],r[0])
items=sorted(agg.items(),key=lambda kv:kv[1][0])[:18]
ax2.barh(range(len(items)),[v[0] for _,v in items],
         color=["#C0392B" if v[1]=="H-bond" else "#2E86C1" for _,v in items])
ax2.set_yticks(range(len(items))); ax2.set_yticklabels([k for k,_ in items],fontsize=7)
ax2.invert_yaxis(); ax2.set_xlabel("min distance to Kaempferol (Å)",fontsize=8)
ax2.axvline(3.3,color="#C0392B",ls="--",lw=0.8)
ax2.set_title(f"B) Top polar contacts (n={len(polar)}; H-bond red)\nCB-Dock2 score = -5.9 kcal/mol",fontsize=9.5,loc="left")
ax2.tick_params(labelsize=7)
plt.savefig(os.path.join(FIG,"Fig4_docking.png"),dpi=160,bbox_inches="tight")
print("Fig4 ->", os.path.join(FIG,"Fig4_docking.png"))
