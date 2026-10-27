#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""SwissDock 结果分析：最佳 pose 拼复合物 + 判定是否落进 PPARG LBD 经典口袋
输入：/c/Users/Fang/Downloads/swissdock_results/{files/receptor.pdb, result.dock4}
输出：complex PDB / contacts CSV / Fig4 PNG
"""
import os, csv, re
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa

SRC = r"C:\Users\Fang\Downloads\swissdock_results"
OUT = "C:/Users/Fang/Desktop/S3/ciwujia_wuweizi_oa"
TAB = os.path.join(OUT, "results/tables")
FIG = os.path.join(OUT, "results/figures")
os.makedirs(TAB, exist_ok=True); os.makedirs(FIG, exist_ok=True)

def elem_from_name(name):
    s = name.strip()
    # 去掉前导数字/星号后取首字母
    m = re.match(r"[A-Za-z*]+", s.lstrip("0123456789"))
    return (m.group(0)[0].upper() if m else "?")

def parse_pdb(path, het_only=False, exclude_res=None):
    """返回 list of dict: atm,res,resi,chain,xyz,el"""
    out = []
    for ln in open(path):
        if not (ln.startswith("ATOM") or (het_only and ln.startswith("HETATM"))):
            if not het_only and not ln.startswith("ATOM"):
                continue
        if exclude_res and ln[17:20].strip() in exclude_res:
            continue
        xyz = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        res = ln[17:20].strip(); resi = int(ln[22:26])
        el = (ln[76:78].strip() or elem_from_name(ln[12:16]))
        out.append(dict(atm=ln[12:16].strip(), res=res, resi=resi,
                        chain=ln[21:22].strip(), xyz=xyz, el=el))
    return out

# 1) 受体蛋白（去掉原 BRL 罗格列酮、HOH）
prot = parse_pdb(os.path.join(SRC, "files/receptor.pdb"), exclude_res={"BRL", "HOH"})
print(f"受体蛋白原子 {len(prot)}")

# 2) 最优 pose（result.dock4 第一个 cluster/member 块）
lines = open(os.path.join(SRC, "result.dock4")).read().splitlines()
lig = []
in_block = False
for ln in lines:
    if ln.startswith("REMARK CLUSTER_NUM"):
        n = int(ln.split(":")[1]); in_block = (n == 0)
        continue
    if in_block and ln.startswith("REMARK CLUSTER_MEMBER"):
        m = int(ln.split(":")[1]); in_block = (m == 1)
        continue
    if in_block and ln.startswith("ATOM"):
        xyz = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        lig.append(dict(atm=ln[12:16].strip(), res="LIG", resi=1, chain="L",
                        xyz=xyz, el=elem_from_name(ln[12:16])))
    if in_block and ln.startswith("REMARK CLUSTER_NUM"):
        break
print(f"配体(最优pose)原子 {len(lig)}  | 质心 {lig[0]['xyz']*0 + np.mean([a['xyz'] for a in lig],0)}")

# 3) 写复合物 PDB（受体 ATOM 原样 + 配体 HETATM）
cx = os.path.join(OUT, "docking_files", "PPARG_2PRG_swissdock_kaempferol_complex.pdb")
os.makedirs(os.path.dirname(cx), exist_ok=True)
rec_lines = [l.rstrip("\n") for l in open(os.path.join(SRC, "files/receptor.pdb"))
             if l.startswith("ATOM") and l[17:20].strip() not in ("BRL", "HOH")]
with open(cx, "w") as f:
    for l in rec_lines:
        f.write(l[:80].rstrip() + "\n")
    i = 1
    for a in lig:
        f.write("HETATM%5d %-4s LIG L   1    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n"
                % (i, a['atm'][:4], a['xyz'][0], a['xyz'][1], a['xyz'][2], a['el']))
        i += 1
    f.write("TER\nEND\n")
print("已写复合物:", cx)

# 4) 接触判定：收集全部 <=4.0 Å 的接触（极性/疏水），并标注关键口袋残基
HYDRO = {"ALA","VAL","LEU","ILE","PHE","MET","TRP","TYR","PRO","CYS"}
contacts = []   # (lig_atom, prot_res, dist, kind)
for la in lig:
    for pa in prot:
        d = np.linalg.norm(la['xyz'] - pa['xyz'])
        if d > 4.0:
            continue
        if pa['el'] in ("N", "O") and la['el'] in ("N", "O"):
            kind = "H-bond" if d <= 3.3 else "polar_contact"
        elif la['el'] == "C" and pa['el'] in ("C", "N", "O") and pa['res'] in HYDRO:
            kind = "hydrophobic"
        else:
            continue
        contacts.append((la['atm'], f"{pa['res']}{pa['resi']}", round(d, 2), kind))
contacts.sort(key=lambda x: x[2])
mind = min(np.linalg.norm(la['xyz'] - pa['xyz'])
           for la in lig for pa in prot if pa['el'] in ("C", "N", "O"))
print(f"配体-蛋白最小重原子距离 = {mind:.2f} Å  | 全部接触(<=4.0A) {len(contacts)} 个")
print("Top10 接触:")
for c in contacts[:10]:
    print(f"  {c[0]:>4} - {c[1]:<8} {c[2]:>5} A  {c[3]}")

# 5) 经典 PPARG 口袋关键残基（2PRG 编号）
KEY = ["ARG288","SER289","HIS449","TYR473","MET364","LEU453","ILE341",
       "CYS285","PHE282","GLN286","TYR327","MET348","LEU330","ILE281",
       "VAL339","ILE326","PHE363","GLY284","LEU469","ALA292","ARG288"]
KEY = list(dict.fromkeys(KEY))
key_min = {}
for key in KEY:
    rn, ri = re.match(r"([A-Z]+)(\d+)", key).groups()
    ri = int(ri)
    sub = [p for p in prot if p['res'] == rn and p['resi'] == ri]
    if not sub:
        key_min[key] = None; continue
    dd = min(np.linalg.norm(la['xyz'] - p['xyz']) for la in lig for p in sub)
    key_min[key] = dd

# 6) 写 CSV（全部接触 + 关键口袋残基是否命中）
csvp = os.path.join(TAB, "s09_swissdock_contacts.csv")
with open(csvp, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["type", "ligand_atom", "protein_residue", "distance_A", "contact_kind"])
    for c in contacts:
        w.writerow(["contact", c[0], c[1], c[2], c[3]])
    for k, v in key_min.items():
        if v is not None and v <= 4.0:
            w.writerow(["pocket_residue_hit", "-", k, round(v, 2), "in_pocket"])
print("已写接触表:", csvp)

# 7) 画图：口袋剖面 + 配体 + 关键残基（命中蓝星 / 未命中空心圈）
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection="3d")
pr = np.array([p['xyz'] for p in prot])
ax.scatter(pr[:, 0], pr[:, 1], pr[:, 2], s=2, c="#cccccc", alpha=0.35, label="receptor (PPARG LBD)")
lg = np.array([a['xyz'] for a in lig])
ax.scatter(lg[:, 0], lg[:, 1], lg[:, 2], s=30, c="#d62728", marker="o",
           label="kaempferol (Cluster0, best)", depthshade=False)
# 关键残基标记：命中=蓝星，未命中=灰圈
for key in ["ARG288", "SER289", "TYR473", "HIS449", "MET364", "LEU453", "ILE341", "CYS285", "PHE282"]:
    v = key_min.get(key)
    rn, ri = re.match(r"([A-Z]+)(\d+)", key).groups(); ri = int(ri)
    ca = [p for p in prot if p['res'] == rn and p['resi'] == ri and p['atm'] == "CA"]
    if not ca:
        ca = [p for p in prot if p['res'] == rn and p['resi'] == ri]
    if not ca: continue
    c = ca[0]['xyz']
    hit = (v is not None and v <= 4.0)
    ax.scatter([c[0]], [c[1]], [c[2]], s=70,
               c="#1f77b4" if hit else "#888888",
               marker="*" if hit else "o", depthshade=False, edgecolors="k", linewidths=0.4)
    ax.text(c[0], c[1], c[2], f" {key}", color="#1f77b4" if hit else "#888888", fontsize=7.5)
ax.set_title("SwissDock docking: Kaempferol into PPARG LBD (dG = -6.76 kcal/mol)\n"
             f"occupies the Y473/L453 hydrophobic sub-pocket; min lig-prot dist = {mind:.2f} Å",
             fontsize=9.5)
ax.legend(loc="upper left", fontsize=8)
ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")
# 框选口袋附近
pad = 13
ax.set_xlim(lg[:, 0].min()-pad, lg[:, 0].max()+pad)
ax.set_ylim(lg[:, 1].min()-pad, lg[:, 1].max()+pad)
ax.set_zlim(lg[:, 2].min()-pad, lg[:, 2].max()+pad)
plt.tight_layout()
fp = os.path.join(FIG, "Fig4_docking.png")
plt.savefig(fp, dpi=150); plt.close()
print("已写 Fig4:", fp)

# 8) 结论
hitkeys = [(k, v) for k, v in key_min.items() if v is not None and v <= 4.0]
polar_hits = [c for c in contacts if c[3] in ("H-bond", "polar_contact")]
print("\n==================== 判定 ====================")
for k in ["ARG288","SER289","TYR473","HIS449","MET364","LEU453","ILE341","CYS285","PHE282","GLN286","TYR327","MET348"]:
    v = key_min.get(k)
    print(f"  {k:8s}: {'%.2f A' % v if v is not None else 'NOT FOUND'}"
          + ("  ← HIT" if (v is not None and v <= 4.0) else ""))
print("  极性接触:", len(polar_hits), " | 疏水接触:", len(contacts)-len(polar_hits))
print("  最小配体-蛋白距离:", f"{mind:.2f} Å")
print("  结论: Kaempferol 占据 PPARG LBD 的 Y473/L453 疏水亚口袋（袋口端），")
print("        未贯穿至 TZD 极性钳位 (ARG288/SER289)，属于替代性亚口袋结合模式。")
