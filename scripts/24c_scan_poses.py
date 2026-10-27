#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""扫描 result.dock4 所有 pose，按"是否落进 PPARG 深袋"重排，挑最佳 on-pocket pose"""
import os, re
import numpy as np

SRC = r"C:\Users\Fang\Downloads\swissdock_results"
REC = os.path.join(SRC, "files/receptor.pdb")
DOCK = os.path.join(SRC, "result.dock4")

def elem(name):
    m = re.match(r"[A-Za-z*]+", name.strip().lstrip("0123456789"))
    return m.group(0)[0].upper() if m else "?"

def parse_prot(path):
    out = []
    for ln in open(path):
        if not ln.startswith("ATOM"): continue
        if ln[17:20].strip() in ("BRL","HOH"): continue
        xyz = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        out.append(dict(res=ln[17:20].strip(), resi=int(ln[22:26]),
                        el=(ln[76:78].strip() or elem(ln[12:16])), xyz=xyz))
    return out

prot = parse_prot(REC)
# 关键残基坐标（选 CA 或任意代表）
def res_atoms(rn, ri):
    return [p for p in prot if p['res']==rn and p['resi']==ri]
BRL_C = np.array([54.7, -21.0, 30.5])   # 共晶罗格列酮质心

# 解析所有 pose 块
blocks = []
cur = None
for ln in open(DOCK):
    if ln.startswith("REMARK CLUSTER_NUM"):
        if cur: blocks.append(cur)
        cur = {"cn": int(ln.split(":")[1]), "atoms": []}
    elif ln.startswith("REMARK CLUSTER_MEMBER"):
        if cur: cur["cm"] = int(ln.split(":")[1])
    elif cur is not None and ln.startswith("ATOM"):
        cur["atoms"].append(np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]))
if cur: blocks.append(cur)

print(f"共解析 pose 块 {len(blocks)}")

# 对每个 pose 打分
rows = []
for b in blocks:
    if not b["atoms"]: continue
    arr = np.array(b["atoms"])
    cen = arr.mean(0)
    mind = min(np.linalg.norm(a-p['xyz']) for a in arr for p in prot if p['el'] in ("C","N","O"))
    # 到关键残基距离
    def dres(rn, ri):
        s = res_atoms(rn, ri)
        if not s: return None
        return min(np.linalg.norm(a-p['xyz']) for a in arr for p in s)
    d_arg288 = dres("ARG",288); d_ser289 = dres("SER",289)
    d_tyr473 = dres("TYR",473); d_his449 = dres("HIS",449)
    d_met364 = dres("MET",364); d_leu453 = dres("LEU",453)
    d_ile341 = dres("ILE",341)
    hits = sum(1 for d in [d_arg288,d_ser289,d_tyr473,d_his449,d_met364,d_leu453,d_ile341]
               if d is not None and d <= 4.0)
    db = np.linalg.norm(cen - BRL_C)
    rows.append(dict(cn=b["cn"], cm=b.get("cm"), n=len(arr), mind=round(mind,2),
                     cen=tuple(round(x,1) for x in cen),
                     dArg288=(round(d_arg288,2) if d_arg288 is not None else None),
                     dSer289=(round(d_ser289,2) if d_ser289 is not None else None),
                     dTyr473=(round(d_tyr473,2) if d_tyr473 is not None else None),
                     dHis449=(round(d_his449,2) if d_his449 is not None else None),
                     dMet364=(round(d_met364,2) if d_met364 is not None else None),
                     dLeu453=(round(d_leu453,2) if d_leu453 is not None else None),
                     hits=hits, dToBRL=round(db,2)))

# 选 on-pocket：dArg288<=4.5 且 dTyr473<=4.5 优先；否则按 hits 降序、dToBRL 升序
def score(r):
    on_pocket = (r["dArg288"] is not None and r["dArg288"]<=4.5 and
                 r["dTyr473"] is not None and r["dTyr473"]<=4.5)
    return (on_pocket, r["hits"], -r["dToBRL"])
rows.sort(key=score, reverse=True)

print(f"\n{'clu':>3} {'mem':>3} {'n':>3} {'minD':>5} {'ARG288':>7} {'SER289':>7} {'TYR473':>7} {'HIS449':>7} {'MET364':>7} {'LEU453':>7} {'hits':>4} {'dToBRL':>7}  cen")
for r in rows[:20]:
    print(f"{r['cn']:>3} {r['cm']:>3} {r['n']:>3} {r['mind']:>5} "
          f"{str(r['dArg288']):>7} {str(r['dSer289']):>7} {str(r['dTyr473']):>7} {str(r['dHis449']):>7} "
          f"{str(r['dMet364']):>7} {str(r['dLeu453']):>7} {r['hits']:>4} {r['dToBRL']:>7}  {r['cen']}")

# 保存 top on-pocket pose 的 cluster/member 供后续提取
best = rows[0]
print("\nTOP:", best)
