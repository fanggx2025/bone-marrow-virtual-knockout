# -*- coding: utf-8 -*-
"""Parse SwissDock AC (EADock) result.dock4 for SchisandrolB -> PPARG LBD chainA."""
import math, re
from collections import defaultdict

BASE = r"C:\Users\Fang\Desktop\S3\ciwujia_wuweizi_oa\docking_files\wuweizi_vina\web_submission\swissdock_results\job5_AC"
BRL_CENTROID = (49.72, -36.98, 19.29)  # chainA rosiglitazone centroid

# ---- 1. parse poses ----
poses = []  # list of dict(cluster, spdG, atoms=[(name,resn,resi,chain,x,y,z)])
cur = None
with open(BASE + r"\result.dock4", "r", errors="ignore") as f:
    for line in f:
        if line.startswith("REMARK SP-dG"):
            spdG = float(line.split(":")[1].strip())
            if cur is None:
                cur = {"spdG": spdG, "atoms": []}
            else:
                cur["spdG"] = min(cur.get("spdG", spdG), spdG)
        elif line.startswith("REMARK CLUSTER_NUM"):
            cn = int(line.split(":")[1].strip())
            if cur is None:
                cur = {"atoms": []}
            cur["cluster"] = cn
        elif line.startswith(("ATOM  ", "HETATM")):
            if cur is None:
                cur = {"atoms": []}
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            cur["atoms"].append((line[12:16].strip(), x, y, z))
        elif line.startswith("TER"):
            if cur and cur.get("atoms"):
                poses.append(cur)
            cur = None
if cur and cur.get("atoms"):
    poses.append(cur)

print(f"total pose blocks: {len(poses)}")

# ---- 2. cluster stats ----
clusters = defaultdict(list)
for p in poses:
    clusters[p.get("cluster", -1)].append(p)
stat = []
for cn, ps in clusters.items():
    vals = [p["spdG"] for p in ps if "spdG" in p]
    if not vals:
        continue
    stat.append((min(vals), cn, len(ps)))
stat.sort()
print("top clusters (best SP-dG, cluster#, nPoses):")
for s in stat[:5]:
    print("  ", s)

# ---- 3. best pose overall: distance to BRL centroid + contacts ----
best = min(poses, key=lambda p: p.get("spdG", 0) if "spdG" in p else 1e9)
cx = sum(a[1] for a in best["atoms"]) / len(best["atoms"])
cy = sum(a[2] for a in best["atoms"]) / len(best["atoms"])
cz = sum(a[3] for a in best["atoms"]) / len(best["atoms"])
dist = math.dist((cx, cy, cz), BRL_CENTROID)
print(f"\nBEST pose: cluster {best.get('cluster')}, SP-dG {best.get('spdG')}")
print(f"pose centroid: ({cx:.2f}, {cy:.2f}, {cz:.2f}), dist to chainA-BRL centroid = {dist:.2f} A")

# ---- 4. contacts with receptor (<4.5 A heavy atoms) ----
rec = []
with open(BASE + r"\files\receptor.pdb", "r", errors="ignore") as f:
    for line in f:
        if line.startswith("ATOM  "):
            rec.append((line[17:20].strip(), int(line[22:26]),
                        float(line[30:38]), float(line[38:46]), float(line[46:54])))
contacts = defaultdict(float)
for rn, ri, x, y, z in rec:
    dmin = min(math.dist((a[1], a[2], a[3]), (x, y, z)) for a in best["atoms"])
    if dmin < 4.5:
        contacts[f"{rn}{ri}"] = round(dmin, 1)
top = sorted(contacts.items(), key=lambda kv: kv[1])
print("contacts (<4.5A):", ", ".join(f"{k}({v})" for k, v in top))

# ---- 5. key pocket residues check ----
for key in ["ARG288", "SER289", "CYS285", "MET364", "HIS323", "HIS449", "TYR473", "LEU330", "ILE341", "MET348"]:
    hit = contacts.get(key)
    if hit:
        print(f"  {key}: {hit} A  <-- contacted")
