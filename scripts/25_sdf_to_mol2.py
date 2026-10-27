# -*- coding: utf-8 -*-
"""SDF (V2000) -> Tripos MOL2 转换器（用于 SwissDock 上传）
显式 H 原子剔除（对接服务器会自行加氢/质子化）。

原子类型推断（Kekulé 感知）:
  C 环内 且 与环内原子成双键  -> C.ar   （芳环碳，含 Kekulé 形式）
  C 双键连氧(羰基)           -> C.2
  其余 C                     -> C.3
  O 双键(C=O)                -> O.2
  O 在环内(醚氧, 如吡喃环)    -> O.2
  其余 O (羟基)              -> O.3
键类型: 1->1, 2->2, 3->3, 4->ar
"""
import sys, os, glob

def parse_sdf(path):
    lines = open(path, encoding='utf-8', errors='ignore').read().splitlines()
    na = int(lines[3][0:3]); nb = int(lines[3][3:6])
    atoms = []
    for i in range(na):
        ln = lines[4 + i]
        atoms.append((ln[31:34].strip(), float(ln[0:10]), float(ln[10:20]), float(ln[20:30])))
    bonds = []
    for i in range(nb):
        ln = lines[4 + na + i]
        bonds.append((int(ln[0:3]), int(ln[3:6]), int(ln[6:9])))
    return lines[0].strip() or os.path.basename(path), atoms, bonds

def convert(path, out):
    name, atoms, bonds = parse_sdf(path)
    # 剔除显式 H
    keep = [i + 1 for i, a in enumerate(atoms) if a[0] != 'H']
    remap = {old: new for new, old in enumerate(keep, 1)}
    ats = [atoms[o - 1] for o in keep]
    bnds = [(remap[a], remap[b], t) for a, b, t in bonds
            if a in remap and b in remap and atoms[a-1][0] != 'H' and atoms[b-1][0] != 'H']
    # 双键集合
    dbl = set()
    for a, b, t in bnds:
        if t in (2, 4):
            dbl.add(a); dbl.add(b)
    # 环感知：反复剪掉度为 1 的原子，剩下的就是环内原子
    ring = set()
    for a, b, _ in bnds:
        ring.add(a); ring.add(b)
    changed = True
    while changed:
        changed = False
        for x in list(ring):
            d = sum(1 for a, b, _ in bnds if (a == x and b in ring) or (b == x and a in ring))
            if d <= 1:
                ring.discard(x); changed = True
    # 芳香碳 = 环内 且 与环内原子成双键（Kekulé 芳环 / type-4 均覆盖）
    ringpair = set()
    for a, b, t in bnds:
        if t in (2, 4) and a in ring and b in ring:
            ringpair.add(a); ringpair.add(b)
    elem = {new: a[0] for new, a in enumerate(ats, 1)}
    def atype(idx, el):
        if el == 'C':
            if idx in ringpair:            return 'C.ar'
            if idx in dbl:                 return 'C.2'
            return 'C.3'
        if el == 'O':
            if idx in dbl or idx in ring:  return 'O.2'
            return 'O.3'
        if el == 'N':
            if idx in ringpair:            return 'N.ar'
            if idx in dbl:                 return 'N.2'
            return 'N.3'
        if el == 'S':
            return 'S.o2' if idx in dbl else 'S.3'
        return el
    with open(out, 'w') as f:
        f.write('@<TRIPOS>MOLECULE\n%s\n %d %d 0 0 0\nNO_CHARGES\n\n' % (name, len(ats), len(bnds)))
        f.write('@<TRIPOS>ATOM\n')
        for new, (el, x, y, z) in enumerate(ats, 1):
            f.write('%5d %-8s %10.4f %10.4f %10.4f %-6s %3d\n' % (new, el + str(new), x, y, z, atype(new, el), 1))
        f.write('@<TRIPOS>BOND\n')
        for i, (a, b, t) in enumerate(bnds, 1):
            bt = {1: '1', 2: '2', 3: '3', 4: 'ar'}.get(t, '1')
            f.write('%5d %5d %5d %s\n' % (i, a, b, bt))
    print('OK %-42s -> %d 原子 / %d 键' % (os.path.basename(out), len(ats), len(bnds)))

if __name__ == '__main__':
    base = os.path.dirname(os.path.abspath(__file__))
    ldir = os.path.normpath(os.path.join(base, '..', 'docking_files', 'ligands'))
    for sdf in sorted(glob.glob(os.path.join(ldir, '*_3D.sdf'))):
        convert(sdf, sdf[:-4] + '.mol2')
