#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
21_wuweizi_etcm_supplement.py
=============================================================================
【五味子木脂素靶点补充（ETCM 不可程序化访问时的等价方案）】

背景：
  - 2026-10-09 用 ChEMBL 拉成分-靶点，18 个五味子木脂素里只有 6 个（Schisandrin B /
    Schisantherin A / Schisanhenol / Gomisin A / Anwulignan / Schisandrol B）返回了
    靶点，其余 12 个（Schisandrin / Deoxyschisandrin / Schisandrin C / Schisantherin B /
    Gomisin G / Gomisin J,N / Pregomisin / Angeloylgomisin H / Tigloylgomisin P /
    Schisandrol A / Benzoylgomisin Q）活性=0。
  - ETCM(tcmip.cn) 网页可访问，但其 herb→ingredient→target 走 ThinkPHP 隐藏接口，
    WebFetch/urllib 均拿不到结构化数据；SwissTargetPrediction 异步任务服务端不返回
    jobId、TCMSP/HERB 下载页仅重定向。故改用**已发表、可引用**的五味子网络药理学
    文献作为补充来源（与 ETCM/HERB 的"文献+预测"性质一致），并在 manuscript 中明确标注。

数据来源（可引用）：
  [1] 基于网络药理学和分子对接技术的五味子治疗代谢相关脂肪性肝病研究.
      中国农业期刊集群. 结论：五味子酯乙/乙素/甲素/酚/醇乙/戈米辛J/当归酰戈米辛H
      为核心成分，关键靶点为 Akt1 与 PPAR-γ（Akt1、PPARG），并富集 AMPK、长寿调节等通路。
  [2] Schisandrin B Alleviates Lipid Metabolism Disorders ... via PPARγ-PCK1 ...
      Gao et al., Pharmaceuticals 2025, 18(10):1441.  (Schisandrin B → PPARG，实验验证：
      分子对接 + PPARG 选择性抑制剂 T0070907 反向验证)
  [3] Schisandrin B enhances autophagy ... through the PPARG-mediated PI3K/AKT/mTOR
      signaling pathway. Li et al., Free Radic Biol Med 2026. (同上，独立验证)

本脚本产出：
  - results/tables/s05b_wuweizi_literature_targets.csv  （7 核心木脂素 × {PPARG, AKT1}）
  - 把该层合并回 s05_compound_targets.csv（新增 source / reference 两列，ChEMBL 行标 ChEMBL）
=============================================================================
"""
import csv, os, shutil

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(BASE, "results", "tables")

# 7 个核心木脂素（见文献[1]），均与本项目 s05_compounds.csv 里的 18 味一一对应
CORE = ["Schisantherin B", "Schisandrin B", "Deoxyschisandrin", "Schisanhenol",
        "Schisandrol B", "Gomisin J", "Angeloylgomisin H"]

REF_NP = ("Schisandra chinensis NP (MAFLD network pharmacology, agrijournal 2024): "
          "7 core lignans target PPARG & AKT1")
REF_EXP = ("Gao et al. Pharmaceuticals 2025,18(10):1441 + Li et al. Free Radic Biol Med 2026 "
           "(Schisandrin B -> PPARG, docking + T0070907 inhibitor validation)")

# 读 s05_compounds.csv 取 cn_name / chembl_id（便于溯源）
comp = {}
with open(os.path.join(TAB, "s05_compounds.csv"), encoding="utf-8") as f:
    for r in csv.DictReader(f):
        comp[r["query"]] = r

def info(name):
    c = comp.get(name, {})
    return c.get("cn_name", ""), c.get("chembl_id", "")

# 构造文献层行
targets = {
    "PPARG": ("Peroxisome proliferator-activated receptor gamma", "SINGLE PROTEIN"),
    "AKT1":  ("RAC-alpha serine/threonine-protein kinase", "SINGLE PROTEIN"),
}
rows = []
for cmp in CORE:
    cn, cid = info(cmp)
    for g, (tname, ttype) in targets.items():
        # Schisandrin B 有实验级 PPARG 证据，其余为核心网络药理学预测
        ref = REF_EXP if (cmp == "Schisandrin B" and g == "PPARG") else REF_NP
        rows.append({
            "herb": "五味子",
            "compound": cmp,
            "cn_name": cn,
            "chembl_id": cid,
            "gene_symbol": g,
            "best_nM": "",                       # 文献层无 ChEMBL 定量，留空→NA
            "assay_types": "NetworkPharmacology",
            "target_name": tname,
            "target_type": ttype,
            "source": "Literature_NP",
            "reference": ref,
        })

out = os.path.join(TAB, "s05b_wuweizi_literature_targets.csv")
fields = ["herb", "compound", "cn_name", "chembl_id", "gene_symbol", "best_nM",
          "assay_types", "target_name", "target_type", "source", "reference"]
with open(out, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader(); w.writerows(rows)
print("✓ 写出文献层:", out, " (%d 行 = %d 成分 × %d 靶点)" % (len(rows), len(CORE), len(targets)))

# ---- 合并回 s05_compound_targets.csv（先备份） ----
src = os.path.join(TAB, "s05_compound_targets.csv")
bak = os.path.join(TAB, "s05_compound_targets.before_wuweizi.csv")
shutil.copy(src, bak)
print("✓ 已备份原表:", bak)

with open(src, encoding="utf-8") as f:
    old = list(csv.DictReader(f))

# 给原 ChEMBL 行补 source/reference
for r in old:
    r.setdefault("source", "ChEMBL")
    r.setdefault("reference", "ChEMBL snapshot 2026-10-09")

merged = old + rows
# 统一列顺序
allfields = ["herb", "compound", "cn_name", "chembl_id", "gene_symbol", "best_nM",
             "assay_types", "target_name", "target_type", "source", "reference"]
with open(src, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=allfields)
    w.writeheader()
    for r in merged:
        w.writerow({k: r.get(k, "") for k in allfields})

# 统计
genes = sorted({r["gene_symbol"] for r in merged})
wu = [r for r in merged if r["source"] == "Literature_NP"]
print("合并后：成分-靶点对 %d 行 | 去重靶点 %d 个" % (len(merged), len(genes)))
print("其中五味子文献层新增：%d 行，涉及靶点 %s" %
      (len(wu), sorted({r["gene_symbol"] for r in wu})))
print("ChEMBL 原行：%d | 五味子文献层：%d" % (len(old), len(wu)))
