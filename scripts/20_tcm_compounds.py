#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
20_tcm_compounds.py —— 刺五加 / 五味子 成分与靶点获取
=============================================================================
【为什么绕开 NCBI / PubChem / UniProt】
  2026-10-08 实测：
    不通：PubChem(NCBI)、UniProt REST、CTD、KEGG REST、STITCH
    通  ：ChEMBL(EBI)、mygene.info、STRING、SwissTargetPrediction、HERB、Enrichr
  所以：化合物走 ChEMBL，UniProt accession → gene symbol 走 mygene.info。

【流程】
  1. 文献整理的成分清单（两味药各 ~18 个主要活性成分）
  2. ChEMBL 按名检索 → molecule_chembl_id / canonical_smiles / 理化性质
  3. ChEMBL activity → 人源靶点（standard_type ∈ IC50/Ki/Kd/EC50，value ≤ 10 µM）
  4. ChEMBL target → UniProt accession → mygene.info → gene symbol
  5. 输出 CSV，供后续与虚拟敲除结果取交集（Fig6 三方共识）
"""

import json, time, urllib.parse, urllib.request, urllib.error, os, sys, csv

BASE = "https://www.ebi.ac.uk/chembl/api/data"
MYGENE = "https://mygene.info/v3"
UA = {"User-Agent": "Mozilla/5.0 (research script)", "Accept": "application/json"}
TIMEOUT = 30


def get(url, tries=5):
    """【坑】ChEMBL 在高频请求下会间歇性返回 500，短重试没用，要拉长退避。"""
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
                return json.loads(r.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            last = e
            if e.code in (500, 502, 503, 429):
                time.sleep(3 + 4 * i)
                continue
            sys.stderr.write("  ! GET failed: %s  (%s)\n" % (url[:110], e))
            return None
        except Exception as e:
            last = e
            time.sleep(2 + 2 * i)
    sys.stderr.write("  ! GET failed after %d tries: %s  (%s)\n" % (tries, url[:110], last))
    return None


# ---------------------------------------------------------------- 成分清单
CIWUJIA = [
    ("Syringin",             "刺五加苷B/紫丁香苷"),
    ("Eleutheroside E",      "刺五加苷E"),
    ("Eleutheroside A",      "刺五加苷A"),
    ("Isofraxidin",          "异嗪皮啶"),
    ("Hyperoside",           "金丝桃苷"),
    ("Rutin",                "芦丁"),
    ("Chlorogenic acid",     "绿原酸"),
    ("Caffeic acid",         "咖啡酸"),
    ("Quercetin",            "槲皮素"),
    ("Kaempferol",           "山柰酚"),
    ("Sesamin",              "芝麻素"),
    ("beta-Sitosterol",      "β-谷甾醇"),
    ("Stigmasterol",         "豆甾醇"),
    ("Daucosterol",          "胡萝卜苷"),
    ("Syringaresinol",       "丁香树脂酚"),
    ("Liriodendrin",         "鹅掌楸苷"),
    ("Protocatechuic acid",  "原儿茶酸"),
    ("Scopoletin",           "东莨菪内酯"),
]

WUWEIZI = [
    ("Schisandrin",          "五味子醇甲"),
    ("Schisandrol B",        "五味子醇乙"),
    ("Deoxyschisandrin",     "五味子甲素"),
    ("Schisandrin B",        "五味子乙素"),
    ("Schisandrin C",        "五味子丙素"),
    ("Schisantherin A",      "五味子酯甲"),
    ("Schisantherin B",      "五味子酯乙"),
    ("Schisanhenol",         "五味子酚"),
    ("Gomisin A",            "戈米辛A"),
    ("Gomisin G",            "戈米辛G"),
    ("Gomisin J",            "戈米辛J"),
    ("Gomisin N",            "戈米辛N"),
    ("Anwulignan",           "安五脂素"),
    ("Pregomisin",           "前戈米辛"),
    ("Angeloylgomisin H",    "当归酰戈米辛H"),
    ("Tigloylgomisin P",     "巴豆酰戈米辛P"),
    ("Schisandrol A",        "五味子醇甲(别名)"),
    ("Benzoylgomisin Q",     "苯甲酰戈米辛Q"),
]

HERBS = [("Acanthopanax_senticosus", "刺五加", CIWUJIA),
         ("Schisandra_chinensis",    "五味子", WUWEIZI)]


# ---------------------------------------------------------------- ChEMBL 检索
def find_molecule(name):
    """按名检索 ChEMBL 分子，优先命中 pref_name / 同义词"""
    q = urllib.parse.quote(name)
    d = get("%s/molecule/search?q=%s&limit=10" % (BASE, q))
    if not d or not d.get("molecules"):
        return None
    low = name.lower()
    for m in d["molecules"]:
        pn = (m.get("pref_name") or "").lower()
        # 【坑】ChEMBL 的 molecule_synonyms 是 dict 列表（带 molecule_synonym /
        #      syn_type 字段），不是字符串列表，直接 join 会 TypeError。
        syn_parts = []
        for s in (m.get("molecule_synonyms") or []):
            if isinstance(s, dict):
                syn_parts.append(str(s.get("molecule_synonym") or ""))
            else:
                syn_parts.append(str(s))
        syn = " ".join(syn_parts).lower()
        if pn == low or low in pn or low in syn:
            return m
    return d["molecules"][0]


def molecule_props(chembl_id):
    d = get("%s/molecule/%s.json" % (BASE, chembl_id))
    if not d:
        return {}
    p = d.get("molecule_properties") or {}
    s = d.get("molecule_structures") or {}
    return {
        "smiles": s.get("canonical_smiles"),
        "mw": p.get("mw_freebase") or p.get("full_mwt"),
        "alogp": p.get("alogp"),
        "hba": p.get("hba"), "hbd": p.get("hbd"),
        "psa": p.get("psa"), "rtb": p.get("rtb"),
        "aromatic_rings": p.get("aromatic_rings"),
        "heavy_atoms": p.get("heavy_atoms"),
        "qed": p.get("qed_weighted"),
        "max_phase": d.get("max_phase"),
        "first_approval": d.get("first_approval"),
        "natural_product": d.get("natural_product"),
    }


def activities(chembl_id, max_n=1000):
    """人源靶点的活性记录（IC50/Ki/Kd/EC50 ≤ 10 µM）"""
    url = ("%s/activity.json?molecule_chembl_id=%s&target_organism=Homo%%20sapiens"
           "&limit=%d" % (BASE, chembl_id, max_n))
    rows = []
    while url and len(rows) < 4000:
        d = get(url)
        if not d or not d.get("activities"):
            break
        rows.extend(d["activities"])
        nxt = (d.get("page_meta") or {}).get("next")
        url = ("https://www.ebi.ac.uk" + nxt) if nxt else None
        if not nxt:
            break
        time.sleep(0.4)          # 【坑位46】分页之间限速，避免触发 429/500 限流
    out = []
    for a in rows:
        v = a.get("standard_value")
        t = a.get("standard_type")
        if v is None or t not in ("IC50", "Ki", "Kd", "EC50", "Potency", "AC50"):
            continue
        try:
            v = float(v)
        except (TypeError, ValueError):
            continue
        if v <= 10000:          # ≤ 10 µM
            out.append({"target_chembl_id": a.get("target_chembl_id"),
                        "type": t, "value_nM": v,
                        "pchembl": a.get("pchembl_value")})
    return out


_target_cache = {}


def targets_info_batch(tids, batch=50):
    """【坑】逐个查 target：槲皮素有几百个靶点，300+ 次请求要跑十几分钟。
    ChEMBL 的 target 列表接口支持 target_chembl_id__in 批量过滤，
    50 个/次 → 300 次请求压缩成 6 次。"""
    todo = [t for t in tids if t and t not in _target_cache]
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        q = urllib.parse.quote(",".join(chunk))
        url = "%s/target.json?target_chembl_id__in=%s&limit=%d" % (BASE, q, batch)
        while url:
            d = get(url)
            if not d:
                break
            for t in (d.get("targets") or []):
                info = {"pref_name": t.get("pref_name"),
                        "type": t.get("target_type"),
                        "organism": t.get("organism"),
                        "accessions": [c.get("accession")
                                       for c in (t.get("target_components") or [])
                                       if c.get("accession")]}
                _target_cache[t.get("target_chembl_id")] = info
            nxt = (d.get("page_meta") or {}).get("next")
            url = ("https://www.ebi.ac.uk" + nxt) if nxt else None
        for t in chunk:
            _target_cache.setdefault(t, None)
        time.sleep(0.3)
    return {t: _target_cache.get(t) for t in tids}


def target_info(tid):
    if tid not in _target_cache:
        targets_info_batch([tid])
    return _target_cache.get(tid)


_gene_cache = {}


def acc_to_symbol_batch(accs):
    """【坑】逐个 accession 查 mygene 太慢：槲皮素这类热门化合物有上千个靶点，
    一条一条请求要跑几十分钟。mygene.info 支持 POST 批量查询（最多 1000 条/次）。"""
    todo = [a for a in accs if a and a not in _gene_cache]
    for i in range(0, len(todo), 200):
        chunk = todo[i:i + 200]
        body = urllib.parse.urlencode({
            "q": ",".join(chunk), "scopes": "uniprot,accession",
            "fields": "symbol,query", "species": "human",
        }).encode("utf-8")
        try:
            req = urllib.request.Request("%s/query" % MYGENE, data=body, headers={
                "User-Agent": UA["User-Agent"],
                "Content-Type": "application/x-www-form-urlencoded",
            })
            with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
                res = json.loads(r.read().decode("utf-8"))
        except Exception as e:
            sys.stderr.write("  ! mygene batch failed: %s\n" % e)
            res = []
        for row in res:
            q = row.get("query")
            if q and row.get("symbol"):
                _gene_cache[q] = row["symbol"]
        for a in chunk:
            _gene_cache.setdefault(a, None)
        time.sleep(0.2)
    return {a: _gene_cache.get(a) for a in accs}


# ---------------------------------------------------------------- ADME 粗筛
def adme_flag(p):
    """Lipinski / Veber 粗筛（TCMSP 的 OB 拿不到时的替代口径）"""
    try:
        mw = float(p.get("mw") or 0); logp = float(p.get("alogp") or 0)
        hba = float(p.get("hba") or 0); hbd = float(p.get("hbd") or 0)
        rtb = float(p.get("rtb") or 0); psa = float(p.get("psa") or 0)
    except (TypeError, ValueError):
        return "NA", []
    viol = []
    if mw > 500: viol.append("MW>500")
    if logp > 5: viol.append("ALogP>5")
    if hba > 10: viol.append("HBA>10")
    if hbd > 5:  viol.append("HBD>5")
    if rtb > 10: viol.append("RotB>10")
    if psa > 140: viol.append("TPSA>140")
    return ("PASS" if len(viol) <= 1 else "FAIL"), viol


# ---------------------------------------------------------------- 主流程
def log(*a):
    # 重定向到文件时必须 flush，否则进度一直看不到
    import builtins
    builtins.print(*a, flush=True)


def main():
    outdir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          "results", "tables")
    os.makedirs(outdir, exist_ok=True)

    comp_rows, targ_rows = [], []

    # 【坑位46】断点续跑：若上次跑到一半被限流打断，已落盘的化合物直接跳过
    done_comp = set()
    for fn in (os.path.join(outdir, "s05_compounds.csv"),):
        if os.path.exists(fn):
            try:
                with open(fn, encoding="utf-8") as f:
                    for row in csv.DictReader(f):
                        if row.get("found") == "True" or row.get("found") is True:
                            done_comp.add(row.get("query"))
            except Exception:
                pass
    if done_comp:
        log("⟳ 断点续跑：已跳过 %d 个已完成成分" % len(done_comp))

    def flush():
        fc = os.path.join(outdir, "s05_compounds.csv")
        ft = os.path.join(outdir, "s05_compound_targets.csv")
        if comp_rows:
            with open(fc, "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=list(comp_rows[0].keys()))
                w.writeheader(); w.writerows(comp_rows)
        if targ_rows:
            with open(ft, "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=list(targ_rows[0].keys()))
                w.writeheader(); w.writerows(targ_rows)

    for herb_latin, herb_cn, comps in HERBS:
        log("\n########## %s (%s) ##########" % (herb_cn, herb_latin))
        for name, cn in comps:
            if name in done_comp:
                log("  ↷ 跳过 %s（已完成）" % name); continue
            time.sleep(1.2)        # 【坑位46】成分之间限速
            m = find_molecule(name)
            if not m:
                log("  ✗ %-22s (%s) 未命中" % (name, cn))
                comp_rows.append({"herb": herb_cn, "herb_latin": herb_latin,
                                  "query": name, "cn_name": cn,
                                  "chembl_id": None, "found": False})
                continue
            cid = m.get("molecule_chembl_id")
            p = molecule_props(cid) or {}
            ok, viol = adme_flag(p)

            acts = activities(cid)
            # 先批量取靶点信息，再批量换 gene symbol（两次批量，避免逐条请求）
            tids = sorted({a["target_chembl_id"] for a in acts})
            targets_info_batch(tids)

            need_acc = []
            for t in tids:
                ti = _target_cache.get(t)
                if ti and ti.get("organism") == "Homo sapiens":
                    need_acc.extend(ti["accessions"])
            acc2sym = acc_to_symbol_batch(sorted(set(need_acc)))

            syms = {}
            for a in acts:
                ti = _target_cache.get(a["target_chembl_id"])
                if not ti or ti.get("organism") != "Homo sapiens":
                    continue
                for acc in ti["accessions"]:
                    s = acc2sym.get(acc)
                    if not s:
                        continue
                    rec = syms.setdefault(s, {"best_nM": 1e18, "types": set(),
                                              "target_name": ti["pref_name"],
                                              "target_type": ti["type"]})
                    if a["value_nM"] < rec["best_nM"]:
                        rec["best_nM"] = a["value_nM"]
                    rec["types"].add(a["type"])

            for s, rec in syms.items():
                targ_rows.append({
                    "herb": herb_cn, "compound": name, "cn_name": cn,
                    "chembl_id": cid, "gene_symbol": s,
                    "best_nM": round(rec["best_nM"], 2),
                    "assay_types": ";".join(sorted(rec["types"])),
                    "target_name": rec["target_name"],
                    "target_type": rec["target_type"],
                })

            log("  ✓ %-22s %-14s MW=%-7s ALogP=%-5s 靶点=%-4d ADME=%s %s"
                  % (name, cid, str(p.get("mw"))[:7], str(p.get("alogp"))[:5],
                     len(syms), ok, ("(" + ",".join(viol) + ")") if viol else ""))

            comp_rows.append({
                "herb": herb_cn, "herb_latin": herb_latin, "query": name,
                "cn_name": cn, "chembl_id": cid, "found": True,
                "pref_name": m.get("pref_name"),
                "smiles": p.get("smiles"), "mw": p.get("mw"), "alogp": p.get("alogp"),
                "hba": p.get("hba"), "hbd": p.get("hbd"), "psa": p.get("psa"),
                "rtb": p.get("rtb"), "qed": p.get("qed"),
                "max_phase": p.get("max_phase"),
                "natural_product": p.get("natural_product"),
                "adme_pass": ok, "adme_violations": ";".join(viol),
                "n_targets": len(syms),
            })
            flush()        # 每完成一个成分就落盘，断点续跑不丢数据

    fc = os.path.join(outdir, "s05_compounds.csv")
    ft = os.path.join(outdir, "s05_compound_targets.csv")
    with open(fc, "w", newline="", encoding="utf-8") as f:
        if comp_rows:
            w = csv.DictWriter(f, fieldnames=list(comp_rows[0].keys()))
            w.writeheader(); w.writerows(comp_rows)
    with open(ft, "w", newline="", encoding="utf-8") as f:
        if targ_rows:
            w = csv.DictWriter(f, fieldnames=list(targ_rows[0].keys()))
            w.writeheader(); w.writerows(targ_rows)

    log("\n✓ 已保存：", fc)
    log("✓ 已保存：", ft)
    log("\n成分命中：%d / %d" % (sum(1 for r in comp_rows if r.get("found")), len(comp_rows)))
    log("化合物-靶点对：", len(targ_rows))
    allg = sorted(set(r["gene_symbol"] for r in targ_rows))
    log("去重靶点数：", len(allg))
    log("靶点：", ", ".join(allg[:60]))


if __name__ == "__main__":
    main()
