# -*- coding: utf-8 -*-
"""Build results/Supplementary_Information.pdf:
   cover + 10 supplementary figures (vector-embedded, each with caption)
   + 15 supplementary tables (S1-S15, landscape A4).
"""
import os, csv, tempfile, re
import fitz  # PyMuPDF
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_LEFT
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                TableStyle, PageBreak)
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase.cidfonts import UnicodeCIDFont

ROOT = "C:/Users/Fang/Desktop/S3/ciwujia_wuweizi_oa"
FIG = os.path.join(ROOT, "results/figures")
TAB = os.path.join(ROOT, "results/tables")
OUT = os.path.join(ROOT, "results/Supplementary_Information.pdf")
ARIAL = "C:/Windows/Fonts/arial.ttf"
ARIALB = "C:/Windows/Fonts/arialbd.ttf"
CJK = re.compile(r"[\u4e00-\u9fff]")

# ---------------------------------------------------------------------------
FIGS = [
    ("S1", "FigS1_dualmethod_convergence.pdf",
     "Dual-method docking convergence. (a) local Vina \u0394G; (b) SwissDock AC/EADock SP-dG; "
     "(c) pose-into-pocket distances for both methods, all confirming 8/8 occupancy of the "
     "orthosteric pocket (pose-centroid distance to chain-A BRL 2.8\u20134.8 \u00c5)."),
    ("S2", "Fig_S_Kaempferol.pdf",
     "PPARG (2PRG chain-A ligand-binding domain) best Vina pose overlaid on the co-crystallized "
     "rosiglitazone (BRL). Compound: Kaempferol (Eleutherococcus senticosus)."),
    ("S3", "Fig_S_Schisandrin_B.pdf",
     "PPARG 2PRG chain-A best Vina pose overlaid on BRL. Compound: Schisandrin B (Schisandra chinensis)."),
    ("S4", "Fig_S_Schisantherin_B.pdf",
     "PPARG 2PRG chain-A best Vina pose overlaid on BRL. Compound: Schisantherin B (Schisandra chinensis)."),
    ("S5", "Fig_S_Deoxyschisandrin.pdf",
     "PPARG 2PRG chain-A best Vina pose overlaid on BRL. Compound: Deoxyschisandrin (Schisandra chinensis)."),
    ("S6", "Fig_S_Schisanhenol.pdf",
     "PPARG 2PRG chain-A best Vina pose overlaid on BRL. Compound: Schisanhenol (Schisandra chinensis)."),
    ("S7", "Fig_S_Schisandrol_B.pdf",
     "PPARG 2PRG chain-A best Vina pose overlaid on BRL. Compound: Schisandrol B (Schisandra chinensis)."),
    ("S8", "Fig_S_Gomisin_J.pdf",
     "PPARG 2PRG chain-A best Vina pose overlaid on BRL. Compound: Gomisin J (Schisandra chinensis)."),
    ("S9", "Fig_S_Angeloylgomisin_H.pdf",
     "PPARG 2PRG chain-A best Vina pose overlaid on BRL. Compound: Angeloylgomisin H (Schisandra chinensis)."),
    ("S10", "Fig_S_wuweizi_PPARG_convergence.pdf",
     "Wuweizi PPARG convergence bar chart \u2014 Vina \u0394G of the eight compounds stratified across "
     "donors / virtual-knockout signal."),
]

TABLES = [
    ("S1",  "s02_subtype_assignment.csv",
     "Stromal-cell subtype assignment per cell (5 subclusters; n = 10,206)."),
    ("S2",  "s03_gene_pseudotime_corr.csv",
     "Gene\u2013pseudotime Spearman correlation (\u03c1) for the curated osteogenic gene set."),
    ("S3",  "s04_candidate_genes.csv",
     "Curated + data-driven candidate regulator list with detection and pseudotime correlation."),
    ("S4",  "s04_knk_scores.csv",
     "Virtual-knockout scores (med_w_end, p_perturb, p_repress, pct_rank_end) for all 183 regulators."),
    ("S5",  "s05_compounds.csv",
     "36 herb compounds (18 Eleutherococcus senticosus + 18 Schisandra chinensis) with ChEMBL IDs / SMILES."),
    ("S6",  "s05_compound_targets.csv",
     "314 ChEMBL compound\u2013target pairs spanning 196 unique targets."),
    ("S7",  "s05b_wuweizi_literature_targets.csv",
     "Schisandra chinensis lignan NP supplement (PPARG/AKT1; source = Literature_NP)."),
    ("S8",  "s06_consensus_targets.csv",
     "Three-way consensus target set with tier labels and per-layer evidence."),
    ("S9",  "s07_oa_validation.csv",
     "OA validation \u2014 log2FC and FDR per consensus target across stromal fibroblast / osteoblast / CAR cell."),
    ("S10", "s09_swissdock_contacts.csv",
     "SwissDock contacts for the kaempferol global EADock pose."),
    ("S11", "s11_webdock_lignans.csv",
     "SwissDock AC/EADock per-lignan results: SP-dG, best cluster, pose-centroid distance to BRL, key contacts."),
    ("S12", "s12_dualmethod_summary.csv",
     "Vina + AC dual-method summary with pocket-occupancy agreement."),
    ("S13", "s10_wuweizi_vina.csv",
     "Eight-compound PPARG docking under the chain-A protocol (Vina \u0394G, pocket occupancy and key-residue contacts)."),
    ("S14", "s14_consensus_calibration.csv",
     "Null model and calibration of the three evidence layers: layer sizes, exact "
     "hypergeometric null for the Tier-1 intersection, permutation sensitivity "
     "analysis, stratified trajectory-concordance test and compound-layer enrichment test."),
    ("S15", "s14_control_gene_directions.csv",
     "Control-gene direction check: seven control regulators with expected and "
     "observed knockout direction."),
]

# ---------------------------------------------------------------------------
# reportlab styles (Arial TTF for full Latin/Greek coverage; CID STSong for CJK cells)
# ---------------------------------------------------------------------------
pdfmetrics.registerFont(TTFont("Arial", ARIAL))
pdfmetrics.registerFont(TTFont("Arial-Bold", ARIALB))
pdfmetrics.registerFont(UnicodeCIDFont("STSong-Light"))
_cell = ParagraphStyle("cell", fontName="Arial", fontSize=6.5, leading=7.6, alignment=TA_LEFT)
_cellcn = ParagraphStyle("cellcn", fontName="STSong-Light", fontSize=6.5, leading=7.6, alignment=TA_LEFT)
_hcell = ParagraphStyle("hcell", fontName="Arial-Bold", fontSize=6.5, leading=7.6, alignment=TA_LEFT)
_title = ParagraphStyle("title", fontName="Arial-Bold", fontSize=11, leading=13, spaceAfter=2)
_desc = ParagraphStyle("desc", fontName="Arial", fontSize=8.5, leading=10.5, spaceAfter=4)

def cell_style(v):
    return _cellcn if CJK.search(v) else _cell

def render_tables_pdf(path):
    doc = SimpleDocTemplate(path, pagesize=landscape(A4),
                            leftMargin=1.4*cm, rightMargin=1.4*cm,
                            topMargin=1.3*cm, bottomMargin=1.3*cm,
                            title="Supplementary Tables S1-S15")
    usable_w = landscape(A4)[0] - 2.8*cm
    story = [Paragraph("Supplementary Tables S1\u2013S15", _title), Spacer(1, 4)]
    first = True
    for sid, fn, desc in TABLES:
        if not first:
            story.append(PageBreak())
        first = False
        story.append(Paragraph("Supplementary Table %s." % sid, _title))
        story.append(Paragraph(desc, _desc))
        with open(os.path.join(TAB, fn), encoding="utf-8", errors="replace", newline="") as f:
            rows = list(csv.reader(f))
        header, data = rows[0], rows[1:]
        ncol = len(header)
        cw = usable_w / ncol
        cell_rows = [[Paragraph(str(h), _hcell) for h in header]]
        for r in data:
            r = (r + [""] * ncol)[:ncol]
            cell_rows.append([cell_style(str(v).replace("\n", " ")).clone("c", text=str(v))
                              if False else Paragraph(str(v).replace("\n", " "), cell_style(str(v)))
                              for v in r])
        t = Table(cell_rows, colWidths=[cw] * ncol, repeatRows=1)
        t.setStyle(TableStyle([
            ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#9aa0a6")),
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e8eaed")),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("LEFTPADDING", (0, 0), (-1, -1), 2.5),
            ("RIGHTPADDING", (0, 0), (-1, -1), 2.5),
            ("TOPPADDING", (0, 0), (-1, -1), 1.6),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 1.6),
        ]))
        story.append(t)
    doc.build(story)

# ---------------------------------------------------------------------------
def wrap_lines(text, font, size, width):
    words = text.split()
    lines, cur = [], ""
    for w in words:
        trial = (cur + " " + w).strip()
        if font.text_length(trial, fontsize=size) <= width or not cur:
            cur = trial
        else:
            lines.append(cur); cur = w
    if cur:
        lines.append(cur)
    return lines

def build():
    for fp in (ARIAL, ARIALB):
        assert os.path.exists(fp), "missing font " + fp
    tmp = tempfile.mktemp(suffix=".pdf")
    render_tables_pdf(tmp)

    A4w, A4h = fitz.paper_size("a4")
    f_reg = fitz.Font(fontfile=ARIAL)
    f_bold = fitz.Font(fontfile=ARIALB)
    doc = fitz.open()

    # ---- cover ----
    p = doc.new_page(width=A4w, height=A4h)
    p.insert_font(fontname="AR", fontfile=ARIAL)
    p.insert_font(fontname="ARB", fontfile=ARIALB)
    blocks = [
        ("Supplementary Information", 22, "ARB", 0.28, 60),
        ("Companion to the manuscript \u201cSingle-cell virtual knockout identifies PPARG, "
         "HIF1A and RUNX1 as negative regulators of human osteogenesis\u201d", 11, "AR", 0.40, 90),
        ("Target journal: PLOS ONE", 10, "AR", 0.52, 30),
        ("Contents: 10 supplementary figures (S1\u2013S10) and 15 supplementary tables (S1\u2013S15).",
         9.5, "AR", 0.58, 40),
    ]
    y = A4h * 0.26
    for txt, size, fname, _rel, boxh in blocks:
        rect = fitz.Rect(60, y, A4w - 60, y + boxh)
        rc = p.insert_textbox(rect, txt, fontsize=size, fontname=fname,
                              align=fitz.TEXT_ALIGN_CENTER, lineheight=1.4)
        assert rc >= 0, "cover text overflow: " + txt[:40]
        y += boxh + 14

    # ---- figure pages: heading, caption, figure below (vector) ----
    for sid, fn, cap in FIGS:
        src = fitz.open(os.path.join(FIG, fn))
        sp = src[0]
        sw, sh = sp.rect.width, sp.rect.height
        p = doc.new_page(width=A4w, height=A4h)
        p.insert_font(fontname="AR", fontfile=ARIAL)
        p.insert_font(fontname="ARB", fontfile=ARIALB)
        y = 46
        p.insert_text((42, y), "Supplementary Fig. %s." % sid, fontsize=11.5, fontname="ARB")
        y += 10
        lines = wrap_lines(cap, f_reg, 9.5, A4w - 84)
        for ln in lines:
            y += 13.5
            p.insert_text((42, y), ln, fontsize=9.5, fontname="AR")
        area = fitz.Rect(42, y + 16, A4w - 42, A4h - 64)
        if sw / sh > area.width / area.height:
            w = area.width; h = w * sh / sw
        else:
            h = area.height; w = h * sw / sh
        x0 = area.x0 + (area.width - w) / 2
        y0 = area.y0
        p.show_pdf_page(fitz.Rect(x0, y0, x0 + w, y0 + h), src, 0)
        src.close()

    # ---- tables ----
    doc.insert_pdf(fitz.open(tmp))

    n = doc.page_count
    for i in range(n):
        doc[i].insert_textbox((A4w/2 - 160, A4h - 34, A4w/2 + 160, A4h - 18),
                              "Supplementary Information  \u00b7  page %d" % (i + 1),
                              fontsize=8, fontname="helv", align=fitz.TEXT_ALIGN_CENTER)
    doc.save(OUT, deflate=True, clean=True)
    doc.close()
    os.remove(tmp)
    print("WROTE", OUT, "pages=", n)

if __name__ == "__main__":
    build()
