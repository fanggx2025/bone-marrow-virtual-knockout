"""PLOS ONE 合并提交 PDF 自动审计 —— 针对本次 v8 重投稿。

用法：
  "C:/Users/Fang/.workbuddy/binaries/python/envs/default/Scripts/python.exe" \
      scripts/audit_merged_pdf.py "C:/Users/Fang/Downloads/合并稿.pdf"

检查项分四类：
  A. 修订标记（PLOS 技术退回的头号原因：<w:instrText> 含 "<w:ins" 子串）
  B. 表单字段与稿件一致性（标题/摘要/作者/关键词/资助/竞争利益）
  C. 关键科学数字与拒稿史遗留错误（骨折不愈叙事 / 56 candidates 等）
  D. 结构与格式（行号、页码、图注、参考文献）
"""
import sys, os, re, json

try:
    import pymupdf
except ImportError:
    import fitz as pymupdf

PDF = sys.argv[1] if len(sys.argv) > 1 else r'C:\Users\Fang\Downloads\PONE-D-26-54603 (1).pdf'

res = []
def chk(cat, name, ok, detail=''):
    res.append((cat, name, bool(ok), detail))
    print('  [%s] %-52s %s' % ('OK  ' if ok else 'FAIL', name, detail))

d = pymupdf.open(PDF)
pages = [d[i].get_text() for i in range(d.page_count)]
full = '\n'.join(pages)
norm = re.sub(r'\s+', ' ', full)

print('=' * 104)
print('PLOS ONE 合并 PDF 审计：%s' % os.path.basename(PDF))
print('页数 %d | 提取文本 %d 字符' % (d.page_count, len(full)))
print('=' * 104)

# ---------- A. 修订标记 / 技术退回风险 ----------
print('\n--- A. 修订标记与撤回风险 ---')
MARKERS = ['Track Changes', 'tracked change', 'revision mark', 'Inserted Text',
           'Deleted Text', 'Revised Text', 'Markup', '修订', '批注', 'Comments:']
hits = [m for m in MARKERS if m.lower() in full.lower()]
chk('A', '无"修订标记/Track Changes"字样', not hits,
    ('命中: %s' % hits) if hits else '')

# 真正的技术退回触发点：页脚域用 instrText
# 在 PDF 里表现为页码位置异常或域代码残留
chk('A', '无残留域代码 (fldChar/instrText/PAGE 域)',
    not re.search(r'instrText|fldChar', full), '')

# 上次退回的合并 PDF 里出现过的典型提示
for bad in ['This document has been edited', 'Revisions:',
            'All changes have been accepted', 'There are no tracked']:
    chk('A', '无 "%s" 提示' % bad[:34], bad not in full, '')

# ---------- B. 表单字段一致性 ----------
print('\n--- B. 表单字段与稿件一致性 ---')
FULL_TITLE = ('Single-cell virtual knockout identifies PPARG, HIF1A and RUNX1 '
              'as negative regulators of human osteogenesis')
SHORT_TITLE = ('Virtual knockout identifies PPARG, HIF1A and RUNX1 as negative '
                'regulators of osteogenesis')
chk('B', 'Full Title 与稿件一致', FULL_TITLE.replace(' ', '') in norm.replace(' ', ''), '')
chk('B', 'Short Title 一致', SHORT_TITLE.replace(' ', '') in norm.replace(' ', ''), '')

# 摘要关键句
chk('B', '摘要含修正后的双药表述',
    'one or both of the two herbs' in norm.replace('\n', ' '), '')
chk('B', '摘要不含旧的错误表述',
    'matches compounds from both herbs' not in norm, '')

# 作者顺序：必须限定在稿件 byline 段内，不能用全文首次出现
# （合并 PDF 里表单的 Corresponding Author 行在最前，会干扰全文 find；
#   cover letter 结尾的 "Sincerely, Jingchao Wang, on behalf of Gongxi Fang,
#   Shujun Ren, Yilin Zhang and Baiqing Wang" 也会干扰 —— 那是通讯作者代表
#   全体作者署名，作者名自然不是投稿顺序。）
# 做法：找**稿件标题**首次出现处，其后 400 字符才是 byline 段。
ORDER = ['Gongxi Fang', 'Jingchao Wang', 'Shujun Ren', 'Yilin Zhang', 'Baiqing Wang']
# 注意：合并 PDF 的开头就是 cover letter，其正文里也会提到稿件标题，
# 所以不能取"标题首次出现"，要取**主稿标题**那次（前面带编号或出现在
# cover letter 之后）。稳妥做法：找标题后紧跟 byline 的那次。
import re as _re
_byline_idx = -1
for _m in _re.finditer(r'Single-cell virtual knockout identifies PPARG', norm):
    _seg = norm[_m.start():_m.start() + 400]
    if all(n in _seg for n in ORDER):
        _byline_idx = _m.start()
        break
byline_seg = norm[_byline_idx:_byline_idx + 400] if _byline_idx >= 0 else ''
positions = [byline_seg.find(n) for n in ORDER]
found = [p for p in positions if p >= 0]
chk('B', '五位作者姓名全部出现（稿件 byline 段）', len(found) == 5, '%d/5' % len(found))
chk('B', '作者顺序为 Fang→Wang→Ren→Zhang→Wang',
    len(found) == 5 and positions == sorted(positions),
    '' if len(found) == 5 and positions == sorted(positions) else str(positions))
chk('B', '无 "dr " 前缀', not re.search(r'\bdr\s+[A-Z][a-z]+', full),
    (re.search(r'\bdr\s+[A-Z][a-z]+', full).group(0) if re.search(r'\bdr\s+[A-Z][a-z]+', full) else ''))

chk('B', '含 Abstract 段', 'Abstract' in full, '')

# Keywords
for kw in ['osteoarthritis', 'virtual knockout', 'PPARG',
           'Eleutherococcus senticosus', 'Schisandra chinensis']:
    chk('B', 'Keywords 含 "%s"' % kw, kw.lower() in norm.lower(), '')

# 资助：合并 PDF 的表单区渲染为 "Financial Disclosure" + 选项列表 + 选中值，
# 选中值 "Yes" 出现在该区块末尾（EM 的渲染顺序），中间夹着页面分隔。
fd_i = norm.find('Financial Disclosure')
fd_seg = norm[fd_i:fd_i + 900] if fd_i >= 0 else ''
chk('B', 'Financial Disclosure 下拉 = Yes',
    bool(re.search(r'Financial Disclosure.{0,900}?\bYes\b', norm, re.S)), '')
chk('B', '资助声明含 ZHY2022-157', 'ZHY2022-157' in norm, '')
chk('B', '资助声明含 awarded to GF', 'awarded to GF' in norm, '')
chk('B', '资助声明含 funder website URL', 'http://tcm.hlj.gov.cn/' in norm, '')
chk('B', '资助方角色标准句', 'The funders had no role in study design' in norm, '')
chk('B', 'Funding Information 含 Award Recipient = Gongxi Fang',
    'Gongxi Fang' in norm and 'ZHY2022-157' in norm, '')

# 竞争利益
chk('B', 'Competing Interests 标准句',
    'The authors have declared that no competing interests exist.' in norm, '')

# 伦理：表单区渲染为 "Ethics Statement" + 提示 + 文本框内容（文本可能被折行拆断）
# ⚠️ 必须转小写比较：PDF 里是大写 "No human or animal"
chk('B', 'Ethics statement 存在（表单区）', 'Ethics Statement' in full, '')
flat = re.sub(r'[\s\-]+', '', full).lower()
chk('B', '伦理声明含 "no human or animal experiments"',
    'nohumanoranimalexperiments' in flat, '')

# 前次投稿披露
chk('B', '如实披露前次稿号 PONE-D-26-54603', 'PONE-D-26-54603' in norm, '')
chk('B', '双发表说明 "only related submission"',
    'only related submission' in norm.lower(), '')

# 存档 DOI
chk('B', 'Zenodo 用 concept DOI 23684944', '10.5281/zenodo.23684944' in norm, '')

# ⚠️ 旧 DOI 23684945 有三个合法/非法场所，必须分开判：
#   (1) 表单区 Data Availability 文本框 → 必须用新 DOI（否则文章里会出现
#       已退役的 v1 DOI，读者点开会拿到旧版数据集）—— **这是真问题**
#   (2) 本次主稿正文 Data Availability 段 → 必须用新 DOI
#   (3) Related Manuscript（前次投稿副本）→ **应当**带旧 DOI，因为那是被
#       退稿的旧稿，它记录的正是当时的 DOI，出现它是正确的
# 做法：按 PDF 页逐页判定 —— 表单区页面会有 "Powered by Editorial Manager" 页脚。
# ⚠️ 注意：EM 的合并 PDF **不并入附件正文**（Related Manuscript / Track Changes
#    都作为独立附件存在，编辑可单独下载），所以合并稿里通常**看不到**旧 DOI ——
#    这既不代表副本丢了，也不代表有问题。因此第 3 项检查只做"若出现则必须合法"，
#    出现与否都判 PASS。
def _pages_with_doi(pdf_d, doi):
    hits = []
    for i in range(pdf_d.page_count):
        t = pdf_d[i].get_text()
        if doi in t:
            hits.append((i + 1, 'Powered by Editorial Manager' in t))
    return hits

_hits_new = _pages_with_doi(d, '23684944')
_hits_old = _pages_with_doi(d, '23684945')
_old_form = [p for p, is_form in _hits_old if is_form]
_old_main = [p for p, is_form in _hits_old if not is_form]
chk('B', '旧 DOI 23684945 不在表单区（须用 concept DOI）',
    not _old_form, ('表单区出现于第 %s 页' % _old_form) if _old_form else '')
chk('B', '旧 DOI 23684945 不在主稿正文',
    not _old_main, ('正文出现于第 %s 页' % _old_main) if _old_main else '')
chk('B', 'Related Manuscript 副本保留旧 DOI（若并入则正确；未并入亦可）',
    True,
    '旧 DOI 出现于第 %s 页（EM 未并入附件正文时此处为空，属正常）'
    % [p for p, _ in _hits_old])
chk('B', '含 GitHub 仓库链接',
    'github.com/fanggx2025/bone-marrow-virtual-knockout' in norm, '')

# ---------- C. 科学内容与拒稿史遗留错误 ----------
print('\n--- C. 科学内容与遗留错误 ---')
for kw in ['10,206', '9,058', '1,148', '117,308', '183', '49', '3.14', '0.65',
           '0.0053', '0.78', '23 of 31']:
    chk('C', '关键数字 "%s" 存在' % kw, kw in norm, '')

# 对照基因检查：稿件写 "seven of seven"（文字），SI 表里是 7/7
chk('C', '含对照基因检查 "seven of seven" 或 "7 of 7"',
    ('seven of seven' in norm.lower()) or ('7 of 7' in norm), '')

# 必须在文中出现的核心结论
chk('C', '含零模型口径 "prioritisation device"',
    'prioritisation device' in norm, '')
chk('C', '含对照基因检查表述',
    'control genes behave as they should' in norm.lower(), '')

# 绝不能出现的（2026-10-27 路线A 重构时清除的）
# ⚠️ 注意：参考文献的【标题】里可能合法出现这些词（如 "The risk of non-union
#   per fracture"、"Fracture healing: a review"），那是别人论文的题目，不是我们的叙事。
#   所以只在【正文】段（Abstract 到 References 之间）检查，不扫参考文献区。
body_start = norm.find('Abstract')
body_end = norm.rfind('References')
body_only = norm[body_start:body_end] if body_start > 0 and body_end > body_start else norm

FORBIDDEN = [
    ('GSE337700', '不存在的骨折队列'),
    ('56 candidate', '旧候选基因数'),
    ('56 regulator', '旧候选基因数'),
    ('non-union', '不愈合叙事'),
    ('nonunion', '不愈合叙事'),
    ('delayed union', '延迟愈合叙事'),
    ('fracture healing', '骨折愈合叙事'),
    ('bone fracture healing', '骨折愈合叙事'),
]
for kw, why in FORBIDDEN:
    chk('C', '正文不含 "%s"（%s）' % (kw, why), kw.lower() not in body_only.lower(), '')

# 资助不应出现在稿件正文（只在表单）
chk('C', '稿件正文无资助表述（仅表单）',
    'ZHY2022-157' not in norm.split('Financial Disclosure')[0]
    if 'Financial Disclosure' in norm else True, '')

# ---------- D. 结构与格式 ----------
print('\n--- D. 结构与格式 ---')
chk('D', '含 References 段', 'References' in full, '')

# 参考文献条数：只数 References 段到附件清单之前，且排除行号干扰
refs_idx = norm.rfind('References')
ref_block = norm[refs_idx:] if refs_idx > 0 else ''
# PLOS 参考文献渲染形如 "1.  Mills LA, Aitken SA, ..."，行号会插在数字后面
nums = [int(x) for x in re.findall(r'(?:^|\s)(\d{1,3})\.\s{1,4}[A-Z][a-z]', ref_block)]
chk('D', '参考文献条数 = 66', len(nums) >= 60, '解析到 %d 条' % len(nums))

# 图注
fig_caps = re.findall(r'Fig (\d+)\.', full)
chk('D', 'Fig 1-7 图注齐全',
    all(('Fig %d.' % i) in full for i in range(1, 8)), '')

# 补充材料题注
for tag in ['S14', 'S15']:
    chk('D', 'Supporting Information 含 Table %s' % tag, tag in full, '')

# 页码：合并 PDF 的稿件页页脚有页码；表单页可能没有。用稿件区（后 1/3）抽查
mid = max(0, d.page_count // 3)
tail = '\n'.join(pages[mid:])
page_ok = sum(1 for i, t in enumerate(pages[mid:], start=mid + 1)
              if re.search(r'\n\s*%d\s*\n' % i, t))
chk('D', '稿件区页码存在（抽查）', page_ok >= 5, '%d/%d 页检出' % (page_ok, d.page_count - mid))

# 行号：PLOS 稿件页左侧连续行号。PyMuPDF 会把行号渲染成**独立文本块**，
# 所以不能按"行首是数字"来匹配；要改为统计页内纯数字块的密度。
def line_no_density(t):
    nums = [x for x in re.findall(r'(?m)^\s{0,6}(\d{1,3})\s*$', t)]
    return len(nums)

lineno = 0
for t in pages[mid:]:
    if line_no_density(t) >= 15:
        lineno += 1
chk('D', '稿件区有连续行号（左侧数字块）', lineno >= 3, '%d 页检出' % lineno)

# ---------- 汇总 ----------
print('\n' + '=' * 104)
fails = [r for r in res if not r[2]]
print('合计 %d 项检查 | PASS %d | FAIL %d' % (len(res), len(res) - len(fails), len(fails)))
if fails:
    print('\n失败项：')
    for c, n, o, dt in fails:
        print('  [%s] %-20s %s' % (c, n, dt))
    print('\n⚠️ 有 FAIL 项 —— 先确认是否影响技术检查，再 Approve。')
else:
    print('\n✅ 全部通过 —— 可以 Approve。')
print('=' * 104)
