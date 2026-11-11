# -*- coding: utf-8 -*-
"""
fix_footer.py — 修复 Word 重开后页脚回归 instrText 的问题

【背景｜2026-11-11】
用户把稿件重命名为 manuscript11-11.docx（很可能用 Word 打开过）。
Word 保存时会把自己的 PAGE 域写成 <w:fldChar>/<w:instrText> 形式，
而 <w:instrText> **包含 <w:ins 子串** —— 这正是 2026-11-03 第一次
技术退回（PONE-D-26-54603）的根因，PLOS 检查器用 `if "<w:ins" in xml`
这种粗暴子串匹配，把干净的稿件误判成"含修订标记"。

【修法】
把页脚里的 PAGE 域改写为 <w:fldSimple w:instr=" PAGE ">
—— 属性形式，不含 instrText，Word 同样能正常渲染页码。

⚠ 修复后不要再用 Word 打开另存，否则会再次写回 instrText。
"""
import os, re, shutil, zipfile, sys

SUB = r'C:\Users\Fang\Desktop\S3\ciwujia_wuweizi_oa\PLOS_submission'
TARGET = os.path.join(SUB, sys.argv[1] if len(sys.argv) > 1 else 'manuscript11-11.docx')

FLD_RE = re.compile(
    r'<w:r>\s*<w:fldChar[^>]*w:fldCharType="begin"[^>]*/>\s*</w:r>'
    r'<w:r>\s*<w:instrText[^>]*>\s*PAGE\s*</w:instrText>\s*</w:r>'
    r'<w:r>\s*<w:fldChar[^>]*w:fldCharType="separate"[^>]*/>\s*</w:r>'
    r'<w:r>(?:<w:rPr>.*?</w:rPr>)?<w:t[^>]*>\s*\d*\s*</w:t></w:r>'
    r'<w:r>\s*<w:fldChar[^>]*w:fldCharType="end"[^>]*/>\s*</w:r>',
    re.S)

SIMPLE = '<w:fldSimple w:instr=" PAGE "><w:r><w:t>1</w:t></w:r></w:fldSimple>'


def clean_footer(xml):
    """把 fldChar/instrText 形式的 PAGE 域换成 fldSimple"""
    n_before = xml.count('instrText')
    xml = FLD_RE.sub(SIMPLE, xml)
    # 兜底：若还有残留 instrText（非 PAGE），说明结构不同，单独处理
    if 'instrText' in xml:
        # 通用做法：把整个 begin..end 域替换成 fldSimple
        xml = re.sub(
            r'<w:r>\s*<w:fldChar[^>]*w:fldCharType="begin"[^>]*/>\s*</w:r>.*?'
            r'<w:r>\s*<w:fldChar[^>]*w:fldCharType="end"[^>]*/>\s*</w:r>',
            SIMPLE, xml, flags=re.S)
    return xml, n_before - xml.count('instrText')


def main():
    src = zipfile.ZipFile(TARGET)
    parts = {n: src.read(n) for n in src.namelist()}
    infos = src.infolist()
    src.close()

    total_removed = 0
    for n in list(parts):
        if n.endswith('.xml') and ('footer' in n or 'header' in n):
            x = parts[n].decode('utf-8')
            if 'instrText' in x:
                x2, rm = clean_footer(x)
                parts[n] = x2.encode('utf-8')
                total_removed += rm
                print('  %s: instrText %d -> %d' % (n, rm, x2.count('instrText')))

    if total_removed == 0:
        print('页脚本来就是干净的，无需修改')
        return

    tmp = TARGET + '.tmp'
    with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as z:
        for info in infos:
            z.writestr(info, parts[info.filename])
    os.replace(tmp, TARGET)
    print('已修复，共移除 %d 处 instrText' % total_removed)

    # ---------- 验证 ----------
    z = zipfile.ZipFile(TARGET)
    bad = {}
    for n in z.namelist():
        if n.endswith('.xml'):
            x = z.read(n).decode('utf-8', errors='ignore')
            h = {k: x.count(k) for k in
                 ['<w:ins ', '<w:ins>', '<w:ins/', '<w:instrText', '<w:del ',
                  'rPrChange', 'pPrChange', 'w:comment', 'insideH', 'insideV']
                 if x.count(k)}
            if h:
                bad[n] = h
    print('\n最终 dirty 检查:', bad if bad else 'NONE  ✅')

    doc = z.read('word/document.xml').decode('utf-8')
    print('连续行号 lnNumType:', 'lnNumType' in doc)
    for n in z.namelist():
        if 'footer' in n and n.endswith('.xml'):
            f = z.read(n).decode('utf-8')
            print('  %s -> fldSimple:%d instrText:%d' % (n, f.count('fldSimple'), f.count('instrText')))


if __name__ == '__main__':
    main()
