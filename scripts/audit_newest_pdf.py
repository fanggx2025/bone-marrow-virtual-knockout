"""一键审计：自动在 Downloads 里找最新的 PLOS ONE 合并稿 PDF 并跑 audit_merged_pdf.py。

用法（下载完合并 PDF 之后直接跑）：
  "C:/Users/Fang/.workbuddy/binaries/python/envs/default/Scripts/python.exe" \
      "C:/Users/Fang/Desktop/S3/ciwujia_wuweizi_oa/scripts/audit_newest_pdf.py"

也可以手动指定文件：
  ... audit_newest_pdf.py "C:/Users/Fang/Downloads/某个文件.pdf"
"""
import os
import re
import subprocess
import sys

PY = sys.executable
HERE = os.path.dirname(os.path.abspath(__file__))
AUDIT = os.path.join(HERE, 'audit_merged_pdf.py')
DL = os.path.join(os.path.expanduser('~'), 'Downloads')

# 合并稿 PDF 的特征：文件名含 PONE / Submission / manuscript，且是较大体积的 PDF
# ⚠️ 只认 PLOS ONE 的稿件号 PONE-D-26-xxxxx：
#    PONE-S-26-xxxxx 是「补充材料/SI」编号，MD-S-26-xxxxx / CBAC-D-… 是别的刊，
#    它们文件名里同样含 PONE，绝不能被误当成合并稿。
# 排除上次的旧稿（用于"旧稿对照"的基线，不作为新稿审计目标）
EXCLUDE = {'pone-d-26-54603 (1).pdf', 'pone-d-26-54603.pdf'}
MANUSCRIPT_RE = re.compile(r'pone-d-\d{2}-\d{5}', re.I)


def find_candidates():
    if not os.path.isdir(DL):
        return []
    out = []
    for fn in os.listdir(DL):
        if not fn.lower().endswith('.pdf'):
            continue
        if fn.lower() in EXCLUDE:
            continue
        full = os.path.join(DL, fn)
        try:
            size = os.path.getsize(full)
        except OSError:
            continue
        if size < 200 * 1024:      # 合并稿至少几百 KB
            continue
        # 必须文件名里带 PONE-D-26-xxxxx 稿件号，且不是上次的旧稿
        if not MANUSCRIPT_RE.search(fn):
            continue
        if fn.lower() in EXCLUDE:
            continue
        out.append((os.path.getmtime(full), full, size))
    out.sort(reverse=True)
    return out


def main():
    if len(sys.argv) > 1:
        target = sys.argv[1]
        if not os.path.isfile(target):
            print('指定的文件不存在：%s' % target)
            return 2
    else:
        cands = find_candidates()
        if not cands:
            print('在 %s 里没找到像样的合并稿 PDF。' % DL)
            print('请先在 EM 里 Build PDF for Approval 并下载，或把路径作为参数传进来。')
            return 1
        print('在 Downloads 里找到 %d 个候选，最新的是：' % len(cands))
        for mt, full, size in cands[:5]:
            import datetime
            print('   %s  %7.1f KB  %s' % (
                datetime.datetime.fromtimestamp(mt).strftime('%Y-%m-%d %H:%M'),
                size / 1024.0, os.path.basename(full)))
        print('')
        target = cands[0][1]

    print('审计目标：%s' % target)
    print('=' * 104)
    r = subprocess.run([PY, AUDIT, target])
    return r.returncode


if __name__ == '__main__':
    sys.exit(main())
