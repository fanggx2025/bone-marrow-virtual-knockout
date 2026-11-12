"""生成 PLOS ONE 要求的 Track Changes 版稿件（Word Compare）。

用途：PONE-D-26-56816 被退回，编辑要求上传一份标出与前次投稿
      PONE-D-26-54603 差异的 Track Changes 版本（item type = 'Other'）。

原理：Word 的 CompareDocuments（原文档 vs 修订文档）会生成真正的
      <w:ins>/<w:del> 修订标记，不是伪装的"加粗+着色"。

关键风险与对策：
  ⚠️ 坑位58：Word 保存会把页脚 PAGE 域写回 <w:instrText>，而 PLOS 检查器
     `if "<w:ins" in xml` 会因 "<w:ins" 子串误判为"含修订标记"→ 技术退回。
     → 生成后必须重跑 101_fix_word_footer.py 把 PAGE 域改回 <w:fldSimple>。
     （注意：本文件是**故意**含修订标记的，所以 PLOS 的人工/内容检查会看到
       真修订 —— 这是编辑要的；但页脚域代码仍需清理，避免误判。）

用法：
  python scripts/102_make_track_changes.py
"""
import os
import shutil
import subprocess
import sys
import time

import win32com.client as win32

BASE = r'C:\Users\Fang\Desktop\S3\ciwujia_wuweizi_oa'
SUB = os.path.join(BASE, 'PLOS_submission')
SCRIPTS = os.path.join(BASE, 'scripts')

# 基线 = PLOS 前次投稿 PONE-D-26-54603 的正文原样（已实测：与 Related
#        Manuscript 副本逐段完全一致，副本仅多 1 段说明）
ORIGINAL = os.path.join(SUB, 'manuscript11-03.docx')
# 修订 = 本次投的 clean 主稿
REVISED = os.path.join(SUB, 'manuscript11-11.docx')
OUTPUT = os.path.join(SUB, 'manuscript11-11_trackchanges.docx')

# Word 常量
wdCompareTargetNew = 2          # wdCompareDestinationNew
wdGranularityWordLevel = 1      # wdGranularityWordLevel
wdCompareDestinationNew = 2
wdDoNotSaveChanges = 0
wdFormatXMLDocument = 12


def main():
    for p in (ORIGINAL, REVISED):
        if not os.path.isfile(p):
            print('找不到文件：%s' % p)
            return 2

    # 先备份输出文件
    if os.path.exists(OUTPUT):
        shutil.copy2(OUTPUT, OUTPUT + '.bak')

    word = win32.gencache.EnsureDispatch('Word.Application')
    word.Visible = False
    word.DisplayAlerts = 0
    doc_orig = doc_rev = comp = None
    try:
        doc_orig = word.Documents.Open(ORIGINAL, ReadOnly=True, AddToRecentFiles=False)
        doc_rev = word.Documents.Open(REVISED, ReadOnly=True, AddToRecentFiles=False)

        # CompareDocuments 参数顺序（Word 2010+）：
        # OriginalDocument, RevisedDocument, Destination, Granularity,
        # CompareFormatting, CompareCaseChanges, CompareWhitespace,
        # CompareTables, CompareHeaders, CompareFootnotes, CompareTextboxes,
        # CompareFields, CompareComments, CompareMoves, RevisedAuthor,
        # IgnoreAllComparisonWarnings
        #
        # ⚠️ CompareFormatting=False 是关键：若为 True，Word 会把两稿之间的
        #    纯格式差异（字体、字号、删除线、底纹…）也标成修订，产生成百上千
        #    条"设置了格式: 字体: 删除线"之类的噪音气泡，编辑根本看不过来，
        #    反而掩盖真正重要的文字改动。
        #    关闭后只保留**文字层面**的增删，正是编辑想看的。
        comp = word.CompareDocuments(
            doc_orig, doc_rev,
            wdCompareTargetNew,
            wdGranularityWordLevel,
            False,  # CompareFormatting —— 关掉格式噪音
            True,   # CompareCaseChanges
            True,   # CompareWhitespace
            True,   # CompareTables
            True,   # CompareHeaders
            True,   # CompareFootnotes
            True,   # CompareTextboxes
            True,   # CompareFields
            True,   # CompareComments
            True,   # CompareMoves
            'Author',  # RevisedAuthor
        )

        # 确保修订可见（默认即为显示）
        try:
            comp.ActiveWindow.View.ShowRevisionsAndComments = True
            comp.ActiveWindow.View.RevisionsView = 0   # wdRevisionsViewFinal
        except Exception:
            pass

        comp.SaveAs2(OUTPUT, FileFormat=wdFormatXMLDocument,
                     AddToRecentFiles=False)
        comp.Close(wdDoNotSaveChanges)
        comp = None
        doc_orig.Close(wdDoNotSaveChanges)
        doc_rev.Close(wdDoNotSaveChanges)
        doc_orig = doc_rev = None
        print('已生成：%s' % OUTPUT)
    finally:
        try:
            if comp is not None:
                comp.Close(wdDoNotSaveChanges)
        except Exception:
            pass
        try:
            if doc_orig is not None:
                doc_orig.Close(wdDoNotSaveChanges)
        except Exception:
            pass
        try:
            if doc_rev is not None:
                doc_rev.Close(wdDoNotSaveChanges)
        except Exception:
            pass
        word.Quit()

    # Word 退出后文件句柄可能还没释放，稍等再改写（否则 PermissionError）
    for _ in range(20):
        time.sleep(0.5)
        try:
            with open(OUTPUT, 'r+b'):
                break
        except PermissionError:
            continue
    tmp_leftover = OUTPUT + '.tmp'
    if os.path.exists(tmp_leftover):
        try:
            os.remove(tmp_leftover)
        except OSError:
            pass

    # ---- 关键：修页脚域代码（坑位58）----
    fixer = os.path.join(SCRIPTS, '101_fix_word_footer.py')
    if os.path.isfile(fixer):
        print('\n--- 运行 101_fix_word_footer.py 清理页脚域代码 ---')
        r = subprocess.run([sys.executable, fixer, os.path.basename(OUTPUT)],
                           capture_output=True, text=True, cwd=SUB)
        print(r.stdout.strip())
        if r.stderr.strip():
            print('STDERR:', r.stderr.strip()[:600])
    else:
        print('\n!! 未找到 101_fix_word_footer.py，请手动修复页脚域代码！')
    return 0


if __name__ == '__main__':
    sys.exit(main())
