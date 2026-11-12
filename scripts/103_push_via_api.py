"""通过 GitHub API 提交改动（git clone 走代理隧道不稳定时的备用通道）。

用法：
  python scripts/103_push_via_api.py scripts/audit_merged_pdf.py

只处理单个文件；内容从本地 github_upload/ 下的同名路径读取。
"""
import base64
import hashlib
import json
import os
import subprocess
import sys
import urllib.error
import urllib.request

REPO = 'fanggx2025/bone-marrow-virtual-knockout'
API = 'https://api.github.com/repos/%s/git' % REPO
BASE = r'C:\Users\Fang\Desktop\S3\ciwujia_wuweizi_oa'
SRC = os.path.join(BASE, 'github_upload')


def get_token():
    tok = os.environ.get('GITHUB_TOKEN')
    if tok:
        return tok
    out = subprocess.run(['git', 'credential', 'fill'],
                         input='protocol=https\nhost=github.com\n\n',
                         capture_output=True, text=True).stdout
    for line in out.splitlines():
        if line.startswith('password='):
            return line.split('=', 1)[1]
    return None


def api(method, url, token, payload=None):
    req = urllib.request.Request(url, method=method)
    req.add_header('Authorization', 'Bearer ' + token)
    req.add_header('Accept', 'application/vnd.github+json')
    req.add_header('X-GitHub-Api-Version', '2022-11-28')
    data = None
    if payload is not None:
        data = json.dumps(payload).encode()
        req.add_header('Content-Type', 'application/json')
    try:
        with urllib.request.urlopen(req, data, timeout=60) as r:
            return r.status, json.loads(r.read().decode())
    except urllib.error.HTTPError as e:
        return e.code, e.read().decode()[:500]


def main():
    if len(sys.argv) < 2:
        print('用法：python 103_push_via_api.py <相对路径> [<相对路径> ...]')
        return 2
    token = get_token()
    if not token:
        print('拿不到凭据')
        return 1

    # 1. 取当前 ref
    st, ref = api('GET', '%s/ref/heads/main' % API, token)
    if st != 200:
        print('取 ref 失败', st, ref)
        return 1
    sha = ref['object']['sha']

    # 2. 取当前 commit
    st, cm = api('GET', '%s/commits/%s' % (API, sha), token)
    if st != 200:
        print('取 commit 失败', st, cm)
        return 1
    base_tree = cm['tree']['sha']

    # 3. 构造新 tree（只加/改指定文件）
    tree_entries = []
    for rel in sys.argv[1:]:
        local = os.path.join(SRC, rel.replace('/', os.sep))
        if not os.path.isfile(local):
            print('找不到本地文件：%s' % local)
            return 1
        blob = open(local, 'rb').read()
        st, b = api('POST', '%s/blobs' % API, token,
                    {'content': base64.b64encode(blob).decode(),
                     'encoding': 'base64'})
        if st not in (200, 201):
            print('建 blob 失败', rel, st, b)
            return 1
        tree_entries.append({'path': rel, 'mode': '100644',
                             'type': 'blob', 'sha': b['sha']})
        print('  blob OK: %s' % rel)

    st, nt = api('POST', '%s/trees' % API, token,
                 {'base_tree': base_tree, 'tree': tree_entries})
    if st not in (200, 201):
        print('建 tree 失败', st, nt)
        return 1

    # 4. 建 commit
    msg = ('Fix submission-audit checks: evaluate the retired Zenodo DOI per page '
           'and treat the Related Manuscript copy as an optional attachment.\n\n'
           'Editorial Manager does not merge attachment bodies into the combined '
           'submission PDF, so the previous version of this check reported a false '
           'failure for the related-manuscript copy.')
    st, c = api('POST', '%s/commits' % API, token,
                {'message': msg, 'tree': nt['sha'], 'parents': [sha]})
    if st not in (200, 201):
        print('建 commit 失败', st, c)
        return 1

    # 5. 更新 ref
    st, r = api('PATCH', '%s/refs/heads/main' % API, token,
                {'sha': c['sha'], 'force': False})
    if st not in (200, 201):
        print('更新 ref 失败', st, r)
        return 1
    print('\n已推送 commit: %s' % c['sha'][:12])

    # 6. 验证
    st, tree = api('GET', '%s/trees/main?recursive=1' % API, token)
    if st == 200:
        files = {t['path'] for t in tree.get('tree', []) if t['type'] == 'blob'}
        print('远端文件总数: %d' % len(files))
        for rel in sys.argv[1:]:
            print('  %-42s %s' % (rel, 'OK' if rel in files else 'MISSING'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
