import urllib.request, urllib.parse, time, os

BASE = 'docking_files'
UA = {'User-Agent': 'Mozilla/5.0'}

LIGANDS = [
    ('Kaempferol', 'CHEMBL150', 'O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12'),
    ('Quercetin', 'CHEMBL50', 'O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12'),
    ('Protocatechuic_acid', 'CHEMBL37537', 'O=C(O)c1ccc(O)c(O)c1'),
]

def fetch(url, tries=3, timeout=25):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as f:
                return f.read()
        except Exception as e:
            print('   err:', str(e)[:70])
            time.sleep(1.5)
    return None

print('=== 下载 3D SDF 到 ligands/ ===')
for name, cid, smi in LIGANDS:
    u = ('https://cactus.nci.nih.gov/chemical/structure/'
         + urllib.parse.quote(smi, safe='')
         + '/file?format=sdf&get3d=true')
    data = fetch(u)
    if not data:
        print(f'   [FAIL] {name}')
        continue
    # 简单校验 z 坐标非全零 -> 真 3D
    lines = data.decode('utf-8', 'ignore').splitlines()
    zvals = [float(l[20:30]) for l in lines[4:4 + 60] if len(l) > 30 and l[30:34].strip()]
    is3d = any(abs(z) > 1e-6 for z in zvals)
    p = f'{BASE}/ligands/{name}_{cid}_3D.sdf'
    open(p, 'wb').write(data)
    print(f'   [OK] {name}: {len(data)}B, 真3D={is3d} -> {p}')
    time.sleep(0.6)

# 清掉早先的 2D 版（同名 _3D.sdf 已覆盖；删除旧 2D 残留）
old = f'{BASE}/ligands/Protocatechuic_acid_CHEMBL37537_3D.sdf'
print('\n完成。当前 ligands/:')
for f in sorted(os.listdir(f'{BASE}/ligands')):
    print('  ', f, os.path.getsize(f'{BASE}/ligands/{f}'), 'B')
