import os, time, urllib.request, urllib.error

BASE = 'docking_files'
os.makedirs(f'{BASE}/proteins', exist_ok=True)
os.makedirs(f'{BASE}/ligands', exist_ok=True)

UA = {'User-Agent': 'Mozilla/5.0'}

def fetch(url, tries=4, timeout=30):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as f:
                return f.read()
        except Exception as e:
            print(f'   retry {i+1}: {str(e)[:60]}')
            time.sleep(1.5 + 2 * i)
    return None

# ---------------- 蛋白（受体）----------------
PROTEINS = [
    ('PPARG', '2PRG', 'PPARG LBD + rosiglitazone (primary docking target)'),
    ('PPARG', '3ET1', 'PPARG LBD + GW1929 (alternative)'),
    ('HIF1A', '4H6J', 'HIF-1alpha PAS-B domain'),
    ('HIF1A', '5W87', 'HIF-1alpha PAS-B (alternative)'),
    ('RUNX1', '1CCP', 'RUNX1 runt domain + DNA'),
    ('RUNX1', '3L8H', 'RUNX1 runt domain (alternative)'),
]

def clean_pdb(raw: bytes) -> bytes:
    """仅保留 ATOM（蛋白本身），去掉 HETATM（原配体/水/离子）"""
    out = []
    for line in raw.decode('utf-8', 'ignore').splitlines():
        if line.startswith('ATOM') or line.startswith('TER') or line.startswith('END'):
            out.append(line)
    return ('\n'.join(out) + '\n').encode('utf-8')

print('=== 下载蛋白 PDB ===')
for gene, pid, note in PROTEINS:
    raw = fetch(f'https://files.rcsb.org/download/{pid}.pdb')
    if not raw or b'ATOM' not in raw:
        print(f'   [FAIL] {gene}_{pid}')
        continue
    p1 = f'{BASE}/proteins/{gene}_{pid}.pdb'
    p2 = f'{BASE}/proteins/{gene}_{pid}_clean.pdb'
    open(p1, 'wb').write(raw)
    open(p2, 'wb').write(clean_pdb(raw))
    n_atom = sum(1 for l in raw.decode('utf-8', 'ignore').splitlines() if l.startswith('ATOM'))
    print(f'   [OK] {gene}_{pid}: raw {len(raw)//1024}KB, ATOM={n_atom} -> +_clean.pdb(仅受体)')
    time.sleep(0.4)

# ---------------- 配体（3D mol from ChEMBL）----------------
LIGANDS = [
    ('Kaempferol', 'CHEMBL150', 'PPARG 主配体 (14.1 nM) / HIF1A 备选'),
    ('Quercetin', 'CHEMBL50', 'HIF1A 配体 (10 uM)'),
    ('Protocatechuic_acid', 'CHEMBL37537', 'RUNX1 配体 (4 uM)'),
]

print('\n=== 下载配体 3D mol ===')
for name, cid, note in LIGANDS:
    data = fetch(f'https://www.ebi.ac.uk/chembl/api/data/molecule/{cid}.mol3D', tries=3)
    fmt = '3D'
    if not data or b'V2000' not in data and b'V3000' not in data:
        data = fetch(f'https://www.ebi.ac.uk/chembl/api/data/molecule/{cid}.sdf', tries=3)
        fmt = '2D-SDF(网站需自行生成3D)'
    if not data:
        print(f'   [FAIL] {name}_{cid}')
        continue
    ext = 'mol' if fmt == '3D' else 'sdf'
    p = f'{BASE}/ligands/{name}_{cid}_3D.{ext}'
    open(p, 'wb').write(data)
    print(f'   [OK] {name}_{cid}: {len(data)}B ({fmt}) | {note}')
    time.sleep(0.8)

print('\n=== 清单 ===')
for root, _, files in os.walk(BASE):
    for f in sorted(files):
        fp = os.path.join(root, f)
        print(f'   {fp}  ({os.path.getsize(fp)//1024}KB)')
