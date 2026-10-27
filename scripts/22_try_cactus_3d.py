import urllib.request, time

UA = {'User-Agent': 'Mozilla/5.0'}

LIGANDS = [
    ('Kaempferol', 'CHEMBL150', 'O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12'),
    ('Quercetin', 'CHEMBL50', 'O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12'),
    ('Protocatechuic_acid', 'CHEMBL37537', 'O=C(O)c1ccc(O)c(O)c1'),
]

def fetch(url, tries=2, timeout=25):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as f:
                return f.read()
        except Exception as e:
            print('   err:', str(e)[:70])
            time.sleep(1.5)
    return None

for name, cid, smi in LIGANDS:
    smi_q = urllib.parse.quote(smi, safe='')
    u = f'https://cactus.nci.nih.gov/chemical/structure/{smi_q}/file?format=sdf&get3d=true'
    data = fetch(u)
    ok = data and b'V2000' in data or (data and b'V3000' in data)
    print(name, 'OK' if ok else 'FAIL', len(data) if data else 0)
    if ok:
        # 3D SDF 里若坐标全为 0 则其实是 2D
        lines = data.decode('utf-8', 'ignore').splitlines()
        print('   first atom line:', lines[4][:60] if len(lines) > 4 else '?')
    time.sleep(0.5)
