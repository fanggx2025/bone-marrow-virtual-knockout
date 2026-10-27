import urllib.request, urllib.error, time

ids = ['2PRG', '3ET1', '1FM9', '4EMA', '4H6J', '5W87', '1CCP', '3L8H']
for i in ids:
    u = 'https://files.rcsb.org/download/%s.pdb' % i
    try:
        req = urllib.request.Request(u, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=15) as f:
            head = f.read(200).decode('utf-8', 'ignore')
        ok = ('HEADER' in head) or ('ATOM' in head)
        print(i, 'OK' if ok else 'RETURNED-BUT-NO-HEADER')
    except urllib.error.HTTPError as e:
        print(i, 'HTTP', e.code)
    except Exception as e:
        print(i, 'ERR', str(e)[:50])
    time.sleep(0.3)
