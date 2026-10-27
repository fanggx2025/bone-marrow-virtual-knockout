#!/usr/bin/env python3
import urllib.request, urllib.parse, ssl, time

HOST="https://swissdock.ch:8443"
CTX=ssl.create_default_context(); CTX.check_hostname=False; CTX.verify_mode=ssl.CERT_NONE

def call(path, params=None, fields=None):
    url=HOST+path
    if params: url+="?"+urllib.parse.urlencode(params)
    data=None; headers={}
    if fields:
        boundary="----B"+str(int(time.time()*1000)); parts=[]
        for k,v in fields.items():
            parts.append(("--"+boundary).encode()); parts.append(('Content-Disposition: form-data; name="%s"\r\n\r\n%s'%(k,v)).encode())
        parts.append(("--"+boundary+"--\r\n").encode()); data=b"\r\n".join(parts)
        headers={"Content-Type":"multipart/form-data; boundary="+boundary}
    req=urllib.request.Request(url,data=data,headers=headers)
    try:
        with urllib.request.urlopen(req,context=CTX,timeout=90) as r: return r.status,r.read().decode("utf-8","replace")
    except urllib.error.HTTPError as e: return e.code,e.read().decode("utf-8","replace")
    except Exception as e: return -1,str(e)

MINPDB="""ATOM      1  N   ALA A   1      11.104  13.207  10.567  1.00  0.00           N
ATOM      2  CA  ALA A   1      11.622  14.581  10.588  1.00  0.00           C
ATOM      3  C   ALA A   1      13.127  14.692  10.622  1.00  0.00           C
ATOM      4  O   ALA A   1      13.627  15.843  10.859  1.00  0.00           O
ATOM      5  CB  ALA A   1      11.011  15.472   9.490  1.00  0.00           C
END
"""

full=open(r"docking_files/wuweizi_vina/receptor/PPARG_2PRG_chainA_apo.pdb").read()

for label,pdb in [("MINIMAL",MINPDB),("FULL_CHAIN_A",full)]:
    st,sid=call("/opensession"); SID=sid.strip().split("\n")[0].strip()
    st,r=call("/preptarget",params={"sessionNumber":SID},fields={"pdbContent":pdb})
    print(f"[{label}] preptarget ->",st,repr(r[:120]))
    # try setparameters to see if target registered
    st2,r2=call("/setparameters",params={"sessionNumber":SID,"exhaust":"1","method":"AC"})
    print(f"    setparameters ->",repr(r2[:120]))
