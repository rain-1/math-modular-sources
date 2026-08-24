#!/usr/bin/env python3
"""12_window.py -- pull the |lambda_2| < 1 rows out of the window scan and label
them: genuine (five singular points) / disguised, real lambda_2 or not, score."""
import sys, glob, json, importlib.util, os
from fractions import Fraction as Fr
import mpmath as mp
HERE=os.path.dirname(os.path.abspath(__file__))
spec=importlib.util.spec_from_file_location("an",os.path.join(HERE,"06_analyse_deep.py"))
an=importlib.util.module_from_spec(spec); spec.loader.exec_module(an)
mp.mp.dps=60
seen=set(); out=[]
for p in sys.argv[1:]:
    for line in open(p):
        if not line.strip() or line.startswith('#'): continue
        key=an.canon(tuple(int(x) for x in line.split()))
        if key in seen: continue
        seen.add(key)
        pr=an.parse(" ".join(map(str,key)))
        if pr is None: continue
        kind,cls,row,co=pr
        a_,d_,g_=co['a'],co['d'],co['g']
        lams=mp.polyroots([1,-a_,d_,-g_],maxsteps=300,extraprec=300)
        lams=sorted(lams,key=lambda z:-abs(z))
        if abs(lams[1])>=1: continue
        try: r=an.analyse(co,full=True)
        except Exception as ex:
            out.append(dict(cls=cls,row=row,err=str(ex))); continue
        if r is None: continue
        r['kind']=kind; r['cls']=cls; r['row']=row
        if r['disc']==0: r['verdict']='REPEATED'
        elif any(r['apparent']) or r['inf_apparent']: r['verdict']='DISGUISED'
        else: r['verdict']='GENUINE'
        out.append(r)
print("rows with |lambda_2| < 1 :", len(out))
for r in sorted(out,key=lambda z:-(z.get('score') or -99)):
    if 'err' in r: print("  ERR",r['cls'],r['row'],r['err']); continue
    print("  %-22s %-26s |l2|=%.5f real_l2=%-5s k=%s score=%+.4f  %s  disc=%s  xi=%s"%(
        str(r['cls']),str(r['row']),abs(complex(*r['lam'][1])),r['real_l2'],r.get('k'),
        r.get('score') or -99,r['verdict'],r['disc'],(r.get('xi') or 'none')[:34]))
json.dump(out,open(os.path.join(HERE,'out','window_rows.json'),'w'))
