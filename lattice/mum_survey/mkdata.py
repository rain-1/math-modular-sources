import json, re, sys
d = json.load(open('cyops.json'))
def norm(s):
    s = s.replace(' ','').replace('^','^')
    # implicit multiplication: digit or ')' followed by '(' or letter
    s = re.sub(r'(\d)\(', r'\1*(', s)
    s = re.sub(r'\)\(', r')*(', s)
    s = re.sub(r'(\d)X', r'\1*X', s)
    s = re.sub(r'\)X', r')*X', s)
    s = re.sub(r'X\(', r'X*(', s)
    return s
rows=[]
for r in d:
    if r.get('aesz') is None: continue
    if r.get('degtheta') != 4: continue
    pl = r.get('pols_list')
    if not pl: continue
    try:
        pols = [norm(p) for p in pl]
    except Exception:
        continue
    rows.append((r['aesz'], r['nn'], r.get('degz'), pols, r.get('h3'), r.get('c2h'), r.get('c3'),
                 r.get('discriminant_list'), r.get('riemann'), (r.get('sol_list') or [])[:8]))
rows.sort(key=lambda x: int(x[0]) if str(x[0]).isdigit() else 10**6)
with open('ops.gp','w') as f:
    f.write('OPS = [\n')
    for a,nn,dz,pols,h3,c2h,c3,disc,riem,sol in rows:
        sing = []
        if riem:
            for e in riem:
                s=e.get('sing')
                if s not in ('0','infinity'): sing.append(s)
        f.write('["%s","%s",%s,[%s],"%s","%s","%s",[%s],[%s]],\n' % (
            a, nn, dz if dz is not None else -1,
            ','.join(pols),
            h3 if h3 is not None else '', c2h if c2h is not None else '', c3 if c3 is not None else '',
            ','.join('"%s"'%x for x in sing),
            ','.join(str(x) for x in sol)))
    f.write('0];\nOPS = OPS[1..#OPS-1];\n')
print(len(rows),'operators of order 4 with AESZ numbers')
