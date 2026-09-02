import collections
D=[]
for L in open('directions.txt'):
    if not L.startswith('DIR'): continue
    p=L.split()
    D.append(dict(N=int(p[1]),k=int(p[2]),u=int(p[3]),ul=int(p[4]),uo=int(p[5]),uord=int(p[6]),
                  v=int(p[7]),vl=int(p[8]),vo=int(p[9]),vord=int(p[10]),d=int(p[11]),
                  eps=int(p[12]),cls=p[13],nm=p[14]))
by=collections.defaultdict(list); byeps=collections.defaultdict(list)
for r in D:
    by[(r['N'],r['k'])].append(r); byeps[(r['N'],r['k'],r['eps'])].append(r)

# ---- (A) per-(N,k) census
print("### A. census")
print("| N | dim E3(G1) | #eps(odd) | #INT dirs | distinct interesting wt-3 periods | dim E4(G1) | dim E4(G0) | #INT dirs | distinct interesting wt-4 periods |")
print("|---|---|---|---|---|---|---|---|---|")
for N in range(1,61):
    row=[str(N)]
    for k in (3,4):
        L=by.get((N,k),[])
        ints=sorted({r['nm'] for r in L if r['cls']=='I'})
        neps=len({r['eps'] for r in L})
        ni=sum(1 for r in L if r['cls']=='I')
        if k==3: row += [str(len(L)), str(neps), str(ni), ", ".join(ints) or "-"]
        else:
            g0=sum(1 for r in L if r['eps'] in (1,0))
            row += [str(len(L)), str(g0), str(ni), ", ".join(ints) or "-"]
    print("| "+" | ".join(row)+" |")

# ---- (B) quadratic-only reachable sets
print()
print("### B. quadratic interesting periods reachable")
print("| N | wt3: L(2,chi_-D) available | wt4: L(3,chi_D) available |")
print("|---|---|---|")
for N in range(1,61):
    a=sorted({r['u'] for r in by.get((N,3),[]) if r['cls']=='I' and r['uord']==2})
    b=sorted({r['u'] for r in by.get((N,4),[]) if r['cls']=='I' and r['uord']==2 and r['u']>1})
    print(f"| {N} | "+(", ".join(f"-{x}" for x in a) or "-")+" | zeta(3)"+("".join(f", {x}" for x in b))+" |")

# ---- (C) smallest levels
print()
print("### C. smallest level realising each named value")
want3=[7,8,11,15,20,23,24]; want4=[5,8,12,13,17,24]
for Dd in want3:
    lv=[N for N in range(1,61) if any(r['cls']=='I' and r['u']==Dd and r['uord']==2 and r['uo']==1 for r in by.get((N,3),[]))]
    print(f"L(2,chi_-{Dd}): smallest N = {lv[0] if lv else 'none<=60'} ; all N<=60: {lv}")
for Dd in want4:
    lv=[N for N in range(1,61) if any(r['cls']=='I' and r['u']==Dd and r['uord']==2 and r['uo']==0 for r in by.get((N,4),[]))]
    print(f"L(3,chi_{Dd}) : smallest N = {lv[0] if lv else 'none<=60'} ; all N<=60: {lv}")

# ---- (D) within-nebentypus interesting counts
print()
print("### D. within a fixed nebentypus: distinct interesting periods")
mx=collections.Counter()
viol=[]
for (N,k,e),L in byeps.items():
    s={r['nm'] for r in L if r['cls']=='I'}
    mx[(k,len(s))]+=1
    if len(s)>1: viol.append((N,k,e,sorted(s)))
print("counts of (weight, #distinct interesting periods in one nebentypus):", dict(sorted(mx.items())))
print("cases with >=2 :", viol[:20], " total:",len(viol))

# ---- (E) weight-3 (zeta(2), L(2,chi_-D)) same-nebentypus pairs
print()
print("### E. weight 3: pairs (pi^2/zeta(2), L(2,psi)) inside ONE nebentypus")
bad=[]; good=collections.defaultdict(list)
for (N,k,e),L in byeps.items():
    if k!=3: continue
    has_int={r['nm'] for r in L if r['cls']=='I'}
    has_z2=any(r['cls']=='E' and r['u']==1 for r in L)
    if has_int and has_z2:
        for nm in has_int: good[N].append((e,nm))
    elif has_int and not has_z2: bad.append((N,e,sorted(has_int)))
print("levels N where SOME nebentypus carries both zeta(2) and an interesting L(2,psi):", sorted(good))
print("nebentypus components carrying an interesting L(2,psi) but NOT zeta(2):", bad)
# specifically quadratic
print()
print("quadratic version: N such that eps=chi_-D gives {zeta(2), L(2,chi_-D)} in one nebentypus")
for Dd in [3,4,7,8,11,15,19,20,23,24,31,35,39,40,43,47,51,52,55,56,59]:
    lv=[N for N in range(1,61) if any((r['cls']=='I' and r['u']==Dd and r['uord']==2 and r['uo']==1) for r in by.get((N,3),[]))]
    if lv: print(f"  chi_-{Dd}: N = {lv}")

# ---- (F) weight-4 mixed pairs across nebentypus
print()
print("### F. weight 4: (zeta(3), L(3,chi)) mixed pairs across nebentypus components")
for N in range(1,61):
    L=by.get((N,4),[])
    ints=sorted({r['nm'] for r in L if r['cls']=='I'})
    if len(ints)>=2:
        epsofnm={}
        for r in L:
            if r['cls']=='I': epsofnm.setdefault(r['nm'],set()).add(r['eps'])
        same=any(len(v)>1 for v in epsofnm.values())
        print(f"  N={N}: {len(ints)} interesting: {', '.join(ints)}   [each period sits in exactly one nebentypus: {not same}]")
