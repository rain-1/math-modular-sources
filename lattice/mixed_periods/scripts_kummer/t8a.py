from decomp import *
setup(110)
HOSTS = [(3,9,'M'),(3,63,'P'),(3,126,'M'),(4,80,'P'),(4,624,'P'),(5,1025,'M'),(5,7775,'P'),(6,15624,'P')]
for (k,N,fam) in HOSTS:
    h=Host(k,N,fam); n=int(mp.nint(mpf(h.D)**(mpf(1)/k)))
    print(f"\n=== k={k} N={N} {fam} D={h.D}={n}^{k} ===")
    for a in range(1,k):
        P=parts(h,a); cd=h.cD(a)
        chk = abs(P['total']-cd)
        rec = P['LiRe']+P['Bloch']+P['ArgP']+P['ReE']
        print(f" a={a}: decomp check {mp.nstr(chk,3)}   reassembly err {mp.nstr(abs(rec-re(cd)),3)}")
        print(f"      c_D          = {mp.nstr(re(cd),30)}")
        print(f"      Bloch part   = {mp.nstr(P['Bloch'],30)}")
        print(f"      Re Li2 part  = {mp.nstr(P['LiRe'],30)}")
        print(f"      arg*log part = {mp.nstr(P['ArgP'],30)}")
        print(f"      Re E         = {mp.nstr(P['ReE'],30)}")
        # coefficients of Li2: are Im(lam) rational multiples of 1/n^{k-a}?
        ims = sorted(set(mp.nstr(im(l)*n**(k-a),12) for l in P['lams']))
        print(f"      Im(lam)*n^(k-a) values: {ims}")
