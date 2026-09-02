from ident import *
setup(140)
HOSTS = [(3,9,'M'),(3,63,'P'),(3,126,'M'),(4,80,'P'),(4,624,'P')]
G = mp.catalan
CL = im(polylog(2, exp(mpc(0,1)*pi/3)))   # Cl_2(pi/3)
print("Cl2(pi/3) =",mp.nstr(CL,25),"   G =",mp.nstr(G,25))
for (k,N,fam) in HOSTS:
    h=Host(k,N,fam); n=int(mp.nint(mpf(h.D)**(mpf(1)/k)))
    REF = CL if k==3 else G; RNAME = "Cl2(pi/3)" if k==3 else "G"
    print(f"\n########## k={k} N={N} {fam} D={h.D}={n}^{k} ##########")
    print("  prime set:",primeset(h))
    # 1. classify each distinct Bloch point
    P=parts(h,1); pts=P['pts']
    seen=[]
    for z in pts:
        d=BW(z)
        if abs(d)<mpf(10)**-60: cls="D(z)=0"
        else:
            r=pslq([d,REF],tol=mpf(10)**-100,maxcoeff=10**6)
            cls = f"D(z) = {Rational(-int(r[1]),int(r[0]))} * {RNAME}" if r else "*** NOT a rational multiple of "+RNAME
        print(f"    z={mp.nstr(z,14):40s} D(z)={mp.nstr(d,20):26s}  {cls}")
