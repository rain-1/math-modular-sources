from kummer import *
setup(60)
tests = [(2,4,'M'),(2,4,'P'),(3,27,'M'),(3,27,'P'),(3,9,'M'),(3,63,'P'),(4,256,'M'),(4,80,'P'),(5,3125,'M'),(5,1025,'P')]
for (k,N,fam) in tests:
    h = Host(k,N,fam)
    print(f"k={k} N={N} fam={fam} D={h.D} w={mp.nstr(h.w(),5)}")
    for a in range(1,k):
        for kern in ['B','D','L']:
            u = h.per_u(a,kern); t = h.per_t(a,kern)
            print(f"   a={a} {kern}: u-quad={mp.nstr(u,25)}  t-quad={mp.nstr(t,25)}  diff={mp.nstr(abs(u-t),3)}")
