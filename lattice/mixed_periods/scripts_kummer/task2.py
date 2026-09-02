from kummer import *
setup(90)
out=open('periods60.txt','w')
def w(s): out.write(s+"\n"); 
w("Fold periods c^{(a)}[kappa] to 60 significant digits.")
w("Family M: H=(1-Nx)^{-1/k}, delta=+1/N, D=N-1.   Family P: H=(1+Nx)^{-1/k}, delta=-1/N, D=N+1.")
w("kernels  B=1/(1-t)   D=log(1-t)/(1-t)   L=log(1-t)")
w("")
maxdiff=0
for k in [3,4,5,6]:
    for m in [1,2,3]:
        N = k**k*m
        for fam in ['M','P']:
            h=Host(k,N,fam)
            w(f"--- k={k} m={m} N=k^k*m={N} fam={fam} D={h.D} ---")
            for a in range(1,k):
                vals={}
                for kern in ['B','D','L']:
                    v = {'B':h.cB,'D':h.cD,'L':h.cL}[kern](a)
                    assert abs(im(v))<mpf(10)**-70, (k,m,fam,a,kern,im(v))
                    v=re(v); vals[kern]=v
                    q = h.per_u(a,kern)
                    maxdiff=max(maxdiff,float(abs(v-q)/abs(v)))
                w(f"  a={a}  B = {mp.nstr(vals['B'],60)}")
                w(f"       D = {mp.nstr(vals['D'],60)}")
                w(f"       L = {mp.nstr(vals['L'],60)}")
w("")
w(f"max relative difference (closed form vs u-quadrature) over the whole table: {maxdiff:.3e}")
out.close()
print(open('periods60.txt').read()[:2000])
print("...")
print("maxdiff", maxdiff)
