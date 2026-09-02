from mpmath import mp, mpf, mpc, sqrt, log, cos, sin, pi, quad, inf
import pickle
mp.dps = 90
NDIG = 60

def sortab(alist):
    A = sorted([mpf(a) for a in alist])          # ascending a  -> descending root
    R = sorted([1/mpf(a) for a in alist])        # ascending root
    Aroot = sorted([mpf(a) for a in alist], reverse=True)  # a matching ascending root
    return Aroot, R

def I0r1(alist, kern):
    A,R = sortab(alist); a1=A[0]; r1=R[0]; rest=A[1:]
    S=sqrt(r1)
    def g(s):
        t=r1-s*s; Q=mpf(1)
        for z in rest: Q*=(1-z*t)
        return 2*kern(t)/(sqrt(a1)*sqrt(Q))
    return quad(g,[0,S],maxdegree=13)

def Ibet(alist,i,j,kern):
    A,R = sortab(alist); ri,rj=R[i],R[j]; ai,aj=A[i],A[j]
    oth=[A[k] for k in range(len(A)) if k not in (i,j)]
    def g(phi):
        u=(1-cos(phi))/2; t=ri+(rj-ri)*u; Q=mpf(1)
        for z in oth: Q*=abs(1-z*t)
        return kern(t)/(sqrt(ai*aj)*sqrt(Q))
    return quad(g,[0,pi],maxdegree=13)

KB=lambda t:1/(1-t); KD=lambda t:log(1-t)/(1-t); KL=lambda t:log(1-t); K1=lambda t:mpf(1)

HOSTS={}
for m in [1,2,3,4,5,6,12]: HOSTS['E_%d'%m]=[1,9,4*m]
for m in [1,2,3]: HOSTS['F_%d'%m]=[1,25,4*m]
HOSTS['G']=[4,8,12]

def FP_branchpoint(alist):
    """E_m,F_m: 1 is a root. finite part of int_0^1 H/(1-t) dt, principal continuation from UHP."""
    A,R=sortab(alist)     # A[0]>A[1]>A[2]=1 ; R[2]=1
    ahi,amid=A[0],A[1]; r1,r2=R[0],R[1]
    g1 = -1/sqrt((ahi-1)*(amid-1))     # g(1)
    u = sqrt(1-r2)
    def integ(s):
        t=r2+s*s
        twosg = -2/sqrt(amid*(ahi*t-1))     # = 2s*g(t)
        return (twosg - 2*s*g1)*(1-t)**mpf('-1.5')
    tail = quad(integ,[0,u],maxdegree=14)
    FPlast = tail - 2*g1/sqrt(1-r2)
    p0 = I0r1(alist,KB)
    p1 = Ibet(alist,0,1,KB)
    return mpc(p0 + FPlast, p1), g1, p0, p1, FPlast

def FP_log(alist):
    """G: 1 not a root. FP = lim [ int_0^x H/(1-t) + H(1) log(1-x) ]."""
    A,R=sortab(alist); r1,r2,r3=R
    h1 = 1/sqrt((A[0]-1)*(A[1]-1)*(A[2]-1))   # |P(1)|^{-1/2}
    def h(t): return 1/sqrt(abs((1-A[0]*t)*(1-A[1]*t)*(1-A[2]*t)))
    u=sqrt(1-r3)
    def integ(s):
        t=r3+s*s
        # 2s*h(t) = 2/sqrt(A? ) : the vanishing factor is (1-A[2]*t)? root r3 <-> a = A[2]
        a3=A[2]
        two_s_h = 2/sqrt(abs((1-A[0]*t)*(1-A[1]*t))*a3)
        return (two_s_h - 2*s*h1)/(1-t)
    tail=quad(integ,[0,u],maxdegree=14)
    FPlast = -mpc(0,1)*(tail + h1*log(1-r3))
    p0=I0r1(alist,KB); p1=Ibet(alist,0,1,KB); p2=Ibet(alist,1,2,KB)
    return mpc(p0,0)+mpc(0,1)*p1 - p2 + FPlast, h1, p0,p1,p2,FPlast

res={}
for name,al in HOSTS.items():
    A,R=sortab(al)
    d={'alist':al,'roots':R}
    for tag,k in [('B',KB),('D',KD),('L',KL),('1',K1)]:
        d['A'+tag]=I0r1(al,k); d['B'+tag]=Ibet(al,0,1,k)
        if name=='G': d['C'+tag]=Ibet(al,1,2,k)
    # egg period for E/F : between r2 and r3 with k=1 (r3=1)
    if name!='G':
        d['egg']=Ibet(al,1,2,K1)
    else:
        d['egg']=None
    if name=='G':
        fp,h1,p0,p1,p2,fpl=FP_log(al); d['FP']=fp; d['FPtype']='log'; d['H1']=h1
    else:
        fp,g1,p0,p1,fpl=FP_branchpoint(al); d['FP']=fp; d['FPtype']='sqrt'; d['g1']=g1
    res[name]=d

with open('res.pkl','wb') as f: pickle.dump({k:{kk:(str(vv) if not isinstance(vv,(list,type(None),str)) else vv) for kk,vv in v.items()} for k,v in res.items()},f)

def S(x,n=NDIG): return mp.nstr(x,n)
for name,d in res.items():
    print('='*80)
    print(name,' a=',d['alist'],' roots=',[S(r,20) for r in d['roots']])
    print('  Omega_fold  = int_{r1}^{r2} dx/sqrt|P|  =', S(d['B1']))
    if d['egg'] is not None:
        print('  Omega_egg   = int_{r2}^{r3} dx/sqrt P   =', S(d['egg']))
    else:
        print('  Omega_23    = int_{r2}^{r3} dx/sqrt P   =', S(d['C1']))
        print('  Omega_01    = int_{0}^{r1}  dx/sqrt P   =', S(d['A1']))
    print('  FP int_0^1 H/(1-t) dt  [%s-type] ='%d['FPtype'])
    print('      Re =', S(d['FP'].real))
    print('      Im =', S(d['FP'].imag))
