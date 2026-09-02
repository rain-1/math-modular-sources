from mpmath import mp, mpf, mpc, sqrt, log, cos, sin, pi, quad, pslq
mp.dps=110
exec(open('master.py').read().split('KB=lambda')[0])
K1=lambda t:mpf(1); Kx=lambda t:t; Kxx=lambda t:t*t
KB=lambda t:1/(1-t); KD=lambda t:log(1-t)/(1-t); KL=lambda t:log(1-t)
HOSTS={}
for m in [1,2,3,4,5,6,12]: HOSTS['E_%d'%m]=([1,9,4*m],m,9)
for m in [1,2,3]: HOSTS['F_%d'%m]=([1,25,4*m],m,25)
HOSTS['G']=([4,8,12],None,None)
LV={}
for ln in open('lval.out'):
    f=ln.strip().split('|')
    if len(f)<9: continue
    LV[f[0]]=dict(N=int(f[1]),Lp0=mpf(f[2]),L2=mpf(f[3]),om1min=mpf(f[4]),om2min=mpf(f[5]),
                  om1=mpf(f[6]),om2=mpf(f[7]),rank=int(f[8]))
def build():
    D={}
    for name,(al,m,a) in HOSTS.items():
        d={'alist':al,'m':m}
        for tag,k in [('1',K1),('x',Kx),('xx',Kxx),('B',KB),('D',KD),('L',KL)]:
            d['A'+tag]=I0r1(al,k); d['B'+tag]=Ibet(al,0,1,k)
            if name=='G': d['C'+tag]=Ibet(al,1,2,k)
        d.update(LV[name]); D[name]=d
    return D
