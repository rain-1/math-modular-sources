from mpmath import mp, mpf, mpc, sqrt, log, cos, sin, pi, quad, pslq
mp.dps=90
exec(open('master.py').read().split('KB=lambda')[0])
K1=lambda t:mpf(1); Kx=lambda t:t
KB=lambda t:1/(1-t); KD=lambda t:log(1-t)/(1-t); KL=lambda t:log(1-t)
HOSTS={}
for m in [1,2,3,4,5,6,12]: HOSTS['E_%d'%m]=[1,9,4*m]
for m in [1,2,3]: HOSTS['F_%d'%m]=[1,25,4*m]
HOSTS['G']=[4,8,12]
import pickle
data={}
for name,al in HOSTS.items():
    d={}
    for tag,k in [('1',K1),('x',Kx),('B',KB),('D',KD),('L',KL)]:
        d['A'+tag]=I0r1(al,k)
        d['B'+tag]=Ibet(al,0,1,k)
        d['C'+tag]=Ibet(al,1,2,k)
    data[name]=d
    print('---',name)
    O,E = d['B1'],d['Bx']
    r=pslq([d['BB'],O,E], maxcoeff=10**8, maxsteps=10**6, tol=mpf(10)**-70)
    print('  PSLQ  B_B vs (Omega,eta)   ->', r)
    r2=pslq([d['AB'],d['A1'],d['Ax'],mpf(1)], maxcoeff=10**8, maxsteps=10**6, tol=mpf(10)**-70)
    print('  PSLQ  A_B vs (A_1,A_x,1)   ->', r2)
    r3=pslq([d['CB'],d['C1'],d['Cx']], maxcoeff=10**8, maxsteps=10**6, tol=mpf(10)**-70)
    print('  PSLQ  C_B vs (C_1,C_x)     ->', r3)
pickle.dump({k:{kk:str(vv) for kk,vv in v.items()} for k,v in data.items()},open('data.pkl','wb'))
