from mpmath import mp, mpf, mpc, sqrt, log, pi, pslq, nstr
import pickle
mp.dps=110
raw=pickle.load(open('D.pkl','rb'))
D={}
for k,v in raw.items():
    D[k]={kk:(mpf(vv) if isinstance(vv,str) else vv) for kk,vv in v.items()}

def trial(label, names, vals, maxcoeff=10**6, tol=None, ndig=70):
    if tol is None: tol=mpf(10)**(-ndig)
    r=pslq(vals, maxcoeff=maxcoeff, maxsteps=10**6, tol=tol)
    if r is None:
        print('   %-46s -> none'%label); return None
    res=sum(mpf(c)*v for c,v in zip(r,vals))
    ok = abs(res) < mpf(10)**(-ndig+5)
    s=' + '.join('%d*%s'%(c,n) for c,n in zip(r,names) if c!=0)
    print('   %-46s -> %s  = 0   (resid %s) %s'%(label,s,nstr(abs(res),3),'OK' if ok else '??'))
    return r

L2=log(mpf(2)); L3=log(mpf(3)); L5=log(mpf(5)); PI=pi
for name in D:
    d=D[name]; m=d['m']
    O=d['B1']; E=d['Bx']; Lp=d['Lp0']
    print('='*70); print(name, d['alist'], ' rank',d['rank'],' N',d['N'])
    print('  B_D =',nstr(d['BD'],30),'  B_L =',nstr(d['BL'],30))
    # 1) B_L in Omega,eta with log coefficients
    trial('B_L ~ (O,e)(1,log2,log3)',['B_L','O','e','Ol2','el2','Ol3','el3'],
          [d['BL'],O,E,O*L2,E*L2,O*L3,E*L3])
    # 2) B_L with L'(E,0)
    trial('B_L ~ O,e,Lp',['B_L','O','e','Lp'],[d['BL'],O,E,Lp])
    trial('B_L ~ O,e,Lp,pi',['B_L','O','e','Lp','pi'],[d['BL'],O,E,Lp,PI])
    # 3) B_D
    trial('B_D ~ (O,e)(1,log2,log3)',['B_D','O','e','Ol2','el2','Ol3','el3'],
          [d['BD'],O,E,O*L2,E*L2,O*L3,E*L3])
    trial('B_D ~ O,e,Lp',['B_D','O','e','Lp'],[d['BD'],O,E,Lp])
    trial('B_D ~ B_L,O,e',['B_D','B_L','O','e'],[d['BD'],d['BL'],O,E])
    trial('B_D,B_L ~ O,e,pi',['B_D','B_L','O','e','pi'],[d['BD'],d['BL'],O,E,PI])
