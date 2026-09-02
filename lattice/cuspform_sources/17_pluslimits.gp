\\ Are the W_N = +1 (non fold-regular) Apery limits identifiable?  60 digits available.
default(parisizemax, 8000000000);
default(realprecision, 60);
read("lib.gp");
NQ = 400;
lev = [5,6,7,8,9,10,12,18];
for(i=1,#lev, my(N=lev[i], mf, d, Bb, ai, M, Kp, g, cv, Th, qc, xi, L1,L2,L3); mf=mfinit([N,4],1); d=mfdim(mf); if(d==0, next); Bb=mfbasis(mf); ai=mfatkininit(mf,N); M=matrix(d,d); for(j=1,d, my(vv=mftobasis(mf,mfatkin(ai,Bb[j]))); for(t=1,d, M[t,j]=vv[t])); Kp=matker(M-1); qc=exp(-2*Pi/sqrt(N*1.0)); for(j=1,matsize(Kp)[2], g=mflinear(mf, primvec(Kp[,j]~)~); cv=mfcoefs(g,NQ-1); if(cv[2]<0, cv=-cv; g=mflinear(mf,-primvec(Kp[,j]~)~)); Th=thetaser(cv,3,NQ); xi = subst(truncate(Th),q,qc) + (2*Pi/sqrt(N*1.0))*subst(truncate(q*deriv(Th,q)),q,qc); L1=lfun(lfunmf(mf,g),1); L2=lfun(lfunmf(mf,g),2); L3=lfun(lfunmf(mf,g),3); print("N=",N," Phi=",vector(6,t,cv[t]),"  xi=",xi); print("     lindep[xi,L1,L2,L3]        = ", lindep([xi,L1,L2,L3],40)); print("     lindep[xi,L3,L2/Pi^2]      = ", lindep([xi,L3,L2/Pi^2],45))));
quit;
