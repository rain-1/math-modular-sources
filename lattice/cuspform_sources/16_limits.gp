\\ The Apery limit of every weight-4 cusp form source, by level (host-independent, Lemma A):
\\ xi(Phi) = Theta(q_c) + (2 Pi/sqrt N) DTheta(q_c),  q_c = exp(-2 Pi/sqrt N).
default(parisizemax, 8000000000);
default(realprecision, 60);
read("lib.gp");
NQ = 400;
lev = [5,6,7,8,9,10,12,18];
for(i=1,#lev, my(N=lev[i], mf, d, Bb, ai, M, Km, Kp, g, cv, Th, qc, xi, K, nm, sg); mf=mfinit([N,4],1); d=mfdim(mf); if(d==0, next); Bb=mfbasis(mf); ai=mfatkininit(mf,N); M=matrix(d,d); for(j=1,d, my(vv=mftobasis(mf,mfatkin(ai,Bb[j]))); for(t=1,d, M[t,j]=vv[t])); Km=matker(M+1); Kp=matker(M-1); qc=exp(-2*Pi/sqrt(N*1.0)); print(""); print("### level ", N, "   dim S_4 = ", d, "   dim(W=-1) = ", matsize(Km)[2]); for(sg=1,2, K=if(sg==1,Km,Kp); nm=if(sg==1,"W=-1 (fold-regular)","W=+1 (not fold-regular)"); for(j=1,matsize(K)[2], g=mflinear(mf, primvec(K[,j]~)~); cv=mfcoefs(g,NQ-1); if(cv[2]<0, cv=-cv; g=mflinear(mf,-primvec(K[,j]~)~)); Th=thetaser(cv,3,NQ); xi = subst(truncate(Th),q,qc) + (2*Pi/sqrt(N*1.0))*subst(truncate(q*deriv(Th,q)),q,qc); print("  ", nm, "  Phi = ", vector(10,t,cv[t])); print("     xi       = ", xi); print("     L(Phi,3) = ", lfun(lfunmf(mf,g),3)))));
quit;
