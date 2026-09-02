\\ (a) exact inhomogeneous recurrence sum_j P_j(n)B_{n+j} = binom(2n+2,n+1)(lam2/4)^(n+1)
\\ (b) the closed form of the source: Phi = eta * x F^2 sqrt(1-lam1 x) = eta * Phi_0/sqrt(1-lam2 x)
default(parisizemax, 24000000000);
read("lib.gp");
read("hosts.gp");
foldregC(N) = my(mf,d,Bb,ai,M); mf=mfinit([N,4],1); d=mfdim(mf); if(d==0, return([0,0])); Bb=mfbasis(mf); ai=mfatkininit(mf,N); M=matrix(d,d); for(j=1,d, my(vv=mftobasis(mf,mfatkin(ai,Bb[j]))); for(t=1,d, M[t,j]=vv[t])); [mf, matker(M+1)];
NQ = 160; NA = 150;
doh(h) = my(N,C,B,dv,r,tag,us,Fs,xs,a,fr,mf,V,g,cv,Th,b,ker,P,res,lam1,lam2,Phi0,S1,rhs,ok,d1,d2,Ph); N=h[1]; C=h[2]; B=h[3]; dv=h[4]; r=h[5]; tag=h[6]; fr=foldregC(N); if(fr[1]==0, return(0)); mf=fr[1]; V=fr[2]; if(matsize(V)[2]==0, return(0)); lam1=B+2*sqrtint(C); lam2=B-2*sqrtint(C); print(""); print("======== ", tag, "  lam1=", lam1, " lam2=", lam2); us=useries(dv,r,NQ); Fs=Fseries(dv,r,NQ); xs=us/(1+B*us+C*us^2); a=peel2(Fs,xs,NA,NQ); ker=fitrecR(a,2,3); P=normPR(vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); Phi0 = Fs*q*deriv(xs,q); S1 = sqrt(1 - lam1*xs + O(q^NQ)); for(j=1,matsize(V)[2], g=mflinear(mf, primvec(V[,j]~)~); cv=mfcoefs(g,NQ-1); if(cv[2]<0, cv=-cv); Th=thetaser(cv,3,NQ); b=peel2(Fs*Th,xs,NA,NQ); ok=1; for(n=0,NA-3, res = subst(P[1],x,n)*b[n+1]+subst(P[2],x,n)*b[n+2]+subst(P[3],x,n)*b[n+3]; rhs = binomial(2*n+2,n+1)*(lam2/4)^(n+1); if(res != rhs && res != -rhs, ok=0; break)); print("-- Phi[",j,"]  exact RHS binom(2n+2,n+1)(lam2/4)^(n+1) for n<=", NA-3, " ?  ", if(ok,"YES","no")); Ph = sum(i=1,NQ-1, cv[i+1]*q^i) + O(q^NQ); d1 = Ph - xs*Fs^2*S1; d2 = Ph + xs*Fs^2*S1; print("   Phi - x F^2 sqrt(1-lam1 x)  vanishes to O(q^",NQ,") ? ", if(d1==O(q^NQ),"YES", concat("no, first nonzero at q^", valuation(d1,q)))); print("   Phi + x F^2 sqrt(1-lam1 x)  vanishes to O(q^",NQ,") ? ", if(d2==O(q^NQ),"YES", concat("no, first nonzero at q^", valuation(d2,q)))));
for(i=1,#HOSTS, doh(HOSTS[i]));
quit;
