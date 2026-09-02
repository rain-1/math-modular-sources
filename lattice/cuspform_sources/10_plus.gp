\\ The W_N = +1 cusp forms: NOT fold-regular.
\\ (i) the Apery limit still exists but the convergence is only O(1/n);
\\ (ii) it equals the fold formula xi_gen = Theta(tau_c) + F(tau_c) DTheta(tau_c)/DF(tau_c)
\\     evaluated at the near fold tau_c = i/sqrt(N)  --  and NOT L(Phi,3);
\\ (iii) so the obstruction is the whole anti-invariant part of Theta, not a scalar.
default(parisizemax, 24000000000);
default(realprecision, 60);
read("lib.gp");
read("hosts.gp");
NQ = 220; NA = 200;
richfit(rv, n0, m) = my(M, V); M = matrix(m, m); for(i=1,m, my(n=n0+i-1); for(j=1,m, M[i,j] = 1.0/n^(j-1))); V = vector(m, i, rv[i])~; (matsolve(M, V))[1];
Wmat(mf, N) = my(d, Bb, ai, M); d=mfdim(mf); Bb=mfbasis(mf); ai=mfatkininit(mf,N); M=matrix(d,d); for(j=1,d, my(vv=mftobasis(mf,mfatkin(ai,Bb[j]))); for(t=1,d, M[t,j]=vv[t])); M;
doh(h) = my(N,C,B,dv,r,tag,us,Fs,xs,a,mf,d,M,Kp,Km,g,cv,Th,b,qc,xc,lam1,lam2,Fc,DFc,Tc,DTc,xigen,rv,xifit,e1,e2,K,nm,sg); N=h[1]; C=h[2]; B=h[3]; dv=h[4]; r=h[5]; tag=h[6]; mf=mfinit([N,4],1); d=mfdim(mf); if(d==0, return(0)); lam1=B+2*sqrt(C*1.0); lam2=B-2*sqrt(C*1.0); print(""); print("======== ", tag, "  lam1=", lam1, " lam2=", lam2); us=useries(dv,r,NQ); Fs=Fseries(dv,r,NQ); xs=us/(1+B*us+C*us^2); a=peel2(Fs,xs,NA,NQ); qc = exp(-2*Pi/sqrt(N*1.0)); xc = subst(truncate(xs), q, qc); print("   x(i/sqrt N) = ", xc, "   1/lam1 = ", 1.0/lam1, "   1/lam2 = ", 1.0/lam2); M=Wmat(mf,N); Km=matker(M+1); Kp=matker(M-1); Fc = subst(truncate(Fs),q,qc); DFc = subst(truncate(q*deriv(Fs,q)),q,qc); for(sg=1,2, K=if(sg==1,Km,Kp); nm=if(sg==1,"MINUS","PLUS "); for(j=1,matsize(K)[2], g=mflinear(mf, primvec(K[,j]~)~); cv=mfcoefs(g,NQ-1); if(cv[2]<0, cv=-cv; g=mflinear(mf, -primvec(K[,j]~)~)); Th=thetaser(cv,3,NQ); b=peel2(Fs*Th,xs,NA,NQ); Tc = subst(truncate(Th),q,qc); DTc = subst(truncate(q*deriv(Th,q)),q,qc); xigen = Tc + Fc*DTc/DFc; rv = vector(10, i, b[NA-9+i]*1.0/a[NA-9+i]); xifit = richfit(rv, NA-9, 10); e1 = rv[10]-rv[9]; e2 = rv[9]-rv[8]; print("-- ", nm, "[",j,"] Phi = ", vector(8,i,cv[i])); print("     B_n/A_n at n=", NA, "        = ", rv[10]); print("     Richardson limit         = ", xifit); print("     fold formula at tau_c    = ", xigen); print("     |Richardson - fold|      = ", abs(xifit-xigen)); print("     L(Phi,3)                 = ", lfun(lfunmf(mf,g),3)); if(e2==0, print("     converged to full precision (fold-regular)"), print("     log|d_n|/log|d_(n-1)| ratio (1 => power law) : ", e1/e2, "   expected (n/(n-1))^-2 = ", ((NA-1.0)/NA)^2))));
for(i=1,#HOSTS, doh(HOSTS[i]));
quit;
