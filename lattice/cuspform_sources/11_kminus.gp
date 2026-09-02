\\ The linear form w_n = xi*A_n - B_n for the fold-regular cusp sources:
\\ rate lambda_2^n n^(-3/2), constant K_- (Richardson), and the Casoratian test.
default(parisizemax, 30000000000);
default(realprecision, 900);
read("lib.gp");
read("hosts.gp");
NQ = 140; NA = 130; NIT = 800;
foldregC(N) = my(mf,d,Bb,ai,M); mf=mfinit([N,4],1); d=mfdim(mf); if(d==0, return([0,0])); Bb=mfbasis(mf); ai=mfatkininit(mf,N); M=matrix(d,d); for(j=1,d, my(vv=mftobasis(mf,mfatkin(ai,Bb[j]))); for(t=1,d, M[t,j]=vv[t])); [mf, matker(M+1)];
findrec(b) = my(ker, ord, dg, P); for(ord=3,5, for(dg=3,8, ker=fitrecR(b,ord,dg); if(matsize(ker)[2]==1, P=normPR(vector(ord+1,j,sum(e=0,dg,ker[(j-1)*(dg+1)+e+1,1]*x^e))); return([ord,P])))); [0,0];
iterhom(P, ord, init, nb) = my(v=vector(nb+ord+1)); for(i=1,ord, v[i]=init[i]); for(n=0,nb, my(s=0, lead=subst(P[ord+1],x,n)); if(lead==0, break); for(j=1,ord, s += subst(P[j],x,n)*v[n+j]); v[n+ord+1] = -s/lead); v;
richfit(rv, n0, m) = my(M, V); M = matrix(m, m); for(i=1,m, my(n=n0+i-1); for(j=1,m, M[i,j] = 1.0/n^(j-1))); V = vector(m, i, rv[i])~; (matsolve(M, V))[1];
doh(h) = my(N,C,B,dv,r,tag,us,Fs,xs,a,fr,mf,V,g,cv,Th,b,ker,P,lam1,lam2,fr2,ord,Q,av,bv,xi,wn,Kv,Km,Kp,Kpv,kappa,d0); N=h[1]; C=h[2]; B=h[3]; dv=h[4]; r=h[5]; tag=h[6]; fr=foldregC(N); if(fr[1]==0, return(0)); mf=fr[1]; V=fr[2]; if(matsize(V)[2]==0, return(0)); lam1=B+2*sqrt(C*1.0); lam2=B-2*sqrt(C*1.0); print(""); print("======== ", tag, "  lam1=", lam1, "  lam2=", lam2); us=useries(dv,r,NQ); Fs=Fseries(dv,r,NQ); xs=us/(1+B*us+C*us^2); a=peel2(Fs,xs,NA,NQ); ker=fitrecR(a,2,3); P=normPR(vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); d0 = polcoeff(P[1],0)/polcoeff(P[1],3); kappa = if(polcoeff(P[1],3)==0, 0, my(cc=-polcoeff(P[1],3), dd=-polcoeff(P[1],0)+0*cc); 0); av = iterhom(P, 2, [a[1],a[2]], NIT+6); Kpv = vector(12, i, av[NIT-11+i]*1.0*(NIT-12+i)^(3/2)/lam1^(NIT-12+i)); Kp = richfit(Kpv, NIT-12, 12); print("   K_+ measured        = ", Kp); print("   K_+ closed form     = ", sqrt(N*1.0)/(2*Pi^(3/2))*sqrt(lam1/(lam1-lam2))); for(j=1,matsize(V)[2], g=mflinear(mf, primvec(V[,j]~)~); cv=mfcoefs(g,NQ-1); if(cv[2]<0, cv=-cv; g=mflinear(mf,-primvec(V[,j]~)~)); Th=thetaser(cv,3,NQ); b=peel2(Fs*Th,xs,NA,NQ); fr2=findrec(b); ord=fr2[1]; Q=fr2[2]; if(ord==0, print("-- Phi[",j,"] no recurrence found"); next); bv = iterhom(Q, ord, vector(ord,i,b[i]), NIT+6); my(okk=1); for(i=1,NA-2, if(bv[i]!=b[i], okk=0; break)); print("-- Phi[",j,"] order-", ord, " recurrence reproduces B_n for n<=", NA-3, " : ", if(okk,"YES","NO")); xi = lfun(lfunmf(mf,g),3); Kv = vector(12, i, my(n=NIT-12+i); (xi*av[n+1]-bv[n+1])*n^(3/2)/lam2^n); Km = richfit(Kv, NIT-12, 12); print("   K_- (Richardson)    = ", Km); print("   K_+ K_- (lam1-lam2) = ", Kp*Km*(lam1-lam2)));
for(i=1,#HOSTS, doh(HOSTS[i]));
quit;
