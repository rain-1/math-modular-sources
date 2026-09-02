\\ B_n is holonomic but not a solution of the host's 3-term recurrence:
\\ find the minimal recurrence, and identify the inhomogeneity exactly.
default(parisizemax, 24000000000);
read("lib.gp");
read("hosts.gp");
foldregC(N) = my(mf,d,Bb,ai,M); mf=mfinit([N,4],1); d=mfdim(mf); if(d==0, return([0,0])); Bb=mfbasis(mf); ai=mfatkininit(mf,N); M=matrix(d,d); for(j=1,d, my(vv=mftobasis(mf,mfatkin(ai,Bb[j]))); for(t=1,d, M[t,j]=vv[t])); [mf, matker(M+1)];
NQ = 200; NA = 180;
findrec(b) = my(ker, ord, dg, P); for(ord=3,5, for(dg=3,8, ker=fitrecR(b,ord,dg); if(matsize(ker)[2]==1, P=normPR(vector(ord+1,j,sum(e=0,dg,ker[(j-1)*(dg+1)+e+1,1]*x^e))); print("   minimal recurrence: order ", ord, " degree ", dg); for(j=1,ord+1, print("     Q",j-1," = ", P[j])); print("     char poly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); return(P)))); print("   no recurrence found with order<=5, degree<=8"); 0;
doh(h) = my(N,C,B,dv,r,tag,us,Fs,xs,a,fr,mf,V,g,cv,Th,b,ker,P,res,lam2,pred); N=h[1]; C=h[2]; B=h[3]; dv=h[4]; r=h[5]; tag=h[6]; fr=foldregC(N); if(fr[1]==0, return(0)); mf=fr[1]; V=fr[2]; if(matsize(V)[2]==0, return(0)); lam2=B-2*sqrtint(C); print(""); print("======== ", tag, "  lambda2 = ", lam2); us=useries(dv,r,NQ); Fs=Fseries(dv,r,NQ); xs=us/(1+B*us+C*us^2); a=peel2(Fs,xs,NA,NQ); ker=fitrecR(a,2,3); P=normPR(vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); for(j=1,matsize(V)[2], g=mflinear(mf, primvec(V[,j]~)~); cv=mfcoefs(g,NQ-1); if(cv[2]<0, cv=-cv); Th=thetaser(cv,3,NQ); b=peel2(Fs*Th,xs,NA,NQ); res=vector(20, n, subst(P[1],x,n-1)*b[n]+subst(P[2],x,n-1)*b[n+1]+subst(P[3],x,n-1)*b[n+2]); pred=vector(20, n, binomial(2*n,n)*lam2^n/2^n*if(1,1,1)); print("-- Phi[",j,"]  residual  = ", res); print("            c*[x^n](1-lam2 x)^(-1/2), n>=1 : ", vector(20,n, binomial(2*n,n)*(lam2/4)^n)); print("            ratio res_n / that   : ", vector(8,n, res[n]*1.0/(binomial(2*n,n)*(lam2/4.0)^n))); findrec(b));
for(i=1,#HOSTS, doh(HOSTS[i]));
quit;
