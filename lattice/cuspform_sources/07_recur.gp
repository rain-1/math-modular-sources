\\ Does the cusp-form companion B_n satisfy the host's own three-term recurrence?
\\ Fit the row recurrence, apply it to B, and report the inhomogeneity.
default(parisizemax, 24000000000);
default(realprecision, 60);
read("lib.gp");
read("hosts.gp");
foldregC(N) = my(mf,d,Bb,ai,M,Km); mf=mfinit([N,4],1); d=mfdim(mf); if(d==0, return([0,0])); Bb=mfbasis(mf); ai=mfatkininit(mf,N); M=matrix(d,d); for(j=1,d, my(vv=mftobasis(mf,mfatkin(ai,Bb[j]))); for(t=1,d, M[t,j]=vv[t])); [mf, matker(M+1)];
NQ = 200; NA = 180;
doh(h) = my(N,C,B,dv,r,tag,us,Fs,xs,a,fr,mf,V,g,cv,Th,b,ker,P,res,Phi0,cv0,Th0,b0); N=h[1]; C=h[2]; B=h[3]; dv=h[4]; r=h[5]; tag=h[6]; fr=foldregC(N); if(fr[1]==0, return(0)); mf=fr[1]; V=fr[2]; if(matsize(V)[2]==0, return(0)); print(""); print("======== ", tag); us=useries(dv,r,NQ); Fs=Fseries(dv,r,NQ); xs=us/(1+B*us+C*us^2); a=peel2(Fs,xs,NA,NQ); ker=fitrecR(a,2,3); P=normPR(vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); print("   P0 = ", P[1]); print("   P1 = ", P[2]); print("   P2 = ", P[3]); print("   charpoly = ", sum(j=1,3, polcoeff(P[j],3)*y^(j-1))); for(j=1,matsize(V)[2], g=mflinear(mf, primvec(V[,j]~)~); cv=mfcoefs(g,NQ-1); if(cv[2]<0, cv=-cv); Th=thetaser(cv,3,NQ); b=peel2(Fs*Th,xs,NA,NQ); res=vector(NA-2, n, subst(P[1],x,n-1)*b[n]+subst(P[2],x,n-1)*b[n+1]+subst(P[3],x,n-1)*b[n+2]); print("-- Phi[",j,"] residual of the host 3-term recurrence on B, n=0..12: ", vector(13,i,res[i])); my(lastnz=0); for(i=1,#res, if(res[i]!=0, lastnz=i)); print("   last nonzero residual at n = ", lastnz-1, "  (out of ", #res-1, ")"));
for(i=1,#HOSTS, doh(HOSTS[i]));
quit;
