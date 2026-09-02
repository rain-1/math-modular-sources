\\ Apery periods of every holomorphic fold-regular source on every Fricke host,
\\ read off the row (B_n/A_n converges like (lambda2/lambda1)^n) and identified by lindep.
default(parisizemax, 24000000000);
default(realprecision, 400);
read("lib.gp");
read("hosts.gp");
Z3 = zeta(3);
Lm3 = lfun(-3, 3);
Lm4 = lfun(-4, 3);
Lf5 = lfun(lfunmf(mfinit([5,4],0), mfeigenbasis(mfinit([5,4],0))[1]), 3);
Lf6 = lfun(lfunmf(mfinit([6,4],0), mfeigenbasis(mfinit([6,4],0))[1]), 3);
Lf9 = lfun(lfunmf(mfinit([9,4],0), mfeigenbasis(mfinit([9,4],0))[1]), 3);
BASIS = [Z3, Lm3, Lm4, Lf5, Lf6, Lf9];
NAM = ["zeta(3)", "L(3,chi-3)", "L(3,chi-4)", "L(5.4.a.a,3)", "L(6.4.a.a,3)", "L(9.4.a.a,3)"];
foldreg(N) = my(mf,d,Bb,ai,M,Km,rows,W); mf=mfinit([N,4],4); d=mfdim(mf); Bb=mfbasis(mf); ai=mfatkininit(mf,N); M=matrix(d,d); for(j=1,d, my(vv=mftobasis(mf,mfatkin(ai,Bb[j]))); for(t=1,d, M[t,j]=vv[t])); Km=matker(M+1); if(matsize(Km)[2]==0, return([mf, matrix(d,0)])); rows=vector(matsize(Km)[2], j, mfcoefs(mflinear(mf,Km[,j]),0)[1]); W=matker(Mat(rows)); [mf, Km*W];
doh(h, NQ, NA) = my(N,C,B,dv,r,tag,us,Fs,xs,a,fr,mf,V,g,cv,Th,b,xi,rel,s); N=h[1]; C=h[2]; B=h[3]; dv=h[4]; r=h[5]; tag=h[6]; fr=foldreg(N); mf=fr[1]; V=fr[2]; print(""); print("======== ", tag, "   dim V = ", matsize(V)[2], "   lam2/lam1 = ", (B-2*sqrt(C*1.0))/(B+2*sqrt(C*1.0))); if(matsize(V)[2]==0, print("   no holomorphic fold-regular source"); return(0)); us=useries(dv,r,NQ); Fs=Fseries(dv,r,NQ); xs=us/(1+B*us+C*us^2); a=peel2(Fs,xs,NA,NQ); for(j=1,matsize(V)[2], g=mflinear(mf, primvec(V[,j]~)~); cv=mfcoefs(g,NQ-1); if(cv[2]<0, cv=-cv); Th=thetaser(cv,3,NQ); b=peel2(Fs*Th,xs,NA,NQ); xi=b[NA+1]*1.0/a[NA+1]; s=-log(abs(b[NA+1]/a[NA+1]/(b[NA]/a[NA])-1))/log(10.); print("-- V[",j,"] Phi = ", vector(10,i,cv[i]), "   k(den) = ", denexp(b,NA)); print("   xi = ", xi); print("   stable digits ~ ", s); rel=lindep(concat([xi], BASIS), floor(s*0.7)); print("   lindep  = ", rel));
for(i=1,#HOSTS, doh(HOSTS[i], 300, 280));
quit;
