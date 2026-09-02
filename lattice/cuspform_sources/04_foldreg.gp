\\ The fold-regular cases: Phi in the W_N = -1 eigenspace of S_4(Gamma_0(N)).
\\ Verifies xi = L(Phi,3), the two functional-equation identities, the denominator power,
\\ the decay rate lambda_2^n n^(-3/2) and its constant K_-.
default(parisizemax, 16000000000);
default(realprecision, 250);
read("lib.gp");
read("hosts.gp");
NQ = 260; NA = 240;
doone(h) = my(N,C,B,dv,r,tag,us,Fs,xs,a,mf,d,Bb,ai,M,Km,g,cv,Th,b,L1,L2,L3,lam1,lam2,xi,w,Km2,kk,rat,Kminus,dig); N=h[1]; C=h[2]; B=h[3]; dv=h[4]; r=h[5]; tag=h[6]; mf=mfinit([N,4],1); d=mfdim(mf); if(d==0, return(0)); Bb=mfbasis(mf); ai=mfatkininit(mf,N); M=matrix(d,d); for(j=1,d, my(vv=mftobasis(mf,mfatkin(ai,Bb[j]))); for(t=1,d, M[t,j]=vv[t])); Km=matker(M+1); if(matsize(Km)[2]==0, return(0)); lam1=B+2*sqrt(C*1.0); lam2=B-2*sqrt(C*1.0); print(""); print("======== HOST ", tag, "   lambda1=", lam1, "  lambda2=", lam2); us=useries(dv,r,NQ); Fs=Fseries(dv,r,NQ); xs=us/(1+B*us+C*us^2); a=peel2(Fs,xs,NA,NQ); for(j=1,matsize(Km)[2], g=mflinear(mf, primvec(Km[,j]~)~); cv=mfcoefs(g,NQ-1); if(cv[2]<0, cv=-cv; g=mflinear(mf,-primvec(Km[,j]~)~)); print("-- Phi[",j,"] = ", vector(15,i,cv[i])); L1=lfun(lfunmf(mf,g),1); L2=lfun(lfunmf(mf,g),2); L3=lfun(lfunmf(mf,g),3); print("   L(Phi,1) = ", L1); print("   L(Phi,2) = ", L2, "   (must vanish)"); print("   L(Phi,3) = ", L3); print("   L(Phi,1) + N*L(Phi,3)/(2*Pi^2) = ", L1 + N*L3/(2*Pi^2), "   (must vanish)"); Th=thetaser(cv,3,NQ); b=peel2(Fs*Th,xs,NA,NQ); print("   B_1..6 = ", vector(6,i,b[i])); kk=denexp(b,NA); print("   k(den) = ", kk); xi=L3; w=vector(12,i, a[NA-11+i]*xi - b[NA-11+i]); print("   w_n/w_(n-1) near n=", NA, " : ", w[12]/w[11], " , ", w[11]/w[10], "   (lambda2 = ", lam2, ")"); print("   K_- = w_n n^(3/2)/lambda2^n : ", w[12]*(NA*1.0)^(3/2)/lam2^NA, " , ", w[11]*((NA-1)*1.0)^(3/2)/lam2^(NA-1)); dig=-log(abs(b[NA+1]/a[NA+1]/xi-1))/log(10.); print("   digits of B_n/A_n = L(Phi,3) at n=", NA, " : ", dig); print("   xi (40 digits) = ", xi));
for(i=1,#HOSTS, doone(HOSTS[i]));
quit;
