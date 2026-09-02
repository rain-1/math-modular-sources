\\ For every Fricke host and every W_N-eigenvector in S_4(Gamma_0(N)):
\\ exact A_n = [x^n]F, B_n = [x^n](F * D^{-3}Phi), denominators, and the ratio B_n/A_n.
default(parisizemax, 16000000000);
default(realprecision, 80);
read("lib.gp");
read("hosts.gp");
NQ = 200; NA = 180;
dowork(h) = my(N, C, B, dv, r, us, Fs, xs, a, mf, d, Bb, ai, M, Kp, Km, g, cv, Th, b, tag, sg, K, nm); N=h[1]; C=h[2]; B=h[3]; dv=h[4]; r=h[5]; tag=h[6]; print(""); print("################ HOST ", tag, "  N=",N," C=",C," B=",B); us=useries(dv,r,NQ); Fs=Fseries(dv,r,NQ); xs=us/(1+B*us+C*us^2); a=peel2(Fs,xs,NA,NQ); print("A = ", vector(8,i,a[i])); mf=mfinit([N,4],1); d=mfdim(mf); print("dim S_4(Gamma_0(",N,")) = ", d); if(d==0, print("  no cusp forms"); return(0)); Bb=mfbasis(mf); ai=mfatkininit(mf,N); M=matrix(d,d); for(j=1,d, my(vv=mftobasis(mf, mfatkin(ai,Bb[j]))); for(t=1,d, M[t,j]=vv[t])); Km=matker(M+1); Kp=matker(M-1); print("dim(W=-1) = ", matsize(Km)[2], "  dim(W=+1) = ", matsize(Kp)[2]); for(sg=1,2, K=if(sg==1,Km,Kp); nm=if(sg==1,"MINUS","PLUS "); for(j=1,matsize(K)[2], g=mflinear(mf, primvec(K[,j]~)~); cv=mfcoefs(g, NQ-1); if(cv[2]<0, cv=-cv); print("--- ", nm, "[",j,"] eps=", if(sg==1,-1,1), "  Phi = ", vector(13,i,cv[i])); Th=thetaser(cv,3,NQ); b=peel2(Fs*Th,xs,NA,NQ); print("    B = ", vector(6,i,b[i])); print("    k(den) = ", denexp(b,NA)); print("    B/A at n=", NA, " : ", b[NA+1]*1.0/a[NA+1]); print("    B/A at n=", NA-1, " : ", b[NA]*1.0/a[NA]); print("    B/A at n=", NA-2, " : ", b[NA-1]*1.0/a[NA-1]); write(concat(concat("data_", tag), concat(concat("_", nm), concat(j, ".txt"))), "A=", a); write(concat(concat("data_", tag), concat(concat("_", nm), concat(j, ".txt"))), "B=", b))); 1;
for(i=1,#HOSTS, dowork(HOSTS[i]));
quit;
