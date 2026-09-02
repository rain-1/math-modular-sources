\\ Structural checks and the arithmetic budget.
\\ (a) Lemma: any W_N-antiinvariant weight-2 F has F'(tau_c) = N tau_c F(tau_c) at tau_c = i/sqrt N,
\\     so the Apery limit xi = Theta(tau_c) + Theta'(tau_c)/(N tau_c) does not depend on the host at all.
\\ (b) hence a rapidly convergent series for L(Phi,3) when Phi|W_N = -Phi.
\\ (c) K_- identification attempts; (d) the irrationality budget log|lambda_2| + 3.
default(parisizemax, 24000000000);
default(realprecision, 80);
read("lib.gp");
read("hosts.gp");
NQ = 300;
print("### (a) F'(tau_c) = N tau_c F(tau_c), i.e. DF(q_c) = (sqrt(N)/(2 Pi)) F(q_c)");
for(i=1,#HOSTS, my(h=HOSTS[i], N=h[1], dv=h[4], r=h[5], Fs, qc, Fc, DFc); Fs=Fseries(dv,r,NQ); qc=exp(-2*Pi/sqrt(N*1.0)); Fc=subst(truncate(Fs),q,qc); DFc=subst(truncate(q*deriv(Fs,q)),q,qc); print("   ", h[6], "   DF/F = ", DFc/Fc, "   sqrt(N)/(2Pi) = ", sqrt(N*1.0)/(2*Pi)));
print("");
print("### (b) L(Phi,3) = Theta(q_c) + (2Pi/sqrt N) DTheta(q_c),  q_c = exp(-2Pi/sqrt N),  for Phi|W_N = -Phi");
for(i=1,#[5,6,7,8,9,10,12,18], my(lv=[5,6,7,8,9,10,12,18], N=lv[i], mf, d, Bb, ai, M, Km, g, cv, Th, qc, val, L3); mf=mfinit([N,4],1); d=mfdim(mf); if(d==0, next); Bb=mfbasis(mf); ai=mfatkininit(mf,N); M=matrix(d,d); for(j=1,d, my(vv=mftobasis(mf,mfatkin(ai,Bb[j]))); for(t=1,d, M[t,j]=vv[t])); Km=matker(M+1); for(j=1,matsize(Km)[2], g=mflinear(mf, primvec(Km[,j]~)~); cv=mfcoefs(g,NQ-1); if(cv[2]<0, cv=-cv; g=mflinear(mf,-primvec(Km[,j]~)~)); Th=thetaser(cv,3,NQ); qc=exp(-2*Pi/sqrt(N*1.0)); val = subst(truncate(Th),q,qc) + (2*Pi/sqrt(N*1.0))*subst(truncate(q*deriv(Th,q)),q,qc); L3 = lfun(lfunmf(mf,g),3); print("   N=",N," Phi=",vector(6,t,cv[t]),"   series = ", val); print("                                        L(Phi,3) = ", L3, "   diff = ", abs(val-L3))));
print("");
print("### (d) irrationality budget: need log|lambda_2| + 3 < 0 (k=3 denominators d_n^3)");
for(i=1,#HOSTS, my(h=HOSTS[i], N=h[1], C=h[2], B=h[3], lam2=h[3]-2*sqrt(h[2]*1.0)); print("   ", h[6], "   lambda_2 = ", lam2, "   log|lambda_2|+3 = ", log(abs(lam2))+3));
quit;
