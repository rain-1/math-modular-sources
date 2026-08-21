/* Two level-16 zeta(5) rows from the same source C^full_{5,2}E_6 on two different
   host coordinates:
      row I  : t = x/(8x^2+2x+1)   (Fricke quotient; the canonical row)
      row II : y = x/(1+4x)        (Mobius shift of the plain Gamma_0(16) Hauptmodul)
   Both have archimedean limit (217/1024) zeta(5).  Test: same 2-adic limit
   (7/32) zeta_2(5)?  and does the cross determinant v_2 grow linearly (Conjecture D)?
   Log: lattice/zeta5_two_row/two_row16.log */
default(parisizemax, 10000000000);
\p 900
LOG = "lattice/zeta5_two_row/two_row16.log";
W(s) = write(LOG, s);
if(type(MM)!="t_INT", MM = 250);
M = MM; dd(n)=if(n<1,1,lcm(vector(n,i,i)));
q='q; ee(d)=eta(q^d+O(q^(M+2)));
X = q*ee(2)*ee(16)^2/(ee(1)^2*ee(8)) + O(q^(M+1));
cv=[1,-85,1428,-5440,4096]; dv=[1,2,4,8,16];
g=vector(M); for(n=1,M, s=0; for(i=1,#dv, if(n%dv[i]==0, s+=cv[i]*sigma(n/dv[i],5))); g[n]=s);
PHI=sum(n=1,M,g[n]*q^n)+O(q^(M+1)); THETA=sum(n=1,M,(g[n]/n^5)*q^n)+O(q^(M+1));

mk(T) = my(DT=q*deriv(T,q), F=PHI/DT, QT=serreverse(T), A=subst(F,q,QT), B=subst(F*THETA,q,QT)); \
        [vector(M,i,polcoeff(A,i-1)), vector(M,i,polcoeff(B,i-1))];
R1 = mk(X/(8*X^2+2*X+1));  A1 = R1[1]; B1 = R1[2];
R2 = mk(X/(1+4*X));        A2 = R2[1]; B2 = R2[2];

/* 2-adic zeta */
om(a)=if(a%4==1,1,-1); tw(a)=a/om(a);
KTR = 300; PR = 3*KTR;
z2 = (1/(4*(5-1))) * sum(i=1,2, my(a=[1,3][i]); (tw(a)+O(2^PR))^(-4) * \
        sum(j=0,KTR, binomial(-4,j)*((4/a)+O(2^PR))^j*bernfrac(j)));
xistar = (7/32)*z2;
W(Str("=== two level-16 zeta(5) rows, M=", M, " ==="));
W(Str("zeta_2(5): v_2 = ", valuation(z2,2), ";  predicted xi_2 = (7/32) zeta_2(5)"));

IDX=[25,50,100,150,200,M-1];
L = 217/1024*zeta(5);
rep(tag, AA, BB) = { W(""); W(Str("#### row ", tag)); \
  W(Str("  A_0..A_6 = ", vector(7,i,AA[i]))); \
  W(Str("  A_n integral: ", sum(n=0,M-1, if(denominator(AA[n+1])!=1,1,0))==0, \
        "   d_n^5 B_n integral: ", sum(n=0,M-1, if(denominator(dd(n)^5*BB[n+1])!=1,1,0))==0)); \
  W(Str("  lindep([B_N/A_N, zeta(5), 1]) = ", lindep([BB[M]/AA[M]*1.0, zeta(5), 1], 60))); \
  for(i=1,#IDX, my(n=IDX[i]); if(n<=M-1 && AA[n+1]!=0, \
    W(Str("    n=",n,"  log|A|/n=", log(abs(AA[n+1]*1.0))/n, \
          "   log|B-LA|/n=", log(abs(BB[n+1]*1.0-L*AA[n+1]))/n, \
          "   v2(A_n)=", valuation(AA[n+1],2), \
          "   v2(B/A-(7/32)zeta_2(5))=", valuation(BB[n+1]/AA[n+1]-xistar,2), \
          "   slope incr=", if(n>1&&AA[n]!=0, valuation(BB[n+1]/AA[n+1]-BB[n]/AA[n],2), 0))))); 0; };
rep("I   t = x/(8x^2+2x+1)", A1, B1);
rep("II  y = x/(1+4x)", A2, B2);
W("");
W("--- Conjecture D cross determinant  Delta_{n,m} = A^I_n B^II_m - A^II_m B^I_n  (both limits (217/1024)zeta(5)) ---");
for(i=1,#IDX, my(n=IDX[i], s=""); if(n<=M-1, \
  for(j=1,#IDX, my(m=IDX[j]); if(m<=M-1, s=Str(s,"  m=",m,":",valuation(A1[n+1]*B2[m+1]-A2[m+1]*B1[n+1],2)))); \
  W(Str("  n=",n,s))));
W(""); W("DONE"); quit;
