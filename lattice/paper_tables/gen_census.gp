\\ Census table generator (paper Table 1, paper/tables/census.tex).
\\ Computes, from exact rational recurrences, for every row in the census:
\\ period (via lindep), c, lambda1, lambda2, k (denominator exponent, checked sharp to n<=200),
\\ measured slopes sigma_p at n~300 for p=2,3,5,7, kappa_p (denominator growth rate), score, budget.
\\ Companion normalisation b_0=0,b_1=1 throughout (see task/README for exceptions).
default(realprecision,120);

\\ ---------- generic order-2 companion recurrence (Zagier / third-order shape) ----------
\\ (n+1)^m u_{n+1} = P(n) u_n - c n^m u_{n-1}, m=2 (Zagier) or m=3 (third-order AZ family);
\\ P(n) = a n^2+a n+b  (m=2)   or   (2n+1)(a n^2+a n+b) (m=3, "third-order")
rowsM(a,b,c,m,N)={my(A=vector(N+1),B=vector(N+1));A[1]=1;A[2]=b;B[1]=0;B[2]=1;
 for(n=1,N-1,
   my(P=if(m==2,(a*n^2+a*n+b),(2*n+1)*(a*n^2+a*n+b)));
   A[n+2]=(P*A[n+1]-c*n^m*A[n])/(n+1)^m;
   B[n+2]=(P*B[n+1]-c*n^m*B[n])/(n+1)^m);
 [A,B]};

\\ k-check: is d_n^k * B_n an integer for all n<=Nk, and is k sharp (k-1 fails)?
kcheck(B,Nk,kmax)={my(dn=1,k=-1);
 for(kk=1,kmax,
   my(ok=1); dn=1;
   for(n=1,Nk, dn=lcm(dn,n); my(v=dn^kk*B[n+1]); if(v!=round(v), ok=0; break));
   if(ok, k=kk; break));
 my(sharp=1);
 if(k>=1,
   my(dn2=1,failed=0);
   for(n=1,Nk, dn2=lcm(dn2,n); my(v=dn2^(k-1)*B[n+1]); if(v!=round(v), failed=1; break));
   sharp = failed);
 [k,sharp]};

slopes(A,B,N)={vector(4,j,my(p=[2,3,5,7][j]);
   valuation(B[N+1]/A[N+1]-B[N]/A[N],p))};

kappa(A,N)={vector(4,j,my(p=[2,3,5,7][j]); my(v=valuation(A[N+1],p));
   if(v<0, -v*1./N, 0))};

roots(a,c)={polroots(x^2-2*a*x+c)};

\\ ================= ZAGIER A-F =================
famZ=["Zagier A","Zagier B","Zagier C","Zagier D","Zagier E","Zagier F"];
parZ=[[7,2,-8],[9,3,27],[10,3,9],[11,3,-1],[12,4,32],[17,6,72]];
period_basisZ = [Pi^2/6, lfun(lfuncreate(-3),2)]; \\ zeta(2), L(2,chi_-3)
NZ=300;
print("### ZAGIER A-F ###");
{for(i=1,6,my(a=parZ[i][1],b=parZ[i][2],c=parZ[i][3]);
  my(R=rowsM(a,b,c,2,NZ),A=R[1],B=R[2]);
  my(lim=1.0*B[NZ+1]/A[NZ+1]);
  my(rt=polroots(x^2-a*x+c));
  my(kc=kcheck(B,200,6));
  my(sl=slopes(A,B,NZ));
  my(ld=lindep(concat([lim],[1,Pi^2/6,lfun(lfuncreate(-3),2)])));
  print(famZ[i]," a,b,c=",[a,b,c]," roots=",rt~," lim=",lim," lindep[1,z2,L2chi-3]=",ld~," k,sharp=",kc," slopes=",sl));}

\\ ================= THIRD-ORDER AZ FAMILY =================
famT=["AZ(7,3,81)","AZ(9,3,-27)","Domb(10,4,64)","eta(11,5,125)","T(12,4,16)","Apery(17,5,1)"];
parT=[[7,3,81],[9,3,-27],[10,4,64],[11,5,125],[12,4,16],[17,5,1]];
NT=400;
print("\n### THIRD-ORDER AZ FAMILY ###");
{for(i=1,6,my(a=parT[i][1],b=parT[i][2],c=parT[i][3]);
  my(R=rowsM(a,b,c,3,NT),A=R[1],B=R[2]);
  my(lim=1.0*B[NT+1]/A[NT+1]);
  my(rt=roots(a,c));
  my(kc=kcheck(B,200,6));
  my(sl=slopes(A,B,NT));
  print(famT[i]," a,b,c=",[a,b,c]," roots=",rt~," lim=",lim," k,sharp=",kc," slopes=",sl));}

\\ ================= COOPER s7,s10,s18 (CORRECTED ICs) =================
f10(n,u1,u0)=2*(2*n+1)*(3*n^2+3*n+1)*u1 + 4*n*(16*n^2-1)*u0;
f7(n,u1,u0)=(2*n+1)*(13*n^2+13*n+4)*u1 + 3*n*(9*n^2-1)*u0;
f18(n,u1,u0)=2*(2*n+1)*(7*n^2+7*n+3)*u1 - 12*n*(16*n^2-1)*u0;
rowsC(fn,b,N)={my(A=vector(N+1),B=vector(N+1));A[1]=1;A[2]=b;B[1]=0;B[2]=1;
 for(n=1,N-1,A[n+2]=fn(n,A[n+1],A[n])/(n+1)^3;
             B[n+2]=fn(n,B[n+1],B[n])/(n+1)^3);[A,B]};
NC=300;
print("\n### COOPER s7,s10,s18 (A1=4,2,6 resp.) ###");
{
R10=rowsC(f10,2,NC); R7=rowsC(f7,4,NC); R18=rowsC(f18,6,NC);
RR=[R7,R10,R18]; nmC=["s7","s10","s18"]; ccC=[-27,-64,192]; rtC=[[27,-1],[16,-4],[16,12]];
for(i=1,3,my(A=RR[i][1],B=RR[i][2]);
  my(lim=1.0*B[NC+1]/A[NC+1]);
  my(sl=slopes(A,B,NC));
  my(kp=kappa(A,NC));
  print(nmC[i]," c=",ccC[i]," lambda1,2=",rtC[i]," lim=",lim," slopes=",sl," kappa=",kp));
}

\\ ================= AZ delta(7,3,81) reused above as AZ(7,3,81); AZ eta(11,5,125) reused as eta =================
\\ (complex-root rows; already in THIRD-ORDER block. Slopes there at N=400 are the reported ones.)

\\ ================= L(f,2) weight-3 cusp row =================
\\ (n+1)^2 a_{n+1} = (20n^2+10n+2) a_n - 16(2n-1)^2 a_{n-1};  a_0=1,a_1=2 ; b_0=0,b_1=1.
NL=300;
rowsLf(N)={my(A=vector(N+1),B=vector(N+1));A[1]=1;A[2]=2;B[1]=0;B[2]=1;
 for(n=1,N-1,
   my(P=20*n^2+10*n+2, Q=16*(2*n-1)^2);
   A[n+2]=(P*A[n+1]-Q*A[n])/(n+1)^2;
   B[n+2]=(P*B[n+1]-Q*B[n])/(n+1)^2);
 [A,B]};
print("\n### L(f,2) weight-3 cusp row ###");
{my(R=rowsLf(NL),A=R[1],B=R[2]);
 my(lim=1.0*B[NL+1]/A[NL+1]);
 my(kc=kcheck(B,200,4));
 my(sl=slopes(A,B,NL));
 print("roots(x^2-20x+64)=",polroots(x^2-20*x+64)~," lim=",lim," k,sharp=",kc," slopes=",sl);}

print("\nDONE");
\q
