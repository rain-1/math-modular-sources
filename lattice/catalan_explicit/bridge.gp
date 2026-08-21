/* The 2-adic bridge across the whole moment family.
   For a pair S={j0,j1}: clear denominators with d=lcm, a_i=d A_i, b_i=d B_i in Z.
   The maximal extra divisibility available to a combination is governed by
   h = a_{j0} b_{j1} - a_{j1} b_{j0}  (the mixed minor, cf. CATALAN_AUDIT 4(d) v_2(h_n)~24.06n).
   Half-log-covolume predictor for the cone minimum (CATALAN_POSITIVITY 3.2):
      F = (1/2n)[ log|M_{j0}| + log|M_{j1}| + 2 log d - v2(h) log 2 ] . */
\p 400
G=Catalan;
{
run(n)= my(m=3*n, A=vector(m+1),B=vector(m+1),Mv=vector(m+1));
 for(j=0,m, my(r=nestgen(m,j)); A[j+1]=r[1];B[j+1]=r[2];Mv[j+1]=r[1]*G-r[2]);
 print("\n=== n=",n,"  v2(minor)/n   [rows j0, cols j1>j0] ===");
 my(bF=1e9,bp=0);
 for(j0=0,m-1, my(row=vector(m-j0), rowF=vector(m-j0));
  for(j1=j0+1,m,
   my(d=lcm(lcm(denominator(A[j0+1]),denominator(B[j0+1])),
            lcm(denominator(A[j1+1]),denominator(B[j1+1]))));
   my(h=(d*A[j0+1])*(d*B[j1+1])-(d*A[j1+1])*(d*B[j0+1]));
   my(v=valuation(h,2), F=(log(abs(Mv[j0+1]))+log(abs(Mv[j1+1]))+2*log(1.0*d)-v*log(2))/(2*n));
   row[j1-j0]=round(v*1.0/n*100)/100; rowF[j1-j0]=round(F*1000)/1000;
   if(F<bF, bF=F; bp=[j0,j1,v,d]));
  print("  j0=",j0," v2/n=",row);
  print("        F   =",rowF));
 print("  BEST F pair (",bp[1],",",bp[2],"): F=",bF,"  v2(h)/n=",bp[3]*1.0/n);
 my(d=lcm(lcm(denominator(A[1]),denominator(B[1])),lcm(denominator(A[n+1]),denominator(B[n+1]))));
 my(h=(d*A[1])*(d*B[n+1])-(d*A[n+1])*(d*B[1]));
 print("  Zudilin/Nesterenko pair (0,",n,"): v2(h)/n=",valuation(h,2)*1.0/n,
       "  F=",(log(abs(Mv[1]))+log(abs(Mv[n+1]))+2*log(1.0*d)-valuation(h,2)*log(2))/(2*n));
}
for(n=2,5, run(n));
\q
