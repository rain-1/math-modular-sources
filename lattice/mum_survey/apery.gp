/* ---------------------------------------------------------------
   apery.gp -- Apery pairs of order-4 CY operators (AvSZ definition).
   Operator L = sum_{i=0}^{k} z^i P_i(theta).  Recurrence
        sum_{i=0}^{k} P_i(n-i) u_{n-i} = 0.
   A: A_0=1, A_n from the recurrence (n>=1).
   B: B_n=0 for n<=0, B_1=1, recurrence for n>=2.
   Apery limit  = lim B_n/A_n   (Almkvist-van Straten-Zudilin 2008).
   --------------------------------------------------------------- */

/* build both sequences to length N+1 */
aperyPair(pols, N) =
{ my(k=#pols-1, P=vector(k+1,i,pols[i]), A=vector(N+1), B=vector(N+1), lead);
  A[1]=1; B[1]=0;
  for(n=1,N,
    lead = subst(P[1],'X,n);
    if(lead==0, return(0));
    my(sa=0, sb=0);
    for(i=1,min(k,n),
      my(c=subst(P[i+1],'X,n-i));
      sa += c*A[n-i+1]; sb += c*B[n-i+1]);
    A[n+1] = -sa/lead;
    if(n==1, B[2]=1, B[n+1] = -sb/lead));
  [A,B];
}

/* denominator exponent: minimal k with d_n^k B_n in Z, tested up to N */
denomExp(B,N) =
{ my(kmax=0, dn=1);
  for(n=1,N, dn=lcm(dn,n);
    my(d=denominator(B[n+1]), e=0);
    while(d>1, my(g=gcd(d,dn)); if(g==1, return(-1)); d/=g; e++);
    if(e>kmax, kmax=e));
  kmax;
}

/* p-adic denominator rate of A: v_p(A_n) ~ -kappa_p n  */
kappaRate(A,p,N) = my(n1=N\2, n2=N); -(valuation(A[n2+1],p)-valuation(A[n1+1],p))/(n2-n1)*1.0;

/* p-adic slope: v_p(B_n/A_n - B_{n-1}/A_{n-1}) ~ sigma_p n */
slopeVal(A,B,p,n) = if(A[n+1]==0 || A[n]==0, 0, my(d=B[n+1]/A[n+1]-B[n]/A[n]); if(d==0, 0, valuation(d,p)));
slopeData(A,B,p,ns) = vector(#ns,i, slopeVal(A,B,p,ns[i]));
