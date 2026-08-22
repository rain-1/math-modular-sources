/* lattice/positivity/signs.gp  --- Task 1: the SIGN census of the corpus.
   For every Apery-like row in the ledger we (i) verify the exact Casoratian
        a_n b_{n-1} - a_{n-1} b_n = -c^{n-1}/n^{w+1}       (w+1 = order)
   (ii) deduce the exact tail identity
        a_N xi - b_N = a_N * sum_{n>N} c^{n-1}/(n^{w+1} a_n a_{n-1}),
   whose summands have the constant sign of c^{n-1}; and (iii) check the
   observed sign of a_N xi - b_N against the prediction
        c>0, a_n>0  ==>  a_N xi - b_N > 0 for every N   (fixed-sign kernel free)
        c<0, a_n>0  ==>  (-1)^N (a_N xi - b_N) > 0 provided |terms| decrease.
   No integral is needed: the Casoratian IS the positive kernel.             */
\p 1500

/* order-2 Zagier normalisation */
{
row2(aa,bb,cc,NN) = my(A=vector(NN+1), B=vector(NN+1));
 A[1]=1; A[2]=bb; B[1]=0; B[2]=1;
 for(n=1,NN-1,
  A[n+2]=((aa*n^2+aa*n+bb)*A[n+1]-cc*n^2*A[n])/(n+1)^2;
  B[n+2]=((aa*n^2+aa*n+bb)*B[n+1]-cc*n^2*B[n])/(n+1)^2);
 [A,B];
}
/* order-3 Zagier/Almkvist-Zudilin normalisation */
{
row3(aa,bb,cc,NN) = my(A=vector(NN+1), B=vector(NN+1));
 A[1]=1; A[2]=bb; B[1]=0; B[2]=1;
 for(n=1,NN-1,
  A[n+2]=((2*n+1)*(aa*n^2+aa*n+bb)*A[n+1]-cc*n^3*A[n])/(n+1)^3;
  B[n+2]=((2*n+1)*(aa*n^2+aa*n+bb)*B[n+1]-cc*n^3*B[n])/(n+1)^3);
 [A,B];
}

{
census(nm, ord, aa, bb, cc, xi, NN) =
 my(r=if(ord==2, row2(aa,bb,cc,NN), row3(aa,bb,cc,NN)), A=r[1], B=r[2], w1=ord);
 /* (i) Casoratian check, exact rational arithmetic */
 my(casok=1);
 for(n=1,NN, if(A[n+1]*B[n] - A[n]*B[n+1] != -cc^(n-1)/n^w1, casok=0; break));
 /* positivity of the a_n */
 my(apos=1, azero=0);
 for(n=0,NN, if(A[n+1]==0, azero=1); if(A[n+1]<=0, apos=0));
 /* (iii) observed signs of the linear form */
 my(sg=vector(NN+1), allpos=1, altern=1, s0=sign(xi));
 for(n=0,NN, sg[n+1]=sign(A[n+1]*xi-B[n+1]));
 for(n=0,NN, if(sg[n+1]<=0, allpos=0); if(sg[n+1]!=(-1)^n*s0, altern=0));
 /* decay/growth rate of the form and of a_n */
 my(lf=log(abs(A[NN+1]*xi-B[NN+1]))/NN, la=log(abs(A[NN+1]))/NN);
 printf("%-14s ord=%d c=%6s  Casoratian:%s  a_n>0:%s  form: %s   n^-1 log|form|=%+8.4f  n^-1 log a_n=%+8.4f\n",
   nm, ord, cc, if(casok,"OK ","FAIL"), if(apos,"yes","NO "),
   if(allpos, "STRICTLY POSITIVE, all n<=NN",
     if(altern, "sign = (-1)^n, all n<=NN   ", "oscillating, no fixed pattern")),
   lf, la);
}

/* the tail identity, checked numerically */
{
tailcheck(nm, ord, aa, bb, cc, xi, NN, Ntest) =
 my(r=if(ord==2, row2(aa,bb,cc,NN), row3(aa,bb,cc,NN)), A=r[1], B=r[2], w1=ord);
 my(S=sum(n=Ntest+1, NN, cc^(n-1)/(n^w1*A[n+1]*A[n])));
 printf("   tail check %-12s N=%d:  a_N*xi-b_N = %.30g   a_N*sum = %.30g\n",
   nm, Ntest, 1.*(A[Ntest+1]*xi-B[Ntest+1]), 1.*(A[Ntest+1]*S));
}
