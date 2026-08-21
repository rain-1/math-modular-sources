/* True-denominator search over SOS: P = Q(w)^2 (and w*Q^2), Q in Z^{d+1}, |q_i|<=BOX.
   A_P = sum q_i q_k A_{i+k},  B_P likewise; den = lcm(den A_P, den B_P).
   Objective  obj = (1/n) log( den * |A_P G - B_P| )   -- the linear form  den*A_P*G - den*B_P
   has INTEGER coefficients and is POSITIVE by construction (P>=0 on (0,oo)).
   Invariant under P -> cP for rational c, so integer Q is WLOG. */
\p 400
G=Catalan;
{
srch(n,DMAX,BOX) = my(m=3*n, A=vector(m+1), B=vector(m+1), Mv=vector(m+1));
 for(j=0,m, my(r=nestgen(m,j)); A[j+1]=r[1]; B[j+1]=r[2]; Mv[j+1]=r[1]*G-r[2]);
 print("\n=== n=",n," m=",m," BOX=",BOX);
 for(sh=0,1, for(d=0,min(DMAX,(m-sh)\2),
  my(best=1e9, bq=0, bden=0, bval=0, nq=(2*BOX+1)^(d+1));
  for(idx=0,nq-1,
   my(t=idx, q=vector(d+1)); for(i=1,d+1, q[i]=t%(2*BOX+1)-BOX; t=t\(2*BOX+1));
   if(q[d+1]<0, next); if(vecmax(vector(d+1,i,abs(q[i])))==0, next);
   if(content(vector(d+1,i,q[i]))!=1, next);
   my(AP=0,BP=0,VP=0);
   for(i=1,d+1, for(k=1,d+1, my(c=q[i]*q[k], j=i+k-2+sh);
      AP+=c*A[j+1]; BP+=c*B[j+1]; VP+=c*Mv[j+1]));
   if(VP<=0, next);
   my(dd=lcm(denominator(AP),denominator(BP)), o=log(1.0*dd*VP)/n);
   if(o<best, best=o; bq=q; bden=dd; bval=VP));
  printf("  d=%d sh=%d : best obj = %9.4f   (1/n)log den=%8.4f  (1/n)log val=%8.4f   Q=%s\n",
    d,sh,best,log(1.0*bden)/n,log(bval)/n,bq));
 );
}
srch(2,3,6); srch(3,3,6); srch(4,3,6); srch(5,3,5);
\q
