/* Positive-P optimisation.  P >= 0 on (0,oo)  <=>  P = sigma0(w) + w*sigma1(w), sigma_i SOS.
   I_n(P) = Q^T H Q  (H = Hankel [M_{i+k}])  or  Q^T H1 Q  (H1 = [M_{i+k+1}]).
   With integer Q the linear form is  q G - p  with  q = den * A_P,  den = lcm_{j<=2d} den_j.
   Report  (1/n) log( den * min_{Q in Z^{d+1}\0} Q^T H Q ). */
\p 400
qf(H,q)=sum(i=1,#q,sum(k=1,#q,q[i]*H[i,k]*q[k]));
G=Catalan;
{
run(n) = my(m=3*n, A=vector(m+1), B=vector(m+1), Mv=vector(m+1), dn=vector(m+1));
 for(j=0,m, my(r=nestgen(m,j)); A[j+1]=r[1]; B[j+1]=r[2];
   Mv[j+1]=r[1]*G-r[2]; dn[j+1]=lcm(denominator(r[1]),denominator(r[2])));
 print("\n=== n=",n," m=",m,"  (1/n)log(den_j*M_j) min over j = ",
   vecmin(vector(m+1,j,log(1.0*dn[j]*abs(Mv[j]))/n)));
 print(" d | shift | (1/n)log det H | (1/n)log(min Q^tHQ) | (1/n)log(lcm den) | (1/n)log(den*val) | Q");
 for(sh=0,1,
  for(d=0,(m-sh)\2,
   my(H=matrix(d+1,d+1,i,k,Mv[i+k-2+sh+1]));
   my(sc=10^300, Hi=matrix(d+1,d+1,i,k,round(H[i,k]*sc/Mv[1+sh])));
   my(L=qflllgram(Hi));
   my(best=0,bq=0);
   for(c=1,d+1, my(q=vectorv(d+1,i,L[i,c]), v=qf(H,q));
     if(v>0 && (best==0||v<best), best=v; bq=q));
   /* also try small combos of the two shortest */
   if(d>=1, for(u=-3,3, for(t=-3,3, if(u||t, my(q=vectorv(d+1,i,u*L[i,1]+t*L[i,2]), v=qf(H,q));
      if(v>0&&v<best, best=v; bq=q)))));
   my(dd=lcm(vector(2*d+1,i,dn[i+sh])));
   printf("%3d | %5d | %10.4f | %12.4f | %12.4f | %12.4f | %s\n",
     d, sh, log(abs(matdet(H)))/n, log(best)/n, log(1.0*dd)/n, log(1.0*dd*best)/n, bq~);
  ));
}
for(n=1,5, run(n));
\q
