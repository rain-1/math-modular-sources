/* Canonical formulation, no arbitrary modulus.
   L_S = { p in Z^S : sum p_j A_j in Z  and  sum p_j B_j in Z }  (S a set of moment indices).
   Every p in L_S with p_j >= 0 gives  0 < q G - p' = I_n(P) with q,p' in Z, P = sum p_j w^j >= 0.
   Report  (1/n) log (cone-min)  over all PAIRS S={j0,j1}, and the 2-row Zudilin/Nesterenko
   pair S={0,n} for reference. */
\p 400
G=Catalan;
{ conemin2(B,l0,l1,R) = my(best=0,bv=0);
 for(i=-R,R, my(lo=-oo,hi=oo,ok=1);
  for(r=1,2, my(a=i*B[r,1],b=B[r,2]);
    if(b==0, if(a<0,ok=0), my(t=-a/b); if(b>0,if(t>lo,lo=t),if(t<hi,hi=t))));
  if(ok, my(j0=if(lo==-oo,-R,ceil(lo)), j1=if(hi==oo,R,floor(hi)));
   if(j0<=j1, for(w=0,1, my(j=if(w==0,j0,j1), c0=i*B[1,1]+j*B[1,2], c1=i*B[2,1]+j*B[2,2]);
     if(c0||c1, my(v=c0*l0+c1*l1; ); if(best==0||v<best,best=v;bv=[c0,c1]))))));
 [best,bv]; }
{
run(n,R)= my(m=3*n, A=vector(m+1),B=vector(m+1),Mv=vector(m+1));
 for(j=0,m, my(r=nestgen(m,j)); A[j+1]=r[1];B[j+1]=r[2];Mv[j+1]=r[1]*G-r[2]);
 print("\n=== n=",n," m=",m,"  canonical pair cone-minima, (1/n)log ===");
 my(tab=matrix(m+1,m+1), bestall=0, bp=0);
 for(j0=0,m-1, for(j1=j0+1,m,
   my(d=lcm(lcm(denominator(A[j0+1]),denominator(B[j0+1])),
            lcm(denominator(A[j1+1]),denominator(B[j1+1]))));
   my(a0=d*A[j0+1],a1=d*A[j1+1],b0=d*B[j0+1],b1=d*B[j1+1]);
   my(K=matkerint([a0,a1,d,0,0,0; b0,b1,0,d,0,0]), L=mathnf(matrix(2,#K[1,],i,j,K[i,j])));
   my(l0=Mv[j0+1], l1=Mv[j1+1]);
   my(sc=10^350, LL=L*qflll([round(L[1,1]*l0*sc),round(L[1,2]*l0*sc);round(L[2,1]*l1*sc),round(L[2,2]*l1*sc)]));
   my(r=conemin2(LL,l0,l1,R), v=log(r[1])/n);
   tab[j0+1,j1+1]=v;
   if(bestall==0||v<bestall, bestall=v; bp=[j0,j1,r[2]]);
 ));
 for(j0=0,m-1, print("  j0=",j0,": ",vector(m-j0,i,round(tab[j0+1,j0+i+1]*1000)/1000)));
 print("  BEST pair: j=(",bp[1],",",bp[2],")  (1/n)log cone-min = ",bestall,"   coeffs=",bp[3]);
 print("  Zudilin/Nesterenko pair (0,n): ", tab[1,n+1]);
}
run(2,400); run(3,400); run(4,400); run(5,300);
\q
