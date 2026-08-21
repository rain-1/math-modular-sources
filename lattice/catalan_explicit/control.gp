/* Controls.
   (a) the honest 2-row bridge lattice objective (1/n)log|q_nG-p_n| at the same small n;
   (b) the FULL moment lattice in dimension m+1 with mixed signs (positivity given up),
       modulus g = 2^t, LLL in the metric |M_j|.  Shows how much of the gap is due to
       the positivity restriction rather than to the moment family being weak. */
\p 400
G=Catalan;
{ klat2(X,Y,V,U,M) = my(K=matkerint([X,V,M,0,0,0; Y,U,0,M,0,0]));
  mathnf(matrix(2,#K[1,],i,j,K[i,j])); }
{
tworow(n,k)= my(D=lcm(vector(6*n,i,i)), M=D^2*2^floor(k*n),
   zr=zudrow(n), nr=nestrow(n), B0=klat2(zr[1],zr[2],nr[1],nr[2],M),
   LZ=(zr[1]*G-zr[2])/M, LN=(nr[1]*G-nr[2])/M, best=0);
 for(i=-300,300, for(j=-300,300, my(cz=i*B0[1,1]+j*B0[1,2], cn=i*B0[2,1]+j*B0[2,2]);
   if(cz||cn, my(v=abs(cz*LZ+cn*LN), pos=(sign(cz*LZ)>=0 && sign(cn*LN)>=0));
     if(pos && (best==0||v<best), best=v))));
 log(best)/n;
}
{
fullmom(n,t)= my(m=3*n, A=vector(m+1),B=vector(m+1),Mv=vector(m+1),d,a,b,g=2^t);
 for(j=0,m, my(r=nestgen(m,j)); A[j+1]=r[1];B[j+1]=r[2];Mv[j+1]=r[1]*G-r[2]);
 d=lcm(vector(m+1,j,lcm(denominator(A[j]),denominator(B[j]))));
 a=vector(m+1,j,d*A[j]); b=vector(m+1,j,d*B[j]);
 my(Mm=matrix(2,m+3), K);
 for(j=1,m+1, Mm[1,j]=a[j]; Mm[2,j]=b[j]); Mm[1,m+2]=g; Mm[2,m+3]=g;
 K=matkerint(Mm); K=matrix(m+1,#K[1,],i,j,K[i,j]);
 my(sc=10^350, W=matrix(m+1,#K[1,],i,j,round(K[i,j]*abs(Mv[i])*sc/vecmin(vector(m+1,l,abs(Mv[l]))))));
 my(T=qflll(W), R=K*T, best=0);
 for(c=1,#R[1,], my(v=sum(i=1,m+1,R[i,c]*Mv[i]));
   if(v!=0, my(q=sum(i=1,m+1,R[i,c]*a[i])/g); if(best==0||abs(v)*d/g<best, best=abs(v)*d/g)));
 [log(best)/n, log(1.0*d)/n];
}
for(n=2,5, printf("n=%d  two-row cone-min (k=22.4): %8.4f\n", n, tworow(n,22.4)));
for(n=2,5, for(t=0,3, my(r=fullmom(n,floor(t*8*n))); printf("n=%d t=%d(g=2^%d) full-moment-lattice (1/n)log|qG-p| = %8.4f   (1/n)log d=%7.3f\n", n,t,floor(t*8*n),r[1],r[2])));
\q
