\p 300
default(parisize,8000000000);
G = Catalan;
/* R_n(t) = c * prod_{k=0}^{4n-1}(t-k) / prod_{r=0}^{3n}(t+r+1/2)^2 */
{
nest(n) =
 my(N=3*n, c=(3*n)!*prod(i=0,3*n-2,3/2+i)/(4*n)!,
    num=c*prod(k=0,4*n-1,('t-k)),
    tail=vector(N+2), pre=vector(N+2), M, rhs, sol, A1, A2, B, C, v);
 /* pre[j+1] = prod_{i=1}^{j}(t+i) ; tail[j+1] = prod_{r=j}^{N}(t+r+1/2)^2 */
 pre[1]=1; for(j=1,N, pre[j+1]=pre[j]*('t+j));
 tail[N+2]=1; forstep(j=N,0,-1, tail[j+1]=tail[j+2]*('t+j+1/2)^2);
 /* basis polys */
 my(bs=vector(2*(N+1)), lin=vector(N+2));
 lin[1]=1; for(j=0,N, lin[j+2]=lin[j+1]*('t+j+1/2));  /* lin[k+1]=prod_{r=0}^{k-1}(t+r+1/2) */
 for(j=0,N,
   bs[2*j+1] = pre[j+1]*lin[j+2]*tail[j+2];   /* coeff of A1_j : prod_{r=0}^{j} */
   bs[2*j+2] = pre[j+1]*lin[j+1]*tail[j+2]);  /* coeff of A2_j : prod_{r=0}^{j-1} */
 my(deg=6*n+1, K=2*(N+1));
 M = matrix(deg+1, K, r, cc, polcoeff(bs[cc], r-1, 't));
 rhs = vectorv(deg+1, r, polcoeff(num, r-1, 't));
 sol = matsolve(M, rhs);
 A1 = vector(N+1, j, sol[2*(j-1)+1]);
 A2 = vector(N+1, j, sol[2*(j-1)+2]);
 [A1,A2];
}
{
for(n=1,4,
  my(N=3*n, r=nest(n), A1=r[1], A2=r[2], bad=0);
  for(j=0,N, my(f=2^(-14*n+2*j+1)*(8*n+2*j)!*j!*(6*n)!/((4*n)!*(4*n+j)!*(3*n-j)!^2*(2*j)!^2));
    if(A2[j+1]!=f, bad=1; print("A2 mismatch n=",n," j=",j," got ",A2[j+1]," formula ",f)));
  print("n=",n," A2 formula (2.11) matches partial fractions: ", if(bad,"NO","YES")));
}
\q
