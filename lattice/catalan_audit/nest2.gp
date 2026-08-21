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
nestrow(n) =
 my(N=3*n, r=nest(n), A1=r[1], A2=r[2], B, C, D=lcm(vector(6*n,i,i)), V, U, J);
 B = sum(j=0,N, A2[j+1]);
 C = sum(j=1,N, sum(v=0,j-1, (v!/prod(i=0,v-1,3/2+i))*(A1[j+1]+A2[j+1]/(v+1/2))));
 V = 4^(7*n+1)*D^2*B; U = 4^(7*n)*D^2*C; J = 4*B*G - C;
 [B,C,V,U,J,D];
}
print("=== Nesterenko (4,7) row ===");
Tp=(3303+437*sqrt(57))/144; Tm=(3303-437*sqrt(57))/144;
printf("2logT=%.8f  2logT^-=%.8f\n", 2*log(Tp), 2*log(Tm));
write("nest_exact.txt","");
{
for(n=1,10,
 my(r=nestrow(n), B=r[1], C=r[2], V=r[3], U=r[4], J=r[5]);
 printf("n=%2d intV=%d intU=%d  J_n=%.6e  U/(V G)=%.16f  log B/n=%.6f log J/n=%.6f  log|VG-U|/n=%.6f\n",
   n, denominator(V)==1, denominator(U)==1, J, U/(V*G), log(B)/n, log(J)/n, log(abs(V*G-U))/n);
 write("nest_exact.txt", n, " ", V, " ", U));
}
A2r=14*log(2)+12+2*log(Tp); E2r=14*log(2)+12+2*log(Tm);
printf("paper A_2=%.10f E_2=%.10f\n",A2r,E2r);
\q
