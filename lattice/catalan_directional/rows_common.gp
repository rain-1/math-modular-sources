/* Exact Zudilin (index 3n) and Nesterenko (4,7) integer rows, in the
   normalisation of NearCriticalAssembly.lean / the literal-window paper:
     X_n = 2^{e_{3n}} D_{6n}^2 Q_{3n},  Y_n = 2^{e_{3n}} D_{6n}^2 P_{3n}
     V_n = 4^{7n+1} D_{6n}^2 B_n,       U_n = 4^{7n} D_{6n}^2 C_n
   (adapted from lattice/catalan_audit/{rows.gp,nest2.gp}). */
default(parisize,8000000000);
ee(m) = min(6*m, 4*m+3+logint(2*m-1,2));
{
zud(M) = my(Q=vector(M+1), P=vector(M+1), a,b,c);
  Q[1]=1; Q[2]=7/4; P[1]=0; P[2]=13/8;
  for(m=1,M-1,
    a=(2*m+1)^2*(2*m+2)^2*(20*m^2-8*m+1);
    b=3520*m^6+5632*m^5+2064*m^4-384*m^3-156*m^2+16*m+7;
    c=(2*m-1)^2*(2*m)^2*(20*m^2+32*m+13);
    Q[m+2]=(b*Q[m+1]+c*Q[m])/a;
    P[m+2]=(b*P[m+1]+c*P[m])/a);
  [Q,P];
}
{
nest(n) =
 my(N=3*n, c=(3*n)!*prod(i=0,3*n-2,3/2+i)/(4*n)!,
    num=c*prod(k=0,4*n-1,('t-k)), tail=vector(N+2), pre=vector(N+2), M, rhs, sol);
 pre[1]=1; for(j=1,N, pre[j+1]=pre[j]*('t+j));
 tail[N+2]=1; forstep(j=N,0,-1, tail[j+1]=tail[j+2]*('t+j+1/2)^2);
 my(bs=vector(2*(N+1)), lin=vector(N+2));
 lin[1]=1; for(j=0,N, lin[j+2]=lin[j+1]*('t+j+1/2));
 for(j=0,N, bs[2*j+1]=pre[j+1]*lin[j+2]*tail[j+2];
            bs[2*j+2]=pre[j+1]*lin[j+1]*tail[j+2]);
 my(deg=6*n+1, K=2*(N+1));
 M = matrix(deg+1, K, r, cc, polcoeff(bs[cc], r-1, 't));
 rhs = vectorv(deg+1, r, polcoeff(num, r-1, 't));
 sol = matsolve(M, rhs);
 [vector(N+1,j,sol[2*(j-1)+1]), vector(N+1,j,sol[2*(j-1)+2])];
}
{
nestrow(n) = my(N=3*n, r=nest(n), A1=r[1], A2=r[2], B, C, D=lcm(vector(6*n,i,i)));
 B = sum(j=0,N, A2[j+1]);
 C = sum(j=1,N, sum(v=0,j-1, (v!/prod(i=0,v-1,3/2+i))*(A1[j+1]+A2[j+1]/(v+1/2))));
 [4^(7*n+1)*D^2*B, 4^(7*n)*D^2*C];
}
{
zudrow(n) = my(m=3*n, D=lcm(vector(6*n,i,i)), ZZ);
 ZZ=zud(m); [2^ee(m)*D^2*ZZ[1][m+1], 2^ee(m)*D^2*ZZ[2][m+1]];
}
