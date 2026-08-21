/* Exact moments  M(m,j) = int_0^1 int_0^1 K_Z^{(m)}(x,y) w(x,y)^j dx dy,
   w = xy/(1-xy),  K_Z^{(m)} = x^{m-1/2}(1-x)^m y^m (1-y)^{m-1/2}/(1-xy)^{m+1}.
   = Beukers/Nesterenko integral I(m-1/2+j, m, m+j, m-1/2, m+1+j).
   Generalisation of nest() in lattice/catalan_positivity/rows_common.gp (which is
   the case (m,j) = (3n,n)).  Returns [A,B] with M(m,j) = A*G - B, A,B in Q. */
{
nestgen(m,jj) =
 my(N=m, TOP=m+jj,
    cc=m!*prod(i=0,m-2,3/2+i)/TOP!,
    num=cc*prod(k=0,TOP-1,('t-k)), tail=vector(N+2), pre=vector(N+2), M, rhs, sol);
 pre[1]=1; for(j=1,N, pre[j+1]=pre[j]*('t+j));
 tail[N+2]=1; forstep(j=N,0,-1, tail[j+1]=tail[j+2]*('t+j+1/2)^2);
 my(bs=vector(2*(N+1)), lin=vector(N+2));
 lin[1]=1; for(j=0,N, lin[j+2]=lin[j+1]*('t+j+1/2));
 for(j=0,N, bs[2*j+1]=pre[j+1]*lin[j+2]*tail[j+2];
            bs[2*j+2]=pre[j+1]*lin[j+1]*tail[j+2]);
 my(deg=2*N+1);
 M = matrix(deg+1, 2*(N+1), r, c, polcoeff(bs[c], r-1, 't));
 rhs = vectorv(deg+1, r, polcoeff(num, r-1, 't));
 sol = matsolve(M, rhs);
 my(A1=vector(N+1,j,sol[2*(j-1)+1]), A2=vector(N+1,j,sol[2*(j-1)+2]));
 my(B = sum(j=0,N, A2[j+1]),
    C = sum(j=1,N, sum(v=0,j-1, (v!/prod(i=0,v-1,3/2+i))*(A1[j+1]+A2[j+1]/(v+1/2)))));
 [4*B, C];   /* M(m,j) = 4B*G - C */
}
