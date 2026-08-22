default(parisizemax, 8000000000);
{
nestpf(n) =
 my(NB=3*n, cc=(3*n)!*prod(i=0,3*n-2,3/2+i)/(4*n)!,
    num=cc*prod(k=0,4*n-1,('t-k)), tail=vector(NB+2), pre=vector(NB+2), mat, rhs, sol);
 pre[1]=1; for(j=1,NB, pre[j+1]=pre[j]*('t+j));
 tail[NB+2]=1; forstep(j=NB,0,-1, tail[j+1]=tail[j+2]*('t+j+1/2)^2);
 my(bs=vector(2*(NB+1)), lin=vector(NB+2));
 lin[1]=1; for(j=0,NB, lin[j+2]=lin[j+1]*('t+j+1/2));
 for(j=0,NB, bs[2*j+1]=pre[j+1]*lin[j+2]*tail[j+2];
             bs[2*j+2]=pre[j+1]*lin[j+1]*tail[j+2]);
 my(deg=6*n+1, KK=2*(NB+1));
 mat = matrix(deg+1, KK, r, c2, polcoeff(bs[c2], r-1, 't));
 rhs = vectorv(deg+1, r, polcoeff(num, r-1, 't));
 sol = matsolve(mat, rhs);
 [vector(NB+1,j,sol[2*(j-1)+1]), vector(NB+1,j,sol[2*(j-1)+2])];
}
{
for(k=1,4, my(n=[5,10,15,20][k], t0=getabstime(), r=nestpf(n));
  print("n=",n,"  time_ms=", getabstime()-t0));
}
quit;
