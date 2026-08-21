\p 60
{
for(t=1,6, my(n=[4,8,12,16,20,24][t], D=lcm(vector(6*n,i,i)), S=D^2,
  zr=zudrow(n), nr=nestrow(n), X=zr[1],Y=zr[2],V=nr[1],U=nr[2], aZ=X/S,aN=V/S);
 my(h=aZ*U-aN*Y);
 printf("n=%2d logS/n=%.4f v2(X)=%d v2(V)=%d v2(aZ)=%d v2(aN)=%d v2(Y)=%d v2(U)=%d v2(h)/n=%.3f  gcd(Y,U,S)-> v2=%d\n",
  n, log(S)/n, valuation(X,2),valuation(V,2),valuation(aZ,2),valuation(aN,2),valuation(Y,2),valuation(U,2), valuation(h,2)/n, valuation(gcd(gcd(Y,U),S),2));
 my(T=2^floor(22.4*n), M=S*T);
 my(K=matkerint([Y,U,S,0; aZ,aN,0,T]), B=mathnf(matrix(2,#K[1,],i,j,K[i,j])));
 printf("     log idx/n=%.4f   log M/n=%.4f\n", log(abs(matdet(B)))/n, log(M)/n);
);
}
\q
