\p 60
/* Absolute ceiling: total odd content of gcd(V_n,U_n) not already in gcd(X_n,Y_n),
   INCLUDING primes > 6n (so the restriction l<=6n in the main scan is justified). */
{
for(idx=1,8, my(n=[5,8,12,16,20,26,32,40][idx]);
 my(zr=zudrow(n), nr=nestrow(n), gZ=gcd(zr[1],zr[2]), gN=gcd(nr[1],nr[2]));
 gZ>>=valuation(gZ,2); gN>>=valuation(gN,2);
 my(rest=gN, small=1);
 forprime(l=3,6*n, my(v=valuation(rest,l)); if(v, small*=l^v; rest/=l^v));
 printf("n=%2d  log(oddgcd V,U)/n=%.5f  log(oddgcd X,Y)/n=%.5f  ceiling=%.5f  |  part with l>6n: %.5f/n\n",
  n, log(gN)/n, log(gZ)/n, (log(gN)-log(gZ))/n, log(rest)/n));
}
