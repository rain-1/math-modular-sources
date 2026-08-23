/* Denominator type of the moving-EMN period H in the coordinate x = z/2:
   [x^{n+1}] H = T_n/(n+1), T_n = sum_k C(n,k)/(2k+1).  Claim: type n*[1..2n] (e=1,b=2). */
{
print("n | den(T_n/(n+1)) | is it | (n+1)*lcm(1..2n+1) ? | minimal (e,b) test");
for(n=0,60,
  my(tn = sum(k=0,n, binomial(n,k)/(2*k+1)), c = tn/(n+1), dn = denominator(c),
     L1 = lcm(vector(2*n+1,k,k)), L2 = lcm(vector(max(n,1),k,k)));
  if(n<=15 || n%10==0,
    print("  n=",n,"  den=",dn,
      "   | (n+1)*lcm(1..2n+1)? ",((n+1)*L1) % dn == 0,
      "   | (n+1)*lcm(1..n)? ",((n+1)*L2) % dn == 0,
      "   | lcm(1..2n+1)? ", L1 % dn == 0)));
print();
print("growth: (1/n)*log den");
for(n=10,60,if(n%10==0,
  my(tn=sum(k=0,n,binomial(n,k)/(2*k+1)), dn=denominator(tn/(n+1)));
  print("  n=",n,"  (1/n)log den = ", log(dn*1.)/n)));
}
quit();
