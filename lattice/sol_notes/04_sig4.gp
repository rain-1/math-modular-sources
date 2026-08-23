/* SOL note 2, sections 2-3: signature-four realisation.
   F(t) = 2F1(1/4,3/4;1;t);  claims:
   (S2) 8t^2(1-t)^2 Y''' + 24t(1-t)(1-2t) Y'' + 2(27t^2-27t+4) Y' + 3(2t-1) Y = 0 for Y=F^2
   (S5) F(64x) = sum binomial(4n,2n)binomial(2n,n) x^n
   (S6) [x^n] D(64x) = binomial(4n,2n)binomial(2n,n)*4*sum_{k<n}(1/(4k+1)-1/(4k+3)),
        and the true denominator type of that coefficient. */
NS = 40;
fser = sum(k=0,NS, (1/4+O(1))*0);  \\ placeholder
Fc(n) = binomial(4*n,2*n)*binomial(2*n,n)/64^n;
{
print("== S5. [t^n] 2F1(1/4,3/4;1;t) = (1/4)_n(3/4)_n/(n!)^2  vs  C(4n,2n)C(2n,n)/64^n");
for(n=0,12,
  my(poch = prod(j=0,n-1,(1/4+j)*(3/4+j))/ (n!)^2);
  print("  n=",n,"  pochhammer=",poch,"  binomial form=",Fc(n),"  equal=",poch==Fc(n)));
}
Ft = sum(n=0,NS, Fc(n)*t^n) + O(t^(NS+1));
Y = Ft^2;
D1 = deriv(Y,t); D2 = deriv(D1,t); D3 = deriv(D2,t);
L = 8*t^2*(1-t)^2*D3 + 24*t*(1-t)*(1-2*t)*D2 + 2*(27*t^2-27*t+4)*D1 + 3*(2*t-1)*Y;
{
print();
print("== S2. symmetric-square ODE applied to Y=F^2:");
print("  L(Y) = ", L);
print("  (zero up to the series truncation order means VERIFIED)");
}
{
print();
print("== S6. coefficients of D(64x) and their exact denominators");
print("  n | c_n = C(4n,2n)C(2n,n)*4*sum(1/(4k+1)-1/(4k+3)) | den(c_n) | d_n=lcm(1..n) | d_n*c_n integral? | min k with den | lcm(1..k)");
for(n=1,40,
  my(bn = binomial(4*n,2*n)*binomial(2*n,n),
     sm = 4*sum(k=0,n-1, 1/(4*k+1) - 1/(4*k+3)),
     cn = bn*sm, dn = denominator(cn), dl = lcm(vector(n,k,k)), kk = 0);
  for(k=1,4*n+4, if(lcm(vector(k,j,j)) % dn == 0, kk=k; break));
  if(n<=12 || n%5==0,
    print("  n=",n," den=",dn,"  d_n=",dl,"  d_n*c_n in Z? ",(dl % dn)==0,
          "   min k: lcm(1..k)^1 kills den = ",kk, "   k/n=",kk*1./n)));
}
quit();
