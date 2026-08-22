default(parisizemax, 8000000000);
default(realprecision, 80);
read("/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/lib12.gp");
\\ Wide-window rate profile of the conductor-12 family.
\\ All rates per unit of b, measured as (X(b2)-X(b1))/(b2-b1) with b1=24, b2=48.
raw(aa, bb) =
{ my(r = c12row(aa,bb), q = r[1], p = r[2], dn = denominator(q),
     v2 = valuation(dn,2), v3 = valuation(dn,3));
  [ log(abs(q*1.0)), log(abs(q*Catalan - p)), v2, v3,
    log((dn/2^v2/3^v3)*1.0), log(lcm(dn,denominator(p))*1.0),
    valuation(numerator(q),2), valuation(numerator(q),3) ];
};
{
print("=== conductor-12 family: wide-window rates, b: 24 -> 48 ===");
print("alpha    logLam    loglam    kap2     kap3     nu       eta");
foreach([[1,2],[2,3],[1,1],[3,2],[2,1],[5,2],[3,1]], fr,
  my(nu0 = fr[1], de0 = fr[2], b1 = 24, b2 = 48, u1, u2, d = b2-b1);
  u1 = raw(nu0*b1/de0, b1); u2 = raw(nu0*b2/de0, b2);
  printf("%d/%d    %8.4f  %8.4f  %7.3f  %7.3f  %7.3f  %7.3f\n", nu0,de0,
     (u2[1]-u1[1])/d, (u2[2]-u1[2])/d, (u2[3]-u1[3])/d, (u2[4]-u1[4])/d,
     (u2[5]-u1[5])/d, (u2[6]-u1[6])/d));
print("");
print("=== integrality structure: is den(Q_b) supported at {2,3}? ===");
foreach([[1,1,24],[1,1,36],[1,1,48],[2,1,24],[2,1,36],[1,2,48]], fr,
  my(r = c12row(fr[1]*fr[3]/fr[2], fr[3]), dq = denominator(r[1]),
     dp = denominator(r[2]), fq = factor(dq));
  printf("alpha=%d/%d b=%d: den(Q) has %d distinct primes, largest = %d, v2=%d v3=%d\n",
     fr[1],fr[2],fr[3], #fq[,1], fq[#fq[,1],1], valuation(dq,2), valuation(dq,3)));
}
quit;
