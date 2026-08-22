default(parisizemax, 8000000000);
default(realprecision, 80);
read("/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/lib12.gp");
raw(aa, bb) =
{ my(r = c12rowB(aa,aa,bb), q = r[1], p = r[2], dn = denominator(q),
     v2 = valuation(dn,2), v3 = valuation(dn,3));
  [ log(abs(q*1.0)), log(abs(q*Catalan - p)), v2, v3,
    log((dn/2^v2/3^v3)*1.0), log(lcm(dn,denominator(p))*1.0) ]; };
{
print("=== conductor-12 product form, small alpha; wide window b: 32 -> 64 ===");
print("alpha   logLam    loglam    kap2    kap3     nu       eta     rho_int = logLam+nu+sum kap_p log p");
foreach([[1,4],[1,2],[3,4]], fr,
  my(a1 = fr[1]*32/fr[2], a2 = fr[1]*64/fr[2], u1 = raw(a1,32), u2 = raw(a2,64), d = 32, v);
  v = vector(6, i, (u2[i]-u1[i])/d);
  printf("%d/%d   %8.4f  %8.4f  %6.3f  %6.3f  %7.3f  %7.3f   %8.4f\n", fr[1],fr[2],
    v[1],v[2],v[3],v[4],v[5],v[6], v[1]+v[5]+v[3]*log(2)+v[4]*log(3)));
}
quit;
