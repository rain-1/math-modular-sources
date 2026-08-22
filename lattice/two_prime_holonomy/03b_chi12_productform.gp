default(parisizemax, 8000000000);
default(realprecision, 60);
read("/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/lib12.gp");
{
print("=== conductor-12 Catalan row, PRODUCT form: identity check & denominators ===");
foreach([[1,1,1],[2,2,1],[2,1,2],[2,2,2],[4,4,2],[3,3,3],[6,6,3],[4,4,4],[8,8,4]], v,
  my(r = c12rowB(v[1],v[2],v[3]), q = r[1], p = r[2], dq = denominator(q), fq);
  fq = if(dq==1, matrix(0,2), factor(dq));
  printf("a1=%d a5=%d b=%d:  Q = %s\n", v[1],v[2],v[3], q);
  printf("            den(Q) primes: %s   den(P) primes: %s\n",
     if(dq==1,"1",fq[,1]~), if(denominator(p)==1,"1",factor(denominator(p))[,1]~));
  printf("            Q*G-P = %.30e\n", q*Catalan - p));
}
quit;
