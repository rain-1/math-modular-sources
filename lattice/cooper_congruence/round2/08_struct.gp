\\ 08_struct.gp -- the structural reformulations of Conjecture 4.1, verified exactly.
read("lib.gp");
M = 400;
{
for(k=1,3,
  my(cp, b, g, ok1, ok2, ok3, Phc, flat, prod, e, R, S, T);
  cp = CPvec(k,M);
  b  = Bvec(k,cp);
  g  = vector(M, n, b[n]/n^2);
  \\ (i) gamma integral
  ok1 = 1; for(n=1,M, if(denominator(g[n])!=1, ok1=0));
  \\ (ii) Xi = sum_e psi(e) (D^2 K)(q^e):  c'(n) = sum_{d|n} psi(n/d) d^2 gamma(d)
  ok2 = 1;
  for(n=1,M, my(s=0); fordiv(n,d, s += psin(k,n/d)*d^2*g[d]); if(s!=cp[n], ok2=0));
  \\ (iii) Phi_flat := sum_e mu(e)psi(e) e Phi(q^e) has coefficients n*beta(n) = a(n^2)
  Phc = Cvec(k,M);
  flat = vector(M);
  for(n=1,M, my(s=0); fordiv(n,e, s += moebius(e)*psin(k,e)*e*Phc[n/e]); flat[n]=s);
  ok3 = 1; for(n=1,M, if(flat[n] != n*b[n], ok3=0));
  print("row ",NAM[k],":  gamma integral (n<=",M,"): ",ok1,
        "   Xi = sum_e psi(e)(D^2K)(q^e): ",ok2,
        "   [q^n] sum_e mu(e)psi(e) e Phi(q^e) = n beta(n): ",ok3);
);
}
\\ the q-product g = prod (1-q^n)^{-n gamma(n)}, psi=1 rows
{
for(k=1,3,
  my(cp,b,g,P,MM=120,ok);
  cp = CPvec(k,MM); b = Bvec(k,cp); g = vector(MM,n,b[n]/n^2);
  P = 1 + O('q^MM);
  for(n=1,MM-1, P = P*(1-'q^n+O('q^MM))^(-n*g[n]));
  ok = 1; for(n=0,MM-1, if(denominator(polcoeff(P,n))!=1, ok=0));
  print("row ",NAM[k],"  prod (1-q^n)^{-n gamma(n)} in Z[[q]] (n<",MM,"): ",ok,
        "   first coeffs ", vector(9,i,polcoeff(P,i-1)));
);
}
quit;
