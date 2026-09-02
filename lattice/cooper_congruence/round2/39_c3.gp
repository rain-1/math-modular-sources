\\ 39_c3.gp -- corrected exponent in (C3):  A(p^r n) = lam^r A(n)  mod p^{r+2}  (k=2)
\\ and mod p^{2r+3} (k=3).  This is exactly Cooper's eq:magnetic (i) with lam = psi(p)p^{k-1}.
default(parisize, 6000000000);
NN = 1300; MM = 1250;
{ E4s = 1 + 240*sum(n=1, NN-1, sigma(n,3)*'q^n) + O('q^NN); }
{ E6s = 1 - 504*sum(n=1, NN-1, sigma(n,5)*'q^n) + O('q^NN); }
DEL = (E4s^3 - E6s^2)/1728;
A4a = vector(MM,m,polcoeff(DEL/E4s^2,m));
A4b = vector(MM,m,polcoeff(E4s*DEL/E6s^2,m));
A6  = vector(MM,m,polcoeff(E6s*DEL/E4s^3,m));
{ Ac(V,n) = if(n<1 || n>#V, 0, V[n]); }
{ tst(nam, V, D0, kk) =
  print("");
  print("=== ", nam, " : A(p^r n) = lam^r A(n) mod p^{r+2} (k=2) / p^{2r+3} (k=3) ===");
  forprime(p=2,19,
    my(lam=kronecker(D0,p)*p^(kk-1), res="");
    for(r=1,5, my(b=0, ex=if(kk==2, r+2, 2*r+3));
      for(n=1, MM\p^r, if((Ac(V,p^r*n) - lam^r*Ac(V,n)) % p^ex != 0, b++));
      res = Str(res, " r=",r,":", if(b==0,"PASS",Str("FAIL",b))));
    print("  p=",p," lam=",lam, res));
}
tst("F4a (psi=chi_-3)", A4a, -3, 2);
tst("F4b (psi=chi_-4)", A4b, -4, 2);
tst("F6  (lam=p^2)",    A6,   1, 3);
quit;
