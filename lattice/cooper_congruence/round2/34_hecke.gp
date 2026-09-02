\\ 34_hecke.gp -- Task (3): the exact T_{p^2} eigen-identity for PZ's level-one inputs.
\\ CLAIM (derived from PZ Lemma 2's uniqueness argument, dim S_4(SL_2(Z))=0):
\\    f | T_{p^2}  =  psi(p) * p * f  +  p^3 * c * g_{m0 p^2},
\\ where f = c q^{-m0} + O(q) in S^{!,+}_{5/2}, psi = chi_{D0}, D0 = -m0 the pole
\\ discriminant, and g_M = q^{-M} + O(q) the integral basis element.
default(parisize, 6000000000);
read("33_coeffs.gp");         \\ defines A4A, A4B : vectors of a(n), n = -4 .. NMAXC
{ acof(V,n) = if(n < -4 || n > NMAXC, 0, V[n+5]); }
\\ f|T_{p^2} at index n, k=2 :  a(p^2 n) + p*(n/p)*a(n) + p^3*a(n/p^2)
{ Tp2(V,p,n) = my(s);
  s = acof(V,p^2*n) + p*kronecker(n,p)*acof(V,n);
  if(n%(p^2)==0, s += p^3*acof(V,n/p^2));
  s; }
{ runtest(nam, V, m0, D0, cden) =
  my(NT);
  print("");
  print("=== ", nam, " :  f|T_\{p^2\} = chi_", D0, "(p)*p*f + p^3*(1/", cden, ")*g_\{", m0, "p^2\} ===");
  forprime(p=2, 17,
    my(lam = kronecker(D0,p)*p, ok=1, bad=0, first=0, prin, NT2, lead);
    NT2 = NMAXC \ p^2;
    \\ principal part of E = f|T_{p^2} - lam*f
    prin = [];
    for(n=-m0*p^2, -1, my(e = Tp2(V,p,n) - lam*acof(V,n)); if(e!=0, prin = concat(prin, [[n,e]])));
    lead = Tp2(V,p,-m0*p^2) - lam*acof(V,-m0*p^2);
    \\ integrality of cden*E/p^3 on 0 <= n <= NT2
    for(n=0, NT2, my(e = Tp2(V,p,n) - lam*acof(V,n));
      if(denominator(cden*e/p^3)!=1, bad++; if(first==0,first=n)));
    print("  p=", p, "  lambda=", lam,
          "   principal part of E: ", prin,
          "   E(0)=", Tp2(V,p,0)-lam*acof(V,0),
          "   ", cden, "*E/p^3 integral for 0<=n<=", NT2, ": ", if(bad==0,"YES",Str("NO (",bad," fails, first n=",first,")")),
          "   leading coeff*", cden, "/p^3 = ", cden*lead/p^3));
}
runtest("f4a  (pole rho, disc -3)", A4A, 3, -3, 64);
runtest("f4b  (pole i,   disc -4)", A4B, 4, -4, 108);
print("");
print("=== control: is any OTHER lambda possible?  residual principal part for wrong lambda ===");
{ forprime(p=5,13,
    my(lamgood=kronecker(-3,p)*p);
    print("  f4a p=",p," : E(-3) for lambda=0,", lamgood, ", ", -lamgood, " : ",
      [Tp2(A4A,p,-3), Tp2(A4A,p,-3)-lamgood*acof(A4A,-3), Tp2(A4A,p,-3)+lamgood*acof(A4A,-3)]));
}
quit;
