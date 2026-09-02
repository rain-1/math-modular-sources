\\ 37_thm.gp -- numerical verification of the LEVEL-ONE THEOREM proved in FINDINGS_PZ.md:
\\   S_a(m) := a(p^{2a} m^2) + (1-psi(p)) * sum_{j=0}^{a-1} p^{a-j} a(p^{2j} m^2)
\\           = p^{3a} * c * g_{m0 p^{2a}}(m^2),
\\ in particular v_p(S_a(m)) >= 3a  for every prime p (p-integral input), a>=0, p nmid m.
\\ Uses the half-integral coefficients built in 33_coeffs.gp.
default(parisize, 8000000000);
read("33_coeffs.gp");
{ acof(V,n) = if(n < -4 || n > NMAXC, 0, V[n+5]); }
{ Sa(V, p, a, m, psi) = my(s = acof(V, p^(2*a)*m^2));
  for(j=0, a-1, s += (1-psi)*p^(a-j)*acof(V, p^(2*j)*m^2));
  s; }
{ run(nam, V, D0, cden) =
  my(tot=0, bad=0, mn=99, lst=[]);
  print("");
  print("=== ", nam, " : v_p(S_a(m)) >= 3a ?   (a>=1, p^{2a} m^2 <= ", NMAXC, ") ===");
  forprime(p=2, 53,
    my(psi = kronecker(D0,p), b=0, c=0, mnp=99, amax=0);
    for(a=1, 6,
      if(p^(2*a) > NMAXC, break);
      amax = a;
      for(m=1, sqrtint(NMAXC\p^(2*a)),
        if(m%p==0, next);
        my(s = Sa(V,p,a,m,psi), v);
        c++;
        v = if(s==0, 10^6, valuation(s,p));
        if(v < 3*a, b++);
        if(v - 3*a < mnp, mnp = v-3*a)));
    if(c>0, print("   p=",p," psi=",psi," a<=",amax," tested ",c," pairs : ",
      if(b==0, Str("PASS  (min excess ", mnp, ")"), Str("FAIL ", b, "/", c))));
    if(b>0, bad++));
  print("   -> primes with a failure: ", bad, if(bad>0, "  (expected: only the prime dividing the denominator of c)", ""));
}
run("f4a (m0=3, psi=chi_-3, c=1/64)",  A4A, -3,  64);
run("f4b (m0=4, psi=chi_-4, c=-1/108)", A4B, -4, 108);
print("");
print("=== the exact identity S_a(m) = p^{3a} c g_{m0 p^{2a}}(m^2) : print S_a(m)/p^{3a} ===");
{ forprime(p=3,11,
   for(a=1,2, if(p^(2*a)>NMAXC, break);
     print("  f4a p=",p," a=",a," : 64*S_a(1)/p^(3a) = ", 64*Sa(A4A,p,a,1,kronecker(-3,p))/p^(3*a),
           "   (must be the integer g_{3 p^{2a}}(1))")));
}
quit;
