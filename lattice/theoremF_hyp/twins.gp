/* twins.gp -- the low-level "q -> -q" twin presentations of rows B, eta, delta
   (and the checks that they reproduce the same rows).
   The substitution q -> -q is tau -> tau + 1/2, i.e. conjugation by
   gamma = [1,1/2;0,1]; it is p-adically an automorphism whenever p is odd,
   so it transports hypotheses (b),(c) at p = 3, 5.                        */
PRECSET = 160;
read("lattice/rigidity/lib.gp");

flip(f) = subst(f, q, -q);

/* expand F in powers of t : returns [a_0..a_n] */
expandIn(t, F, n) = { my(r = serreverse(t)); Vec(Pol(subst(F, q, r) + O(q^(n+1))))[1..n+1]; }

rowcheck(name, t, F, arec, n) =
{ my(a = expandIn(t, F, n), ok=1);
  print(name, ":  t = ", Vec(Pol(t+O(q^7))), "   F = ", Vec(Pol(F+O(q^7))));
  print("   a_n from (t,F), n=0..", n, ": ", Vecrev(a)[1..min(n+1,12)]);
  print("   recurrence a_n      , n=0..", n, ": ", arec[1..min(n+1,12)]);
  for(i=1,n+1, if(Vecrev(a)[i] != arec[i], ok=0));
  print("   MATCH: ", ok);
  ok;
}
