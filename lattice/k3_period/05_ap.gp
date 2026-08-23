/* 05_ap.gp -- transcendental Frobenius traces of the elliptic K3.
 *
 * For an elliptic surface pi: X -> P^1 with section, F = R^1 pi_* Q_l,
 *   sum_{t in P^1(F_p)} a_t = - tr(Frob | H^1(P^1, j_* F)),
 * (the local invariants at the bad fibres cancel between H^1_c and the
 * skyscrapers), and H^1(P^1,j_*F) is exactly the orthogonal complement of
 * the trivial lattice in H^2(X), i.e. T (x) Q_l when rho = 20.  Hence
 *   a_p := tr(Frob|T) = - sum_{t in P^1(F_p)} a_t
 *        = sum_{t in F_p} sum_{x in F_p} chi(x^3 - 27 c4(t) x - 54 c6(t)),
 * the fibre at t = infinity being of type IV* (additive, a = 0).
 * Here a_t = p - #affine points is the usual trace with a_t = +1/-1/0 at
 * split/nonsplit multiplicative / additive fibres.
 */
default(parisizemax, 6000000000);

c4f(t) = 864*t^5 + 40*t^4 - 256*t^3 + 96*t^2 - 16*t + 1;
c6f(t) = -(5832*t^8 + 34560*t^7 - 30016*t^6 + 12624*t^5 - 4380*t^4 + 1280*t^3 - 240*t^2 + 24*t - 1);

apk3(p) =
{ my(ch, s, aa, bb, v);
  ch = vector(p, i, kronecker(i-1, p));   /* ch[k+1] = chi(k) */
  s = 0;
  for(tt = 0, p-1,
    aa = (-27*c4f(tt)) % p; bb = (-54*c6f(tt)) % p;
    v = 0;
    for(xx = 0, p-1, v += ch[((xx^3 + aa*xx + bb) % p) + 1]);
    s += v);
  s;
}

forprime(p = 5, 200, print(p, " ", apk3(p)));
quit;
