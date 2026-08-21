/* descent.gp -- hypothesis (b) of Theorem F for the census rows.
   For each modular Apery row (Gamma_0(N), t, F) given as eta quotients:
     * Ligozat cusp divisors of t and F, deg(t) on X_0(N)
     * for every Atkin-Lehner W_Q (Q || N): does W_Q preserve div(t)?
       and (numerically, in the upper half plane) the scalar  t|W_Q / t
       and the weight-w eigenvalue  F|_w W_Q / F.
   Output feeds THEOREM_F_HYPOTHESES.md.                                   */

\\ ---------- Ligozat ----------
\\ ord_{1/c} of prod eta(d tau)^{e_d} on X_0(N), in units of the local
\\ uniformiser (i.e. the order of vanishing as a function on X_0(N)).
ligoz(N, L, c) =
{ my(s = 0);
  for(i=1,#L, my(d=L[i][1], e=L[i][2]); s += gcd(c,d)^2*e/d);
  N/(24*c*gcd(c,N/c)) * s;
}
divisorL(N, L) = { my(v=List()); fordiv(N, c, listput(v,[c, ligoz(N,L,c)])); Vec(v); }
degt(N, L) = { my(s=0); fordiv(N, c, my(o=ligoz(N,L,c)); if(o<0, s -= o)); s; }

\\ ---------- Atkin-Lehner matrices ----------
\\ W_Q = [Q*x, y; N*z, Q*w], det = Q, needs Q*x*w - (N/Q)*y*z = 1.
ALmat(N,Q) =
{ my(M=N/Q, g, u, v);
  if(gcd(Q,M)!=1, error("Q not exact divisor"));
  \\ solve Q*x*w - M*y*z = 1 ; take z=1, y=-1 => Q*x*w + M = 1 impossible.
  \\ use bezout: Q*a - M*b = 1
  my(bz = bezout(Q, M));          \\ bz = [u,v,g] with u*Q+v*M = 1
  u = bz[1]; v = bz[2];
  \\ x=u, w=1, y=-v, z=1 : Q*u*1 - M*(-v)*1 = Q*u + M*v = 1 . OK
  [Q*u, -v; N, Q];
}
actAL(W, tau) = (W[1,1]*tau + W[1,2])/(W[2,1]*tau + W[2,2]);

\\ ---------- cusp permutation of W_Q ----------
\\ denominator c |-> c' with gcd(c',Q) = Q/gcd(c,Q), gcd(c',N/Q)=gcd(c,N/Q)
cuspAL(N,Q,c) = { my(M=N/Q, a=Q/gcd(c,Q), b=gcd(c,M)); a*b; }

\\ ---------- numeric eta quotient ----------
etaN(L, tau) = { my(r=1.0); for(i=1,#L, r *= eta(L[i][1]*tau, 1)^L[i][2]); r; }
\\ full eta quotient value (eta() with flag 1 gives the true eta(tau))

report(name, N, Lt, LF, wt, p) =
{ my(dt = divisorL(N,Lt), dF = divisorL(N,LF), D = degt(N,Lt), tau, W, sc, ef);
  print("=== ", name, "  N=", N, "  wt F=", wt, "  slope prime p=", p);
  print("   div(t): ", dt, "    deg t = ", D);
  print("   div(F): ", dF);
  tau = 0.31 + 0.77*I;
  fordiv(N, Q, if(Q>1 && gcd(Q,N/Q)==1,
    my(okdiv=1, permok=1);
    \\ check divisor preserved
    fordiv(N, c, my(cc=cuspAL(N,Q,c));
      if(ligoz(N,Lt,c) != ligoz(N,Lt,cc), okdiv=0));
    W = ALmat(N,Q);
    sc = etaN(Lt, actAL(W,tau)) / etaN(Lt, tau);
    ef = etaN(LF, actAL(W,tau)) / (etaN(LF,tau) * (W[2,1]*tau+W[2,2])^wt / Q^(wt/2));
    print("   W_", Q, ": div-preserved=", okdiv,
          "   t|W/t = ", sc,
          "   F|_wW/F = ", ef,
          if(p>0, concat(concat("   [p|Q: ", Str(Mod(Q,p)==0)), "]"), ""));
  ));
  print();
}
