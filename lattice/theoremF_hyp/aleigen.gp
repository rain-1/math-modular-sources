/* aleigen.gp -- the Atkin-Lehner oldform eigenvalue check for hypothesis (b).
   Setting: Phi = sum_{d | N/M} c_d V_d E, E primitive of level M, weight k1=w+2,
   nebentypus chi = psi*phi;  Theta = D^{-(w+1)} Phi = sum_d c_d d^{-(w+1)} V_d cE.
   For Q || N with gcd(Q,M)=1 the Atkin-Lehner involution acts on the oldform
   basis at weight k by
        V_d  |-->  lambda * chi(Q/d1) * (Q/d1^2)^{k/2} * V_{(Q/d1) d2},
   d = d1*d2 with d1 the Q-part of d, lambda = chi(w) a global scalar.
   We set lambda = 1 (it cancels between F and Theta, whose nebentypi are
   inverse).  eigW(cs, ds, N, M, Q, k, chi) returns the eigenvalue, or 0 if the
   vector is not an eigenvector.                                             */

qpart(d,Q) = { my(r=1); forprime(l=2, Q, if(Q%l==0, r *= l^valuation(d,l))); r; }

eigW(cs, ds, Q, k, chi) =
{ my(n=#ds, img=vector(n), pos, ev, ok=1);
  for(i=1,n, my(d=ds[i], d1=qpart(d,Q), d2=d/d1, e=(Q/d1)*d2, f=Q/d1^2, sc);
    sc = chi(Q/d1) * f^(k/2);
    pos = 0; for(j=1,n, if(ds[j]==e, pos=j));
    if(pos==0, error("basis not closed: ", e));
    img[pos] += cs[i]*sc);
  ev = 0;
  for(i=1,n, if(cs[i]!=0, ev = img[i]/cs[i]; break));
  for(i=1,n, if(img[i] != ev*cs[i], ok=0));
  if(ok, ev, "NOT-EIGENVECTOR");
}

/* full report for one row */
alrow(name, w, ds, cs, Q, chi) =
{ my(k1=w+2, k2=-w, cs2 = vector(#ds, i, cs[i]/ds[i]^(w+1)));
  print("  ", name, "  W_", Q, ":  eig(Phi, wt ", k1, ") = ", eigW(cs, ds, Q, k1, chi),
        "    eig(Theta, wt ", k2, ") = ", eigW(cs2, ds, Q, k2, chi));
}
