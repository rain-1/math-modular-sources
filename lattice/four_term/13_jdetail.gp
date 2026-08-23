/* 13_jdetail.gp -- for a row that passes the J-map test, recover the rational
 * function Jcal = U/V of degree deg J and read off the fibre configuration:
 * the poles of Jcal (with multiplicity) are the I_n / I_n^* fibres and their n.
 */
default(parisizemax, 8000000000);

jdetail(rn, rd, mm, j1, j2, aa, cc, dd, ff, cr) =
{ my(v, h, gam, ce, uu, vv, dg, fac, tot);
  v = jtest5(rn, rd, mm, j1, j2, aa, cc, dd, ff, cr);
  if(v == 0, print("  no elliptic-surface structure"); return(0));
  h = v[1]; gam = v[3];
  ce = certify(rn, rd, mm, j1, j2, aa, cc, dd, ff, cr, h, gam, DMAX, NTERM);
  uu = ce[4]; vv = ce[5]; dg = ce[1];
  print("  h = ", h, "   deg J = ", dg, "   gamma = ", gam);
  print("  U = ", uu);
  print("  V = ", vv);
  fac = factor(vv);
  print("  poles of J (finite):");
  tot = 0;
  for(i = 1, matsize(fac)[1],
      print("     ", fac[i,1], "  multiplicity ", fac[i,2]);
      tot += poldegree(fac[i,1])*fac[i,2]);
  print("  order of pole at t = infinity : ", dg - poldegree(vv),
        "   (deg U = ", poldegree(uu), ", deg V = ", poldegree(vv), ")");
  print("  zeros of J (j = 0, type II/IV/IV*/II* points): ", factor(uu));
  print("  J - 1728 : ", factor(uu - 1728*vv));
  v;
}
