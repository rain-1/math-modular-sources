{
show(lab, f, ps) =
  my(mf, co, ap, D, u);
  mf = mfinit(f, 0);
  print("--- ", lab, "  dim=", mfdim(mf));
  f = mfeigenbasis(mf)[1];
  co = mfcoefs(f, 30);
  for(i=1,#ps, my(p=ps[i]);
    ap = co[p+1];
    print("   p=",p,"  a_p=",ap);
  );
}
/* Apery zeta(3): eta(2z)^4 eta(4z)^4 in S_4(Gamma_0(8)) */
show("S_4(Gamma_0(8)) [Apery zeta(3), Beukers]", [8,4,1], [5,7,11,13]);
/* weight-3 CM newform level 12 character chi_{-3}: eta(2z)^3 eta(6z)^3 */
show("S_3(Gamma_0(12), chi_-3)", [12,3,-3], [5,7,11,13]);
/* weight-3 level 8 chi_{-4}? and level 5 quartic for row D */
show("S_3(Gamma_0(8),chi_-4)?", [8,3,-4], [5,7,11,13]);
quit
