default(parisizemax, 8000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/lib.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/census_util.gp");
fn = "/home/ubuntu/code/math-modular-sources/lattice/multi_prime/out/census_7_summary.log";
logit(fn, "### wide-window slope estimates  sigma_p = (v(n2)-v(n1))/(n2-n1) ###");
sw(nm, rr, p, n1, n2) = {
  my(qv = rr[1], pv = rr[2],
     v1 = valuation(pv[n1+1]/qv[n1+1] - pv[n1]/qv[n1], p),
     v2 = valuation(pv[n2+1]/qv[n2+1] - pv[n2]/qv[n2], p));
  logit(fn, Str("  ", nm, "  p=", p, "  v(", n1, ")=", v1, " v(", n2, ")=", v2,
                "  sigma=", (v2 - v1)*1.0/(n2 - n1)));
};
NN = 400;
sw("Domb  row3(10,4,64)", row3(10,4,64,NN), 2, 100, 399);
sw("T     row3(12,4,16)", row3(12,4,16,NN), 2, 100, 399);
sw("A row2(7,2,-8)",  row2(7,2,-8,NN),  2, 100, 399);
sw("B row2(9,3,27)",  row2(9,3,27,NN),  3, 100, 399);
sw("C row2(10,3,9)",  row2(10,3,9,NN),  3, 100, 399);
sw("E row2(12,4,32)", row2(12,4,32,NN), 2, 100, 399);
sw("F row2(17,6,72)", row2(17,6,72,NN), 2, 100, 399);
sw("F row2(17,6,72)", row2(17,6,72,NN), 3, 100, 399);
sw("s18 cooper18",    cooper18(NN),     3, 100, 399);
sw("delta row3(7,3,81)",   row3(7,3,81,NN),   3, 100, 399);
sw("zeta  row3(9,3,-27)",  row3(9,3,-27,NN),  3, 100, 399);
sw("eta   row3(11,5,125)", row3(11,5,125,NN), 5, 100, 399);
{ my(zz = zudrow(300));
  sw("zudrow", zz, 2, 100, 299);
  logit(fn, Str("  zudrow kappa_2 wide: ",
     (valuation(denominator(zz[1][300]),2) - valuation(denominator(zz[1][101]),2))*1.0/199)); }
logit(fn, "DONE");
quit;
