read("/home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/cdt_ab.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/cdt_audit.gp");

NMAX = eval(getenv("NMAX")); if(NMAX == 0, NMAX = 60);
{
docase2(r, e, nmax, tag) =
  print("=== ", tag, "  expected sharp LCM exponent = ", e, " ===");
  print("  b_0..b_6 = ", vector(min(7,nmax+1), j, r[2][j]));
  print("  a_0..a_6 = ", vector(min(7,nmax+1), j, r[1][j]));
  audit(r[1], r[2], e, nmax, concat("  ", tag));
}
docase2(run14(NMAX), 2, NMAX, "X_1(4)=X_0(4), chi_{-4}  (2-adic Catalan)");
docase2(run09(NMAX), 2, NMAX, "X_0(9), chi_{-3}  (3-adic analogue)");
docase2(run03(NMAX), 2, NMAX, "X_0(3), chi_{-3}  (level-3 comparison)");
quit;
