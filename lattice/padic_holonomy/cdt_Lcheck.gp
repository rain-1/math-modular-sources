/* Identify eta for the odd-character cases:
      L_p(2,chi) = lim L(-2K,chi) = lim -B_{2K+1,chi}/(2K+1),  2K+1 == -1 mod (p-1)p^N.
   (chi(p)=0 for both (chi_{-4},p=2) and (chi_{-3},p=3), so no Euler factor.)
   Claim: 2*eta = L_p(2,chi)  with eta = -a_m/b_m.                              */
read("/home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/cdt_ab.gp");
{
bchi(n, D) = my(f = abs(D)); f^(n-1)*sum(a = 1, f, kronecker(D, a)*subst(bernpol(n), x, a/f));
}
{
chk(r, p, D, m, nmax, tag) =
  my(et, mm, zz);
  et = -r[1][m+1]/r[2][m+1];
  print(tag, ":  eta = -a_", m, "/b_", m, ";  2*eta should be L_", p, "(2,chi_", D, ")");
  print("     2*eta = ", 2*et);
  for(N = 2, nmax,
    mm = (p-1)*p^N - 1;
    zz = -bchi(mm, D)/mm;
    print("      N=", N, " (2K+1=", mm, ")   v_", p, "(2*eta - L) = ",
          if(2*et - zz == 0, "exact", valuation(2*et - zz, p)),
          "     [wrong sign: v_p(2*eta + L) = ", if(2*et + zz == 0, "exact", valuation(2*et + zz, p)), "]"));
  print("");
}
chk(run14(60), 2, -4, 60, 7, "X_1(4) chi_{-4}  (2-adic Catalan)");
chk(run09(60), 3, -3, 60, 4, "X_0(9) chi_{-3}  (3-adic)");
chk(run03(60), 3, -3, 60, 4, "X_0(3) chi_{-3}  (3-adic, level 3)");
quit;
