default(parisizemax, 8000000000);
default(realprecision, 60);
read("/home/ubuntu/code/math-modular-sources/lattice/euler_criterion/lp.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/lib12.gp");
\\ Theorem-F prediction for the conductor-12 Catalan row:
\\   psi = chi_{-4}, w = 1, r_infty = 1.
\\   p=2: omega = character mod 4 = chi_{-4}, so chi_{-4} omega^{-1} = 1,
\\        Lam_2 = zeta_2(2) = L_2(2,1); E_2(2) = 1 - chi_{-4}(2)/4 = 1; xi_2 = zeta_2(2).
\\   p=3: omega = chi_{-3}, so chi_{-4} omega^{-1} = chi_12 (even, conductor 12),
\\        Lam_3 = L_3(2,chi_12); E_3(2) = 1 - chi_{-4}(3)/9 = 10/9; xi_3 = (9/10) L_3(2,chi_12).
{
my(PR2 = 900, PR3 = 400, z2, l3, x2, x3, ctl);
z2 = Lp(2, triv,  2, PR2);
l3 = Lp(3, chi12, 2, PR3);
print("zeta_2(2)          = ", z2 + O(2^40));
print("L_3(2,chi_12)      = ", l3 + O(3^30));
print("pred xi_2 = zeta_2(2);  pred xi_3 = (9/10) L_3(2,chi_12) = ", (9/10)*l3 + O(3^30));
print("");
print("row     b   v2(xi_b - zeta_2(2))   v3(xi_b - (9/10)L_3(2,chi12))   Cauchy v2   Cauchy v3");
foreach([[1,1],[1,2],[3,2],[2,1]], fr,
  my(nu0 = fr[1], de0 = fr[2]);
  foreach([24,36,48], bb,
    if(bb % de0, next);
    my(aa = nu0*bb/de0, r, rm, x, xm, d2, d3, c2, c3);
    if(aa < 1 || aa > 4*bb, next);
    r  = c12row(aa, bb);      x  = r[2]/r[1];
    rm = c12row(aa - nu0/de0*1, bb - 1); xm = rm[2]/rm[1];
    d2 = valuation((x + O(2^PR2)) - z2, 2);
    d3 = valuation((x + O(3^PR3)) - (9/10)*l3, 3);
    c2 = valuation(x - xm, 2); c3 = valuation(x - xm, 3);
    printf("%d/%d  %3d   %8d   %8d      %8d  %8d\n", nu0,de0,bb, d2, d3, c2, c3)));
print("");
print("--- controls: wrong rational scalar (row alpha=1, b=36) ---");
my(r = c12row(36,36), x = r[2]/r[1]);
foreach([1/2, 1, 2, 4, 5/4, 4/5, 8/5, 5/8, 9/10, 10/9, 3/2, 2/3], sc,
  printf("  v2(xi - %s*zeta_2(2)) = %4d      v3(xi - %s*L_3(2,chi12)) = %4d\n",
    sc, valuation((x+O(2^PR2)) - sc*z2, 2), sc, valuation((x+O(3^PR3)) - sc*l3, 3)));
}
quit;
