LB = 3^(-2)*(zetahurwitz(2,1/3) - zetahurwitz(2,2/3));
print("L(chi_-3,2) = ", LB);

\\ ---------------- V1 : Zagier B ----------------
tcB = (9 - I*sqrt(27))/54;
xiB = dofold("V1 Zagier B", 9,9,3, 27,0,0, tcB, 2500, 900, 2400);
knB = LB/2 + I*2*Pi^2/(27*sqrt(3));
print("V1 known xi_B      = ", knB);
print("V1 |xi - known|    = ", abs(xiB - knB));
print("V1 |conj(xi)-known|= ", abs(conj(xiB) - knB));

\\ ---------------- V2 : Zagier E ----------------
tcE = 1/8 + 0.*I;
xiE = dofold("V2 Zagier E", 12,12,4, 32,0,0, tcE, 2500, 900, 2400);
G = Catalan;
print("V2 known xi_E = G/2 = ", G/2);
print("V2 |xi - G/2|       = ", abs(xiE - G/2));

\\ ---------------- V3 : Zagier D ----------------
tcD = (sqrt(125)-11)/2 + 0.*I;
xiD = dofold("V3 Zagier D", 11,11,3, -1,0,0, tcD, 2500, 900, 2400);
print("V3 known xi_D = zeta(2)/5 = ", zeta(2)/5);
print("V3 |xi - zeta(2)/5|       = ", abs(xiD - zeta(2)/5));
