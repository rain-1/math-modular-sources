/* 07_measures.gp -- irrationality measures and Nesterenko exponents for the
   cusp-move placements.  All numbers reported in CUSP_MOVE_PROGRAM.md Sections 5-6. */
default(realprecision, 40);

/* single-row measure:  mu <= 1 + (k + log lam1)/(-k - log|lam2|)  when score > 0 */
{ mumeas(nm, l1, l2, k) =
  my(sg = k + log(abs(l1)), dl = -k - log(abs(l2)));
  print("  ", nm, ":  lam1=", l1, "  lam2=", l2, "  k=", k);
  print("      sigma = k + log lam1     = ", sg);
  print("      delta = -k - log|lam2|   = ", dl, "   (= the score)");
  if(dl > 0, print("      mu <= 1 + sigma/delta   = ", 1 + sg/dl),
             print("      delta <= 0: no irrationality, no measure"));
  [sg, dl]; }

print("=== 5.1  single-row irrationality measures, positive-score rows ===");
sq5 = sqrt(5); sq2 = sqrt(2);
lamD = (11+5*sq5)/2;  muD = (11-5*sq5)/2;
print(" Apery zeta(2) row = Zagier D (11,3,-1), placement at v=0:");
rD = mumeas("  D.v=0", lamD, muD, 2);
print(" its Galois-conjugate placement at v=mu (over Q(sqrt5)):");
rDp = mumeas("  D.v=mu", lamD-muD, -muD, 2);
print("");
lamB = 4*(17+12*sq2); muB = 4*(17-12*sq2);
print(" Beukers (136,10,4) = sqrt(Apery), placement at v=0:");
rB = mumeas("  B.v=0", lamB, muB, 2);
print(" its Galois-conjugate placement at v=mu (over Q(sqrt2)):");
rBp = mumeas("  B.v=mu", lamB-muB, -muB, 2);

print("");
print("=== 5.2  the two-sequence (merged) bound  mu <= 1 + min_i sigma_i/delta_i ===");
{ print(" Apery zeta(2):   single ", 1+rD[1]/rD[2], "   merged with conjugate ", 1 + min(rD[1]/rD[2], rDp[1]/rDp[2])); }
{ print(" Beukers:         single ", 1+rB[1]/rB[2], "   merged with conjugate ", 1 + min(rB[1]/rB[2], rBp[1]/rBp[2])); }
print(" gain (Apery)  = ", (1+rD[1]/rD[2]) - (1 + min(rD[1]/rD[2], rDp[1]/rDp[2])));
print(" gain (Beukers)= ", (1+rB[1]/rB[2]) - (1 + min(rB[1]/rB[2], rBp[1]/rBp[2])));

print("");
print("=== 5.3  the norm obstruction for the conjugate placement ===");
print(" Apery zeta(2): |lam2| at v=mu and at v=lam are ", abs(muD), " and ", abs(lamD));
print("   product (= |d|) = ", abs(muD*lamD), "   log = ", log(abs(muD*lamD)));
{ print("   norm rate of the linear form with d_n^2 at both places = ", log(abs(muD*lamD)) + 4, "  > 0  : no descent to Q"); }
print(" Beukers: |lam2| at v=mu and at v=lam are ", abs(muB), " and ", abs(lamB-muB));
print("   product = ", abs(muB*(lamB-muB)), "   log = ", log(abs(muB*(lamB-muB))));
{ print("   norm rate with d_n^2 at both places = ", log(abs(muB*(lamB-muB))) + 4, "  > 0 : no descent to Q"); }

print("");
print("=== 5.4  the only same-period placement pair: Zagier C and F' ===");
print(" C  : lam1=9, lam2=1,  k=2,  xi = L(2,chi-3)/2");
print(" F' : lam1=9, lam2=8,  k=2,  xi = -5 L(2,chi-3)/8   (ratio -5/4, exact)");
{ my(sC = 2+log(9), dC = -2-log(1), sF = 2+log(9), dF = -2-log(8));
  print("  sigma_C=", sC, "  delta_C=", dC, "     sigma_F=", sF, "  delta_F=", dF);
  print("  both delta < 0: the merged bound  1 + min(sigma/delta)  is vacuous.");
  print("  deficit to delta_C > 0 : ", -dC, " nats/step;  to delta_F > 0 : ", -dF); }

print("");
print("=== 6  Nesterenko exponent for (1, zeta(2), L(2,chi-3)) from the pair (A, C) ===");
{ my(sA = 2+log(8), dA = -2-log(1), sC = 2+log(9), dC = -2-log(1), sg, tau);
  print(" row A : coefficient rate sigma_A = k + log lam1 = ", sA,
        "   form rate log|L_n|/n = ", -dA);
  print(" row C : coefficient rate sigma_C = ", sC, "   form rate = ", -dC);
  sg = max(sA, sC); tau = min(dA, dC)/sg;
  print(" joint  sigma = ", sg, "   tau = delta/sigma = ", tau);
  print(" Nesterenko: dim_Q <1, zeta(2), L(2,chi-3)> >= 1 + tau = ", 1+tau, "   (vacuous)");
  print(" needed for dim >= 2 : tau > 0, i.e. delta > 0        -> deficit ", -min(dA,dC), " nats/step");
  print(" needed for dim >= 3 : tau >= 2, i.e. delta >= 2 sigma -> deficit ",
        2*sg - min(dA,dC), " nats/step");
  print("");
  print(" best-possible (budget) version: replace |lam2| by 1/lam1 (Fricke-perfect):");
  my(bA = log(8)-2, bC = log(9)-2, tb);
  print("   budget_A = log lam1 - k = ", bA, "   budget_C = ", bC);
  tb = min(bA,bC)/sg;
  print("   tau_budget = ", tb, "   dim >= ", 1+tb, "  (> 1: would give ONE irrationality)");
  print("   still short of dim >= 3 by ", 2 - tb, " in tau, i.e. ",
        (2-tb)*sg, " nats/step"); }

print("");
print("=== 6.2  simultaneous approximation with a common denominator ===");
{ my(s0 = 2 + log(8) + log(9), e1 = 2 + log(9), e2 = 2 + log(8));
  print(" common denominator q_n = d_n^2 a_n^A a_n^C : log q_n / n = ", s0);
  print("   |q_n zeta(2) - p_n^(1)|  rate = ", e1, "  (grows)");
  print("   |q_n L(2,chi-3) - p_n^(2)| rate = ", e2, "  (grows)");
  print("   simultaneous-approximation exponent = -max(e1,e2)/s0 = ", -max(e1,e2)/s0,
        "  (need > 1/2 for the Nesterenko dual criterion in dimension 2)"); }
quit;
