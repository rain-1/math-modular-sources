default(parisizemax, 6000000000);
default(realprecision, 1100);
LOG = "/home/ubuntu/code/math-modular-sources/lattice/multislope/row2_archid.log";
W(s) = write(LOG, s);
read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row2_arch.txt");
XI = XIARCH;
tg = [zeta(3), zeta(5), zeta(7), Pi^2, Pi^3, Pi^4, Pi^5, Pi^6, Pi^2*zeta(3), Catalan, log(2.0), log(2.0)^2, log(2.0)^3, log(13.0), log(17.0), sqrt(17.0), sqrt(13.0), zeta(3)^2, log(2.0)^4];
tn = ["zeta3","zeta5","zeta7","pi^2","pi^3","pi^4","pi^5","pi^6","pi^2 zeta3","G","log2","log2^2","log2^3","log13","log17","sqrt17","sqrt13","zeta3^2","log2^4"];
NM = ["C","D"];
W("=== AESZ207: are xi_inf(C), xi_inf(D) rational multiples of standard constants? (1100 digits) ===");
W("   [heights are max #decimal digits of the lindep coefficients; spurious level ~ 550 for 2 terms]");
{ for(s=2,3,
    for(i=1,#tg,
      my(L = lindep([XI[s], tg[i]]), h = max(#Str(abs(L[1])), #Str(abs(L[2]))));
      if(h < 60, W(Str("  HIT? ", NM[s-1], " vs ", tn[i], ": ", L)))));
  W("  (no HIT lines above = every single-constant test is at spurious height)"); }
W("");
W("--- 3-term: xi = q0 + q1*C for C in the catalogue ---");
{ my(best = 10^9, bi = 0, bs = 0);
  for(s=2,3, for(i=1,#tg,
    my(L = lindep([XI[s], 1.0, tg[i]]), h = vecmax(vector(3,j,#Str(abs(L[j])))));
    if(h < best, best = h; bi = i; bs = s)));
  W(Str("  smallest height over all (companion, constant) pairs: ", best, " digits, at ", NM[bs-1], " vs ", tn[bi]));
  W("  (spurious level for 3 terms at 1100 digits is ~367 digits)"); }
W("");
W("--- algdep of xi(C), xi(D), and of the ratios, degree 2 and 4 ---");
{ for(s=2,3,
    W(Str("  ", NM[s-1], ": algdep deg2 height ", #Str(vecmax(abs(Vec(algdep(XI[s],2))))),
          "   deg4 height ", #Str(vecmax(abs(Vec(algdep(XI[s],4))))))));
  W(Str("  xi(C)/xi(B): algdep deg2 height ", #Str(vecmax(abs(Vec(algdep(XI[2]/XI[1],2))))),
        "   xi(D)/xi(B): ", #Str(vecmax(abs(Vec(algdep(XI[3]/XI[1],2)))))));
  W(Str("  xi(C)/xi(B) = ", XI[2]/XI[1]*1.0));
  W(Str("  xi(D)/xi(B) = ", XI[3]/XI[1]*1.0)); }
W("DONE"); quit;
