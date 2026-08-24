/* 09_prop.gp -- checks for Propositions D1/D2 of consolidation/FOUR_TERM_DEEP.md.
 * D1: a four-term row with |lam_2| < 1 has IRREDUCIBLE characteristic cubic
 *     lam^3 - a lam^2 + d lam - g  (g != 0), hence (Cor F2.1) equal exponents.
 * D2: if in addition lam_2 is real then the cubic is totally real, disc > 0.  */
default(realprecision, 60);
CROWS = [[0,1,2,1,1,154,42,-128,-12,8],[0,1,2,1,1,87,24,-152,-16,16],[1,2,4,5,5,72,6,-192,-60,8],[0,1,3,-1,4,22,6,-72,12,36],[0,1,2,1,1,22,6,-96,-12,24],[0,1,2,1,1,11,3,-73,-8,28],[-1,2,4,-1,-1,16,10,64,20,4],[1,2,4,5,5,16,2,64,20,4],[0,1,2,1,1,60,18,240,12,16],[0,1,2,1,1,20,6,-80,-4,-16],[0,1,2,1,1,3,0,-104,-16,80],[0,1,3,1,2,11,4,37,3,3],[0,1,2,1,1,16,6,96,12,64],[0,1,2,1,1,16,6,76,8,24],[0,1,2,1,1,16,6,96,12,48],[0,1,4,1,3,14,5,97,8,12],[0,1,2,1,1,16,6,16,-4,-96],[0,1,6,-1,7,13,4,432,-24,-48],[0,1,2,1,1,6,2,64,4,-40],[0,1,3,1,2,3,0,-216,-24,108],[0,1,3,1,2,33,12,324,36,108],[0,1,3,-1,4,3,0,-216,48,108],[0,1,2,1,1,17,6,112,8,24],[0,1,4,1,3,28,10,208,12,16],[0,1,6,1,5,17,6,56,0,-12],[0,1,2,1,1,32,12,272,16,64],[0,1,2,1,1,42,15,441,24,100],[-1,3,1,0,0,42,22,480,64,-64],[1,3,1,1,1,42,8,480,64,-64],[0,1,2,1,1,64,24,1030,26,48]];
print("--- census rows with |lam_2| <= 1 ---");
chkrow(v) = my(m=v[3], a=v[6], d=v[8], gg=v[10]*m^2, ch, l, l2, irr, dsc); ch = x^3-a*x^2+d*x-gg; l = vecsort(polroots(ch), z->-abs(z)); l2 = abs(l[2]); irr = polisirreducible(ch); dsc = poldisc(ch); if(l2 <= 1.0000001, print("  |lam2|=", l2, "  irreducible=", irr, "  disc=", dsc, "  row=", v, if(irr, "", concat("  factors ", Str(factor(ch)))))); 0;
for(i=1, #CROWS, chkrow(CROWS[i]));
print();
print("--- D1 exhaustive numerical check over integer cubics ---");
BAD = 0; TOT = 0;
scanbox(AB, DB, GB) = my(ch, l, l2, irr); for(a = -AB, AB, for(d = -DB, DB, for(gg = -GB, GB, if(gg == 0, next); ch = x^3 - a*x^2 + d*x - gg; l = vecsort(polroots(ch), z->-abs(z)); l2 = abs(l[2]); if(l2 < 0.999999, TOT = TOT + 1; irr = polisirreducible(ch); if(!irr, BAD = BAD + 1; print("  COUNTEREXAMPLE a=", a, " d=", d, " g=", gg, " |lam2|=", l2)))))); 0;
scanbox(24, 70, 70);
print("  cubics with |lam_2|<1 in |a|<=24,|d|<=70,|g|<=70 : ", TOT, "   reducible among them: ", BAD);
quit;
