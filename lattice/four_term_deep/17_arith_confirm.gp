default(parisizemax, 24000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/17_arith_lib.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/mum_survey/lpgen.gp");
fn = "/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/arith_confirm.log";
logit(s) = { print(s); write(fn, s); };
NN = 400;
PR = 1200;

/* targets, project normalisation: Lam_p(chi) = L_p(2, chi*om^-1) */
Lm4 = LpG(2, chim4, -1, 2, PR);   \\ = zeta_2(2) = Calegari's L_2(2,chi_-4)
Lm3 = LpG(2, chim3, -1, 2, PR);   \\ = L_2(2,chi_12)
logit(Str("zeta_2(2) = L_2(2,chi-4 om^-1) mod 2^24 : ", lift(Lm4+O(2^24))));
logit(Str("L_2(2,chi12)=L_2(2,chi-3 om^-1) mod 2^24: ", lift(Lm3+O(2^24))));

{PRED = [[1/4, Lm4, "1/4 * zeta_2(2)"], [1/2, Lm4, "1/2 * zeta_2(2)"],
         [-1/2, Lm4, "-1/2 * zeta_2(2)"], [3/8, Lm3, "3/8 * L_2(2,chi12)"],
         [3/4, Lm3, "3/4 * L_2(2,chi12)"], [0, Lm4, "(no limit)"]];}

chk(nm, av, bv, pr) = {
  my(xi, d, vv, cas);
  d = bv[NN+1]/av[NN+1] - bv[NN]/av[NN];
  vv = if(d==0, -1, valuation(d,2));
  logit(Str("--- ", nm, "   v_2(incr_N)=", vv));
  if(pr[1] == 0, logit("    no 2-adic limit predicted"); return);
  xi = bv[NN+1]/av[NN+1] + O(2^min(PR, vv-4));
  my(r = xi - pr[1]*pr[2]);
  logit(Str("    xi_2 - (", pr[3], ") : ", if(r==0, Str("ZERO to ", min(PR,vv-4), " digits"), Str("v_2 = ", valuation(r,2)))));
};

/* 2-adic valuation of Cas_n over a window, both parities */
casprof(nm, av, bv) = {
  my(cas = vector(60));
  for(n=1,60, cas[n] = av[n+1]*bv[n+2] - av[n+2]*bv[n+1]);
  logit(Str("    v_2(Cas_n), n=41..60: ", vector(20,i, valuation(cas[40+i],2))));
  logit(Str("    v_2(Cas_n)-9n/2 or -7n/2 etc: (raw above)"));
};

{for(i=1,#ROWS, my(rr = row4(ROWS[i][2], NN)); chk(ROWS[i][1], rr[1], rr[2], PRED[i]); casprof(ROWS[i][1], rr[1], rr[2]));}
{my(rr = zrow2(12,4,32,NN)); chk("Zagier E (12,4,32)", rr[1], rr[2], [1/2, Lm4, "1/2 * zeta_2(2)"]); casprof("E", rr[1], rr[2]);}
logit("DONE");
quit;
