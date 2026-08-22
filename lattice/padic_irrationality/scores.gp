\\ scores.gp -- assemble the p-adic irrationality score
\\    S_p = (v_p(c) + kappa_p) log p - k - log lambda_1
\\    theta_p = sigma_p log p / (k + kappa_p log p + log lambda_1),  sigma_p = v_p(c) + 2 kappa_p
\\ for every census cell with a proved/verified p-adic Apery limit, plus the
\\ Calegari calibration cells.  Data provenance is in the markdown.
\\ Run: gp -q /home/ubuntu/code/math-modular-sources/lattice/padic_irrationality/scores.gp
default(realprecision,25);

\\ cell = [name, p, sigma_p, kappa_p, k, log lambda_1, xi_p label]
{ CELLS = [
 \\ ---- Calegari 2005 calibration (his own rows) ----
 ["Cal zeta_2(3)  X_0(2)",   2, 12, 0, 3, 6*log(2), "zeta_2(3)/2"],
 ["Cal zeta_3(3)  X_0(3)",   3,  6, 0, 3, 3*log(3), "zeta_3(3)/2"],
 ["Cal zeta_5(3)  X_0(5)",   5,  3, 0, 3, 1.5*log(5), "zeta_5(3)/2"],
 ["Cal zeta_7(3)  X_0(7)",   7,  2, 0, 3, 1.0*log(7), "zeta_7(3)/2"],
 ["Cal zeta_13(3) X_0(13)", 13,  1, 0, 3, 0.5*log(13),"zeta_13(3)/2"],
 ["Cal zeta_2(5)  X_0(2)",   2, 12, 0, 5, 6*log(2), "zeta_2(5)/2"],
 ["Cal zeta_2(7)  X_0(2)",   2, 12, 0, 7, 6*log(2), "zeta_2(7)/2"],
 ["Cal Catalan    X_1(4)",   2,  8, 0, 2, 4*log(2), "zeta_2(2)/2 = L_2(2,chi_-4)/2"],
 \\ ---- Zagier order-2 census (k=2) ----
 ["Zagier A",  2, 3, 0, 2, log(8),      "0  (degenerate)"],
 ["Zagier B",  3, 3, 0, 2, log(27)/2,   "zeta_3(2)/2"],
 ["Zagier C",  3, 2, 0, 2, log(9),      "zeta_3(2)/2"],
 ["Zagier E",  2, 5, 0, 2, log(8),      "zeta_2(2)/2"],
 ["Zagier F",  2, 3, 0, 2, log(9),      "L_2(2,chi_12)/2"],
 ["Zagier F",  3, 2, 0, 2, log(9),      "5 zeta_3(2)/8"],
 \\ ---- AZ order-3 census (k=3) ----
 ["AZ delta",  3, 4, 0, 3, log(9),      "zeta_3(3)/4"],
 ["AZ zeta",   3, 3, 0, 3, log((18+sqrt(18^2+4*27))/2), "0 (degenerate)"],
 ["Domb",      2, 6, 0, 3, log(16),     "zeta_2(3)/3"],
 ["AZ eta",    5, 3, 0, 3, log(125)/2,  "zeta_5(3)/2"],
 ["T = eps",   2, 4, 0, 3, log((24+sqrt(24^2-4*16))/2), "zeta_2(3)/4"],
 \\ ---- Cooper ----
 ["Cooper s18",3, 1, 0, 2, log(16),     "zeta_3(2)/2"],
 \\ ---- other verified cells ----
 ["cusp L(f,2) lvl12", 2, 2, 0, 2, log(16), "0 (degenerate)"],
 ["zeta(5) level 16 I", 2, 1, 0, 5, log(2+4*sqrt(2)), "7 zeta_2(5)/32"],
 ["zeta(5) level 16 II",2, 1, 0, 5, log(4),            "7 zeta_2(5)/32"],
 ["AESZ 207 (rank 4)",  2,12, 0, 4, log(89531.3893), "unidentified"],
 \\ ---- hypergeometric rows with p-power denominators (kappa > 0) ----
 ["Zudilin Catalan",    2, 8, 4, 4, 5*log((1+sqrt(5))/2), "zeta_2(2)"],
 ["Nesterenko (4,7)",   2,28,14,12, 7.650713334239925,    "zeta_2(2)"]
]; }

print("cell | p | sigma_p | kappa_p | k | log lambda_1 | S_p | theta_p");
{
my(rows=List());
for(i=1,#CELLS,
  my(cc=CELLS[i], nm=cc[1], pp=cc[2], sg=cc[3], kap=cc[4], kk=cc[5], l1=cc[6], id=cc[7]);
  my(den = kk + kap*log(pp) + l1);
  my(Sp = sg*log(pp) - den, th = sg*log(pp)/den);
  listput(rows,[Sp, nm, pp, sg, kap, kk, l1, th, id]);
  printf("%-22s | %2d | %2d | %2d | %2d | %8.5f | %+9.5f | %8.6f | %s\n",
         nm, pp, sg, kap, kk, l1, Sp, th, id);
);
print("\n--- sorted by S_p, descending ---");
rows = vecsort(Vec(rows),1,4);
for(i=1,#rows, my(r=rows[i]);
  printf("%+9.5f  theta=%8.6f  %-22s p=%d   xi_p = %s\n", r[1], r[8], r[2], r[3], r[9]));
}

print("\n=== cusp-move placement scan (order-2 systems) ===");
print("row : placement (roots) -> c, lambda_1, S_p at each p | c");
{
my(D=[["A",7,-8],["B",9,27],["C",10,9],["D",11,-1],["E",12,32],["F",17,72]]);
for(i=1,#D,
  my(nm=D[i][1], a=D[i][2], c=D[i][3], disc=a^2-4*c);
  my(lam, mu);
  lam=(a-sqrt(disc+0.))/2; mu=(a+sqrt(disc+0.))/2;
  print("\n", nm, "  roots ", lam, ", ", mu);
  my(PL=[[lam,mu],[-lam,mu-lam],[-mu,lam-mu]]);
  for(j=1,3, my(r=PL[j], cc=r[1]*r[2], l1=max(abs(r[1]),abs(r[2])));
    my(ci=round(real(cc)));
    print("   placement ",j,": c = ",ci,"   lambda_1 = ",l1,
      "   S_p: ", concat(vector(4,t, my(pp=[2,3,5,7][t]);
        if(ci!=0 && valuation(ci,pp)>0,
           Str("  p=",pp,": ",valuation(ci,pp)*log(pp)-2-log(l1)), ""))));
  );
);
}
quit;
