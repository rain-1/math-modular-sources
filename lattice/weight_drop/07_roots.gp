default(parisizemax,4000000000);
NQ = 90;
read("lattice/weight_drop/03_setup.gp"); build();
lam=[1,2,2];
{for(r=1,#DAT, my(nm=DAT[r][1],T=DAT[r][3],F=DAT[r][4]);
  my(g=sqrt(F), PSI=g^3*T/lam[r]);
  print("=== ",nm);
  print("  v_2(den [q^m] g), m=1..60: ", vector(60,m, valuation(denominator(polcoeff(g,m)),2)));
  print("  v_2(den [q^m] Psi_root), m=1..40: ", vector(40,m, valuation(denominator(polcoeff(PSI,m)),2)));
);}
/* s7: Psi_root = eta_1^3 eta_7^3 ? */
{my(E17 = q*eta(q+O(q^NQ))^3*eta(q^7+O(q^NQ))^3, g=sqrt(DAT[1][4]), PSI=g^3*DAT[1][3]);
 print("s7:  Psi_root - eta_1^3 eta_7^3 = ", PSI - E17);
 print("s7:  g - E_1(chi_-7) = ", g - (1+2*sum(n=1,NQ-1, sumdiv(n,d,kronecker(-7,d))*q^n)+O(q^NQ)));}
quit;
