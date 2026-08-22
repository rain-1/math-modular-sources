/* ============================================================
   06_units_and_level12.gp
   (a) identify the linear gauge factors 1-t_C, 1-9t_C as eta quotients
       (hence as modular units, and read off which cusp goes to t=infinity);
   (b) cross-check against the level-12 cover picture of RIGIDITY_PROOF.md.
   ============================================================ */
read("lib.gp");
E1 = etaprod(1); E2 = etaprod(2); E3 = etaprod(3); E6 = etaprod(6);
E4 = etaprod(4); E12 = etaprod(12);
tC = 'q * E1^4 * E6^8 / (E2^8 * E3^4);
FC = E2^6 * E3 / (E1^3 * E6^2);
tFz = 'q * E1^5 * E3 * E4^5 * E6^2 * E12 / E2^14;

/* fit f = prod_{d|6} E_d^{r_d} by matching log-coefficients */
{ etafit(f, dv) =
  my(n = #dv, L = log(f + O('q^40)), rows = 20);
  my(Ms = matrix(rows, n, i, j, polcoef(log(etaprod(dv[j],40)), i, 'q)));
  my(rhs = vectorv(rows, i, polcoef(L, i, 'q)));
  Vec(matinverseimage(Ms, rhs));
}
dv = [1,2,3,6];
r1 = etafit(1-tC, dv);   r2 = etafit(1-9*tC, dv);   r3 = etafit(FC, dv);
print("=== A. eta-quotient exponents (basis E_1,E_2,E_3,E_6) ===");
print("1 - t_C   ~ ", r1);
print("1 - 9 t_C ~ ", r2);
print("F_C       ~ ", r3);
chk(f, r) = f - E1^r[1]*E2^r[2]*E3^r[3]*E6^r[4];
print("check 1-t_C   : ", chk(1-tC, r1));
print("check 1-9 t_C : ", chk(1-9*tC, r2));
print("check F_C     : ", chk(FC, r3));

/* Ligozat order of vanishing of prod eta_d^{r_d} at the cusp 1/c' of Gamma_0(6) */
{ ligozat(r, dv, NL, cp) =
  my(s = 0);
  for(j=1,#dv, my(d=dv[j]); s += gcd(cp,d)^2*r[j]/(gcd(cp,NL\cp)*cp*d));
  NL*s/24;
}
print();
print("=== B. Ligozat orders at the four cusps of X_0(6) (c'=1:cusp 0, 2:1/2, 3:1/3, 6:oo) ===");
{ show(nm, r) = print(nm, ":  ord_0=", ligozat(r,dv,6,1), "  ord_{1/2}=", ligozat(r,dv,6,2),
                       "  ord_{1/3}=", ligozat(r,dv,6,3), "  ord_oo=", ligozat(r,dv,6,6)); }
show("t_C      ", [4,-8,-4,8]);
show("1-t_C    ", r1);
show("1-9t_C   ", r2);
show("F_C      ", r3);
print("cusp widths of Gamma_0(6):  0 -> 6,  1/2 -> 3,  1/3 -> 2,  oo -> 1");

print();
print("=== C. link with the level-12 cover picture of RIGIDITY_PROOF.md ===");
u = tC; v = Vop(2,tC);
print("modular equation  v^2-(9u^2-8u+1)v+u^2   : ", v^2-(9*u^2-8*u+1)*v+u^2);
print("t_F - (v-u)/(9v-1)                        : ", tFz - (v-u)/(9*v-1));
print("t_F' - ( -(v-u)/(9v-1) )|_{tau->tau+1/2}  : ", tC/(1-9*tC) + twist((v-u)/(9*v-1)));
quit;
