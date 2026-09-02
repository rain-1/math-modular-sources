/* 04_polar.gp -- (P1): rho = (2 pi i)^{r+1} a_0(Phi|_{r+1} gamma) / (r! c^{r+1}). */
read("lib.gp");
default(realprecision, 60);
cuspmat(a,c) = my(g=bezout(a,c)); [a, -g[2]; c, g[1]];
G3 = znstar(3,1); G4 = znstar(4,1); G5 = znstar(5,1);
CH3 = [G3,[1]]; CH4 = [G4,[1]]; CHQ5 = [G5,[2]]; PS = [G5,[1]]; PSB = [G5,[3]];
/* our E^{psi,phi} = mfeisenstein(r+1, phi, psi) */
mkE(r, phi, psi) = mfeisenstein(r+1, phi, psi);
mkS(r, phi, psi, dl, cl) = my(cp=vector(#dl, i, mfbd(mkE(r,phi,psi), dl[i]))); mflinear(cp, cl);
chk(nm, src, FF) = {
  my(r=src[2], mf=mfinit(FF), cus, mt, ex, a0, pred, b, ss);
  print("== ", nm);
  cus = mfcusps(mf);
  for(i=1,#cus, my(cu=cus[i]);
    if(type(cu)=="t_INFINITY", next());
    my(aa=numerator(cu), cc=denominator(cu));
    mt = cuspmat(aa,cc);
    ex = mfslashexpansion(mf, FF, mt, 3, 0);
    a0 = ex[1];
    b = cperiod(src, aa, cc);
    pred = (2*Pi*I)^(r+1)*a0/(r!*cc^(r+1));
    printf("   cusp %s/%s : a_0 = %s\n", aa, cc, a0);
    printf("       rho(computed) = %s\n", b[2]);
    printf("       (2 pi i)^{r+1} a_0/(r! c^{r+1}) = %s     diff = %s\n", pred, abs(pred-b[2])));
}
chk("A  (outer chi_-3, level 6)", srcbyname("A"), mkS(2, CH3, 1, [1,2], [1,-1]));
chk("C  (inner chi_-3, level 6)", srcbyname("C"), mkS(2, 1, CH3, [1,2], [1,-8]));
chk("E  (inner chi_-4, level 8)", srcbyname("E"), mkS(2, 1, CH4, [1,2], [1,-8]));
chk("gamma (trivial/trivial, level 6)", srcbyname("gamma"), mkS(3, 1, 1, [1,2,3,6], [1,-28,63,-36]));
chk("zeta (chi_-3 both, level 9)", srcbyname("zeta"), mkS(3, CH3, CH3, [1], [1]));
/* Gamma_1(5): Phi_D = (1/2+i)E^{1,psi4} + (1/2-i)E^{1,psibar4} */
PD = mflinear([mkE(2,PS,1), mkE(2,PSB,1)], [1/2+I, 1/2-I]);
print("== D  (outer w, Gamma_1(5))   coefs: ", mfcoefs(PD,6));
chk("D", srcbyname("D"), PD);
quit;
