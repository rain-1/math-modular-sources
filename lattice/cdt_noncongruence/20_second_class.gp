/* Is there a SECOND fold-regular conditional class on the (Gamma_0(5)+5, t) host?
   t = u/(1+22u+125u^2), u = (eta_5/eta_1)^6, base F = E_{2,5}, lambda = 2.
   For each candidate weight-2 form F' with F'(0)=1: t-expansion, g'=sqrt(F'),
   minimal 2-power denominator lambda', minimal recurrence of a'_n=lambda'^n[t^n]g',
   characteristic roots -> singular t-values, and (if order 2) k' and xi'.
   PARI builtin names psi, M, Phi, S, cmp avoided.                            */
default(parisize, 3000000000);
default(realprecision, 60);
MQ = 80; NROW = 60; NK = 55; PGATE = 2147483647;
DL = vector(NK+2, i, if(i==1, 1, lcm(vector(i-1,j,j))));
etaprod(d, MM) = { my(s = 1 + O('q^MM)); forstep(n=d, MM-1, d, s *= (1 - 'q^n)); s; };
fitminp(A, maxo, maxd) = { my(L = #A - 1, Ap = vector(#A, i, Mod(A[i], PGATE))); for(o = 1, maxo, for(d = 0, maxd, my(nc = (o+1)*(d+1), rows = List(), r); if(L - o + 1 < nc + 6, next); for(n = o, L, r = vector(nc); for(j = 0, o, for(e = 0, d, r[j*(d+1)+e+1] = Mod(n,PGATE)^e * Ap[n+1-j])); listput(rows, r)); my(Mx = matconcat(Vec(rows)~), K = matker(Mx)); if(#K == 1, my(rows2 = List(), r2); for(n = o, L, r2 = vector(nc); for(j = 0, o, for(e = 0, d, r2[j*(d+1)+e+1] = n^e * A[n+1-j])); listput(rows2, r2)); my(M2 = matconcat(Vec(rows2)~), K2 = matker(M2)); if(#K2 == 1, return([o, d, K2[,1]]))))); 0; };
charpol_from(v, o, d) = { my(c = vector(o+1, j, v[(j-1)*(d+1)+d+1])); sum(j=1, o+1, c[j]*'L^(o+1-j)); };
compk(rec, NN) = { my(o=rec[1], d=rec[2], v=rec[3], b=vector(NN+1)); for(i=1,o, b[i]=0); b[o]=1; my(cf=vector(o+1, j, sum(e=0,d, v[(j-1)*(d+1)+e+1]*'n^e))); for(n=o, NN, my(s2=0); for(j=1,o, s2 += subst(cf[j+1],'n,n)*b[n+1-j]); my(lead=subst(cf[1],'n,n)); if(lead==0, return(-1)); b[n+1] = -s2/lead); my(k=0); while(k<=8, my(ok=1); for(n=0,NN, if(denominator(DL[n+1]^k*b[n+1])!=1, ok=0; break)); if(ok, break); k++); [k, b]; };

E1 = etaprod(1, MQ); E5 = etaprod(5, MQ);
uu = 'q*(E5/E1)^6 + O('q^MQ);
tq = uu/(1+22*uu+125*uu^2) + O('q^MQ);
qt = serreverse(tq);
print("t(q) = ", Ser(vector(7,i,polcoeff(tq,i-1)),'q));
print("q(t) integral: ", vecmax(vector(NROW,i,denominator(polcoeff(qt,i))))==1);
BL1 = 44+20*sqrt(5); BL2 = 44-20*sqrt(5);
print("base: lam_1 = ", BL1, ", lam_2 = ", BL2, ", lambda = 2, t_1 = ", 2/BL1, ", t_2 = ", 2/BL2);
print("");

analyse(nm, Fq) = { my(Ft, g, c, m2, lam, a, rec, cp, rts, kk, xi, badp);
  Ft = subst(Fq, 'q, qt);
  if(polcoeff(Ft,0) != 1, print("  ", nm, ": F(0) != 1, skipped"); return(0));
  g = sqrt(Ft); c = vector(NROW+1, i, polcoeff(g, i-1)); m2 = 0; badp = 1;
  for(n = 1, NROW, if(c[n+1], my(dn = denominator(c[n+1]), v2 = valuation(dn,2)); if(dn != 2^v2, badp = max(badp, dn/2^v2)); if(v2 > 0, m2 = max(m2, ceil(v2/n)))));
  if(badp > 1, print("  ", nm, ": denominators not 2-powers (odd part up to ", badp, ") -- not an R1 row"); return(0));
  lam = 2^m2; a = vector(NROW+1, i, lam^(i-1)*c[i]);
  if(vecmax(vector(NROW+1,i,denominator(a[i]))) != 1, print("  ", nm, ": not integral at lambda = ", lam); return(0));
  rec = fitminp(a, 6, 4);
  if(rec == 0, print("  ", nm, ": lambda = ", lam, ", NO recurrence with order<=6, degree<=4"); return(0));
  cp = charpol_from(rec[3], rec[1], rec[2]); rts = polroots(cp);
  print("  ", nm, ": lambda = ", lam, ", minimal recurrence order ", rec[1], ", degree ", rec[2]);
  print("      char roots: ", vector(#rts, i, if(abs(imag(rts[i]))<1e-30, real(rts[i]), rts[i])));
  print("      singular t = lambda/root: ", vector(#rts, i, if(abs(rts[i])>1e-30, lam/rts[i], 0)));
  if(rec[1] == 2, kk = compk(rec, NK); if(type(kk)=="t_VEC", xi = kk[2][NK+1]/a[NK+1]; print("      k' = ", kk[1], ",  xi' ~ ", xi)));
  1; };

print("=== (i) the base form F = E_{2,5} ===");
F5 = 1 + 6*sum(n=1,MQ-1, (sigma(n) - if(n%5==0,5*sigma(n/5),0))*'q^n) + O('q^MQ);
analyse("E_{2,5}", F5);
print("");
print("=== (ii) the gauge family F*R(t)^2, R in Q(t) -- Q(x)-DEPENDENT by construction ===");
PT = 1 - 88*tq - 64*tq^2 + O('q^MQ);
analyse("E_{2,5}*P(t)^2", F5*PT^2);
analyse("E_{2,5}*(1-t)^2 ", F5*(1-tq)^2);
print("");
print("=== (iii) level-5 nebentypus M_2(Gamma_0(5),chi_5), dim 2 ===");
GG = znstar(5,1); CH = znconreychar(GG,4); MFI = mfinit([5,2,[GG,CH]],4); BAS = mfbasis(MFI);
print("  basis size ", #BAS);
for(i = 1, #BAS, print("    basis[", i, "] coefs 0..8: ", mfcoefs(BAS[i], 8)));
CA = mfcoefs(BAS[2], MQ-1); CB = mfcoefs(BAS[1], MQ-1);
CA = vector(MQ, j, CA[j]/CA[1]);
print("  normalised E^{1,chi_5}: ", vector(9, j, CA[j]));
analyse("E^{1,chi_5}", sum(n=0,MQ-1, CA[n+1]*'q^n) + O('q^MQ));
print("");
print("=== (iv) one-parameter family  E^{1,chi_5} + s*E^{chi_5,1}, s in Q ===");
SLIST = [-6,-5,-4,-3,-2,-1,-1/2,-1/3,-1/4,-1/5,1/5,1/4,1/3,1/2,1,2,3,4,5,6,8,-8,10,-10,12,-12,16,-16,20,-20,24,-24,25,-25];
HITS = 0;
for(i = 1, #SLIST, my(sv = SLIST[i], Fq = sum(n=0,MQ-1, (CA[n+1]+sv*CB[n+1])*'q^n) + O('q^MQ)); HITS += analyse(concat("s = ", Str(sv)), Fq));
print("  candidates surviving: ", HITS);
