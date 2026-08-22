/* Candidates at levels 10,15,20,25 on the SAME t-line (t is the Gamma_0(5)+5
   Hauptmodul).  For every rational basis element of M_2(Gamma_0(N)) with
   a_0 != 0, normalise to F'(0)=1 and run the same analysis. */
default(parisize, 3000000000);
default(realprecision, 50);
MQ = 70; NROW = 52; NK = 45; PGATE = 2147483647;
DL = vector(NK+2, i, if(i==1, 1, lcm(vector(i-1,j,j))));
etaprod(d, MM) = { my(s = 1 + O('q^MM)); forstep(n=d, MM-1, d, s *= (1 - 'q^n)); s; };
fitminp(A, maxo, maxd) = { my(L = #A - 1, Ap = vector(#A, i, Mod(A[i], PGATE))); for(o = 1, maxo, for(d = 0, maxd, my(nc = (o+1)*(d+1), rows = List(), r); if(L - o + 1 < nc + 6, next); for(n = o, L, r = vector(nc); for(j = 0, o, for(e = 0, d, r[j*(d+1)+e+1] = Mod(n,PGATE)^e * Ap[n+1-j])); listput(rows, r)); my(Mx = matconcat(Vec(rows)~), K = matker(Mx)); if(#K == 1, my(rows2 = List(), r2); for(n = o, L, r2 = vector(nc); for(j = 0, o, for(e = 0, d, r2[j*(d+1)+e+1] = n^e * A[n+1-j])); listput(rows2, r2)); my(M2 = matconcat(Vec(rows2)~), K2 = matker(M2)); if(#K2 == 1, return([o, d, K2[,1]]))))); 0; };
charpol_from(v, o, d) = { my(c = vector(o+1, j, v[(j-1)*(d+1)+d+1])); sum(j=1, o+1, c[j]*'L^(o+1-j)); };
E1 = etaprod(1, MQ); E5 = etaprod(5, MQ);
uu = 'q*(E5/E1)^6 + O('q^MQ); tq = uu/(1+22*uu+125*uu^2) + O('q^MQ); qt = serreverse(tq);
BT1 = 0.02254248593736856025573354295704764715; BT2 = -2.77254248593736856025573354295704764715;
analyse(nm, Fq) = { my(Ft, g, c, m2, lam, a, rec, cp, rts, badp, sv, extra);
  Ft = subst(Fq, 'q, qt);
  if(polcoeff(Ft,0) != 1, return(0));
  g = sqrt(Ft); c = vector(NROW+1, i, polcoeff(g, i-1)); m2 = 0; badp = 1;
  for(n = 1, NROW, if(c[n+1], my(dn = denominator(c[n+1]), v2 = valuation(dn,2)); if(dn != 2^v2, badp = max(badp, dn/2^v2)); if(v2 > 0, m2 = max(m2, ceil(v2/n)))));
  if(badp > 1, print("    ", nm, ": denominators not 2-powers -- not an R1 row"); return(0));
  lam = 2^m2; a = vector(NROW+1, i, lam^(i-1)*c[i]);
  if(vecmax(vector(NROW+1,i,denominator(a[i]))) != 1, print("    ", nm, ": not integral at lambda = ", lam); return(0));
  rec = fitminp(a, 6, 4);
  if(rec == 0, print("    ", nm, ": lambda = ", lam, ", no recurrence order<=6 degree<=4"); return(0));
  cp = charpol_from(rec[3], rec[1], rec[2]); rts = polroots(cp);
  sv = vector(#rts, i, if(abs(rts[i])>1e-30, real(lam/rts[i]), 0));
  extra = 0; for(i=1,#sv, if(abs(sv[i]-BT1)>1e-12 && abs(sv[i]-BT2)>1e-12, extra=1));
  print("    ", nm, ": lambda = ", lam, ", order ", rec[1], ", deg ", rec[2], ", singular t = ", sv, if(extra, "   <-- EXTRA singular points", "   <-- same singular set"));
  if(!extra && rec[1]==2, 1, 0); };
for(i = 2, 5, my(NN = 5*i, mfi, bas); print("level ", NN, ":"); mfi = mfinit([NN,2,1],4); bas = mfbasis(mfi); print("  dim M_2 = ", #bas); for(j = 1, #bas, my(cf = mfcoefs(bas[j], MQ-1)); if(cf[1] != 0, analyse(concat(concat("N=",Str(NN)), concat(" basis[",concat(Str(j),"]"))), sum(n=0,MQ-1, (cf[n+1]/cf[1])*'q^n) + O('q^MQ)), print("    basis[", j, "]: a_0 = 0, cannot normalise to F(0)=1"))));
