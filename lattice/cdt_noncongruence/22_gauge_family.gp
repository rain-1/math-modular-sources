/* The only weight-2 forms F' on Gamma_0(5)+5 whose sqrt has NO extra singular
   point are F' = F*P(t)^b, P(t) = 1 - 44t - 16t^2 = prod(1 - t/t_i).
   b even -> g' = g*P^{b/2} is a Q(t)-multiple of g: Q(x)-DEPENDENT.
   b odd  -> g' = g*P^{(b-1)/2}*sqrt(P): genuinely new, but sqrt(P) doubles the
   characteristic roots (ROOT_ROWS hypothesis (H4)).  Tested here.
   Also: are the level-10/15/20/25 order-2 hits just E_{2,5} again?           */
default(parisize, 3000000000); default(realprecision, 50);
MQ = 70; NROW = 52; NK = 45; PGATE = 2147483647;
DL = vector(NK+2, i, if(i==1, 1, lcm(vector(i-1,j,j))));
etaprod(d, MM) = { my(s = 1 + O('q^MM)); forstep(n=d, MM-1, d, s *= (1 - 'q^n)); s; };
fitminp(A, maxo, maxd) = { my(L = #A - 1, Ap = vector(#A, i, Mod(A[i], PGATE))); for(o = 1, maxo, for(d = 0, maxd, my(nc = (o+1)*(d+1), rows = List(), r); if(L - o + 1 < nc + 6, next); for(n = o, L, r = vector(nc); for(j = 0, o, for(e = 0, d, r[j*(d+1)+e+1] = Mod(n,PGATE)^e * Ap[n+1-j])); listput(rows, r)); my(Mx = matconcat(Vec(rows)~), K = matker(Mx)); if(#K == 1, my(rows2 = List(), r2); for(n = o, L, r2 = vector(nc); for(j = 0, o, for(e = 0, d, r2[j*(d+1)+e+1] = n^e * A[n+1-j])); listput(rows2, r2)); my(M2 = matconcat(Vec(rows2)~), K2 = matker(M2)); if(#K2 == 1, return([o, d, K2[,1]]))))); 0; };
charpol_from(v, o, d) = { my(c = vector(o+1, j, v[(j-1)*(d+1)+d+1])); sum(j=1, o+1, c[j]*'L^(o+1-j)); };
compk(rec, NN) = { my(o=rec[1], d=rec[2], v=rec[3], b=vector(NN+1)); for(i=1,o, b[i]=0); b[o]=1; my(cf=vector(o+1, j, sum(e=0,d, v[(j-1)*(d+1)+e+1]*'n^e))); for(n=o, NN, my(s2=0); for(j=1,o, s2 += subst(cf[j+1],'n,n)*b[n+1-j]); my(lead=subst(cf[1],'n,n)); if(lead==0, return(-1)); b[n+1] = -s2/lead); my(k=0); while(k<=9, my(ok=1); for(n=0,NN, if(denominator(DL[n+1]^k*b[n+1])!=1, ok=0; break)); if(ok, break); k++); [k, b]; };
E1 = etaprod(1, MQ); E5 = etaprod(5, MQ);
uu = 'q*(E5/E1)^6 + O('q^MQ); tq = uu/(1+22*uu+125*uu^2) + O('q^MQ); qt = serreverse(tq);
F5 = 1 + 6*sum(n=1,MQ-1, (sigma(n) - if(n%5==0,5*sigma(n/5),0))*'q^n) + O('q^MQ);
PQ = 1 - 44*tq - 16*tq^2 + O('q^MQ);
print("P(t) = 1 - 44t - 16t^2, roots t_1,t_2 = ", polroots(1 - 44*'t - 16*'t^2));
analyse(nm, Fq) = { my(Ft, g, c, m2, lam, a, rec, cp, rts, badp, kk); Ft = subst(Fq, 'q, qt); if(polcoeff(Ft,0) != 1, print("    ", nm, ": F(0)!=1"); return(0)); g = sqrt(Ft); c = vector(NROW+1, i, polcoeff(g, i-1)); m2 = 0; badp = 1; for(n = 1, NROW, if(c[n+1], my(dn = denominator(c[n+1]), v2 = valuation(dn,2)); if(dn != 2^v2, badp = max(badp, dn/2^v2)); if(v2 > 0, m2 = max(m2, ceil(v2/n))))); if(badp > 1, print("    ", nm, ": denominators not 2-powers"); return(0)); lam = 2^m2; a = vector(NROW+1, i, lam^(i-1)*c[i]); if(vecmax(vector(NROW+1,i,denominator(a[i]))) != 1, print("    ", nm, ": not integral at lambda=", lam); return(0)); rec = fitminp(a, 6, 4); if(rec == 0, print("    ", nm, ": lambda=", lam, ", no recurrence order<=6 deg<=4"); return(0)); cp = charpol_from(rec[3], rec[1], rec[2]); rts = polroots(cp); print("    ", nm, ": lambda=", lam, ", order ", rec[1], ", deg ", rec[2]); print("        char roots ", vector(#rts,i,if(abs(imag(rts[i]))<1e-25, real(rts[i]), rts[i]))); kk = compk(rec, NK); if(type(kk)=="t_VEC", print("        k' = ", kk[1])); 1; };
print(""); print("=== F * P(t)^b ===");
analyse("b = 1 (F*P)", F5*PQ);
analyse("b = 2 (F*P^2)", F5*PQ^2);
analyse("b = 3 (F*P^3)", F5*PQ^3);
print(""); print("=== sqrt(P) alone (the new multiplier factor) ===");
{ my(sp = sqrt(subst(PQ,'q,qt)), c = vector(NROW+1,i,polcoeff(sp,i-1)), m2=0, badp=1); for(n=1,NROW, if(c[n+1], my(dn=denominator(c[n+1]), v2=valuation(dn,2)); if(dn!=2^v2, badp=max(badp,dn/2^v2)); if(v2>0, m2=max(m2,ceil(v2/n))))); print("    sqrt(P): lambda = ", 2^m2, ", odd denominators ", badp); my(lam=2^m2, a=vector(NROW+1,i,lam^(i-1)*c[i]), rec=fitminp(a,6,4)); if(rec, my(cp=charpol_from(rec[3],rec[1],rec[2])); print("    order ", rec[1], " deg ", rec[2], " roots ", polroots(cp)), print("    no recurrence")); }
print(""); print("=== are the level-10/15/20/25 order-2 hits just E_{2,5}? ===");
for(i = 2, 5, my(NN = 5*i, bas = mfbasis(mfinit([NN,2,1],4)), fnd = 0); for(j = 1, #bas, my(cf = mfcoefs(bas[j], 20)); if(cf[1] != 0, my(nc = vector(21, w, cf[w]/cf[1]), e5 = vector(21, w, polcoeff(F5, w-1))); if(nc == e5, print("  level ", NN, " basis[", j, "] IS E_{2,5}"); fnd = 1))); if(!fnd, print("  level ", NN, ": no basis element equals E_{2,5}; first coeffs of each normalised element:"); for(j=1,#bas, my(cf=mfcoefs(bas[j],8)); if(cf[1]!=0, print("    [",j,"] ", vector(9,w,cf[w]/cf[1]))))));
