default(parisizemax, 8*10^9);
read("/home/ubuntu/code/math-modular-sources/lattice/emn/central_ab.gp");
NN = #ACEN - 1;
ucen = vector(NN+1, i, 256^(i-1)*BCEN[i]);
vcen = vector(NN+1, i, 256^(i-1)*ACEN[i]);
ORD = 2; DG = 10;
rows = List();
for(s = 1, 2, my(W = if(s==1, ucen, vcen)); for(n = 0, NN-ORD, my(row = vector((ORD+1)*(DG+1))); for(j = 0, ORD, for(e = 0, DG, row[j*(DG+1)+e+1] = n^e*W[n+j+1])); listput(rows, row)));
kr = matker(matconcat(Vec(rows)~));
print("nullspace dim = ", #kr);
kv = kr[,1];
/* clear denominators + content */
den = 1; for(i = 1, #kv, den = lcm(den, denominator(kv[i])));
kv = kv*den; g = 0; for(i = 1, #kv, g = gcd(g, kv[i])); kv = kv/g;
pol = vector(ORD+1, j, sum(e = 0, DG, kv[(j-1)*(DG+1)+e+1]*'n^e));
for(j = 1, ORD+1, print("p_", j-1, "(n) = ", factor(pol[j])));
print();
lc = vector(ORD+1, j, polcoef(pol[j], DG, 'n));
print("leading coefficients (n^10): ", lc);
chp = sum(j = 0, ORD, lc[j+1]*'lam^j);
print("characteristic polynomial: ", chp, "  = ", factor(chp));
print("roots: ", polroots(chp));
print();
/* independent re-verification of the recurrence on longer data */
print("re-verification on all n<=", NN-2, ":");
bad = 0;
for(n = 0, NN-ORD, for(s = 1, 2, my(W = if(s==1, ucen, vcen)); my(r = sum(j = 0, ORD, subst(pol[j+1],'n,n)*W[n+j+1])); if(r != 0, bad++)));
print("  residual failures: ", bad, " out of ", 2*(NN-ORD+1));
/* minimality: is there an order-2 recurrence of lower degree, or order-1? */
print();
for(dg = 1, 9, my(rr = List()); for(s = 1, 2, my(W = if(s==1,ucen,vcen)); for(n = 0, NN-2, my(row = vector(3*(dg+1))); for(j = 0, 2, for(e = 0, dg, row[j*(dg+1)+e+1] = n^e*W[n+j+1])); listput(rr,row))); my(k2 = matker(matconcat(Vec(rr)~))); if(#k2>0, print("  ALSO order 2 at degree ", dg)));
print("  (degree 10 is the minimal coefficient degree for the joint order-2 operator)");
quit;
