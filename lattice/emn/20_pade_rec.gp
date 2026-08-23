/* (h)/(2.1) joint holonomic system for the radial Pade pair.
   NOTE: with data only to n=34 the joint fit is under-determined; the decisive
   computation is 21_pade_u.gp (the b-side alone, to n=160).  Kept for the record. */
default(parisizemax, 12*10^9);
read("/home/ubuntu/code/math-modular-sources/lattice/emn/pade_ab.gp");
NN = #APAD - 1;
upad = vector(NN+1, i, 256^(i-1)*BPAD[i]);
vpad = vector(NN+1, i, 256^(i-1)*APAD[i]);
print("n <= ", NN, ";  256^n b_n integral? ", vecmax(vector(NN+1,i,denominator(upad[i])))==1);
print("den(256^n a_n) odd? ", vecmax(vector(NN+1,i,denominator(vpad[i])%2))==1);
fitr(seqs, ord, dg, nmax) = my(rows = List(), unk = (ord+1)*(dg+1)); for(s = 1, #seqs, my(W = seqs[s]); for(n = 0, nmax-ord, my(row = vector(unk)); for(j = 0, ord, for(e = 0, dg, row[j*(dg+1)+e+1] = n^e*W[n+j+1])); listput(rows,row))); matker(matconcat(Vec(rows)~));
print();
print("=== joint (U,V): only well-determined fits are reported (unknowns + 4 <= #equations) ===");
for(ord = 2, 4, for(dg = 1, 18, if((ord+1)*(dg+1) + 4 <= 2*(NN-ord+1), my(k = fitr([upad,vpad], ord, dg, NN)); if(#k == 1, print("  order ", ord, " degree ", dg, ": nullspace dim 1")))));
print("  (no well-determined joint operator of order <= 4 within reach of n <= 34;");
print("   see 21_pade_u.gp for the b-side operator, order 3, degree 16.)");
quit;
