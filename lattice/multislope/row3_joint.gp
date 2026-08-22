default(parisizemax, 6000000000);
LOG = "/home/ubuntu/code/math-modular-sources/lattice/multislope/row3_joint.log";
W(s) = write(LOG, s);
read("/home/ubuntu/code/math-modular-sources/lattice/zeta5_two_row/level16_rows.txt");
read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row3_rec.txt");
NN = #An; QQ = QROW3; RR0 = 16;
W(Str("=== level-16 zeta(5): joint operator scan, terms to n=", NN-1, " ==="));
/* residual of the A-operator applied to B */
{ my(rs = vector(20, k, my(n=RR0+k-1); sum(i=0,RR0, subst(QQ[i+1],'n,n)*Bn[n-i+1])));
  W(Str("  residual R_n = L_A(B)_n for n=16..35 (should be 0 if B solves L_A):"));
  W(Str("    all zero? ", vecmax(vector(20,k,abs(rs[k])))==0));
  W(Str("    R_16 numerator size (digits) = ", #Str(numerator(rs[1])), ",  R_17 = ", #Str(numerator(rs[2])))); }
pp = 2^61-1;
Am = vector(NN, i, Mod(An[i],pp)); Bm = vector(NN, i, Mod(Bn[i],pp));
W("");
W("--- kernel dims mod 2^61-1 for an operator killing BOTH A and B (rows from A and from B) ---");
W("    r \\ D:  11   12   13   14   15");
{ for(r=16, 24,
    my(row = Str("    r=", r, " : "));
    for(dg=11, 15,
      my(nun=(r+1)*(dg+1), nr=(nun\2)+20, kd);
      if(r+nr <= NN-1,
        kd = #matker(matconcat([
          matrix(nr, nun, a, b, my(i=(b-1)\(dg+1), j=(b-1)%(dg+1), n=r+a); Am[n-i+1]*Mod(n,pp)^j);
          matrix(nr, nun, a, b, my(i=(b-1)\(dg+1), j=(b-1)%(dg+1), n=r+a); Bm[n-i+1]*Mod(n,pp)^j)])),
        kd = -1);
      row = Str(row, kd, "   "));
    W(row)); }
W("");
W("--- for comparison: kernel dims for B alone ---");
W("    r \\ D:  11   12   13   14   15");
{ for(r=16, 24,
    my(row = Str("    r=", r, " : "));
    for(dg=11, 15,
      my(nun=(r+1)*(dg+1), nr=nun+20, kd);
      if(r+nr <= NN-1,
        kd = #matker(matrix(nr, nun, a, b, my(i=(b-1)\(dg+1), j=(b-1)%(dg+1), n=r+a); Bm[n-i+1]*Mod(n,pp)^j)),
        kd = -1);
      row = Str(row, kd, "   "));
    W(row)); }
W("DONE"); quit;
