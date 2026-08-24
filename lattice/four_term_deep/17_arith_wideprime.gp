default(parisizemax, 12000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/17_arith_lib.gp");
fn = "/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/arith_wideprime.log";
logit(s) = { print(s); write(fn, s); };
NN = 300;
PL = primes(20);   \\ 2..71
wide(nm, av, bv) = {
  my(out = "", cas = vector(NN+1));
  for(n=0,NN-1, cas[n+1] = av[n+1]*bv[n+2] - av[n+2]*bv[n+1]);
  for(i=1,#PL, my(p=PL[i], v1, v2, w1, w2, k1, k2);
    v1 = my(d=bv[151]/av[151]-bv[150]/av[150]); v1 = if(d==0,0,valuation(bv[151]/av[151]-bv[150]/av[150],p));
    v2 = if(bv[NN+1]/av[NN+1]-bv[NN]/av[NN]==0,0,valuation(bv[NN+1]/av[NN+1]-bv[NN]/av[NN],p));
    w1 = valuation(cas[151],p); w2 = valuation(cas[NN],p);
    k1 = valuation(denominator(bv[151]),p); k2 = valuation(denominator(bv[NN+1]),p);
    out = Str(out, "  p=",p," sig=",(v2-v1)*1.0/150," w=",(w2-w1)*1.0/149," kap=",-(k2-k1)*1.0/150,"\n"));
  logit(Str("=== ", nm, " ===\n", out));
};
{for(i=1,#ROWS, my(rr = row4(ROWS[i][2], NN)); wide(ROWS[i][1], rr[1], rr[2]));}
logit("DONE");
quit;
