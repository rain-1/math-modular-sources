default(parisizemax, 20000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/17_arith_lib.gp");
fn = "/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/arith_two.log";
logit(s) = { print(s); write(fn, s); };
NN = 802;
ns = [100,200,300,400,500,600,700,800];
vv(x,p) = if(x==0, "ZERO", valuation(x,p));

VA = List(); VC = List(); VI = List(); NMS = List();
rep2(nm, av, bv) = {
  my(cas = vector(NN+1), va, vc, vi, z=0);
  for(n=0,NN-1, cas[n+1] = av[n+1]*bv[n+2] - av[n+2]*bv[n+1];
                if(cas[n+1]==0, z=n));
  va = vector(#ns, j, valuation(av[ns[j]+1],2));
  vc = vector(#ns, j, valuation(cas[ns[j]+1],2));
  vi = vector(#ns, j, my(n=ns[j], d=bv[n+1]/av[n+1]-bv[n]/av[n]); if(d==0,-10^9,valuation(d,2)));
  logit(Str("--- ", nm, "  (Cas zero at n=", z, ") ---"));
  logit(Str("  n           : ", ns));
  logit(Str("  v2(a_n)     : ", va, "   rate ", (va[#ns]-va[#ns-1])*1.0/100));
  logit(Str("  v2(Cas_n)   : ", vc, "   w2   ", (vc[#ns]-vc[#ns-1])*1.0/100, "   (vc/n: ", vector(#ns,j,vc[j]*1.0/ns[j]), ")"));
  logit(Str("  v2(incr_n)  : ", vi, "   sig2 ", (vi[#ns]-vi[#ns-1])*1.0/100, "   (vi/n: ", vector(#ns,j,vi[j]*1.0/ns[j]), ")"));
  logit(Str("  v2(a_n) n=40..60 : ", vector(21,i,valuation(av[40+i],2))));
  listput(VA, va); listput(VC, vc); listput(VI, vi); listput(NMS, nm);
  logit("");
};

{for(i=1,#ROWS, my(rr = row4(ROWS[i][2], NN)); rep2(ROWS[i][1], rr[1], rr[2]));}
{my(rr = zrow2(12,4,32,NN)); rep2("Zagier E (12,4,32)", rr[1], rr[2]);}
logit("DONE");
quit;
