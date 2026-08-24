default(parisizemax, 8000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/17_arith_lib.gp");
fn = "/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/arith_main.log";
logit(s) = { print(s); write(fn, s); };
NN = 400;
ns = [50, 100, 200, 300, 399];
PL = [2,3,5,7,11,13];

/* d_n = lcm(1..n) */
{DN = vector(NN+1); DN[1] = 1; for(n=1,NN, DN[n+1] = lcm(DN[n], n));}

slopes(vv, ns) = (vv[#ns]-vv[#ns-1])*1.0/(ns[#ns]-ns[#ns-1]);

report(nm, av, bv) = {
  my(allint = 1, kok = 1, kmax = 0, cas = vector(NN+1));
  logit(Str("=== ", nm, "  (N=", NN, ") ==="));
  for(n=0,NN, if(denominator(av[n+1])!=1, allint=0));
  for(n=0,NN, my(d=denominator(bv[n+1])); if(DN[n+1]^2 % d != 0, kok=0);
              if(d>1, my(e=0,dd=d,g); while(dd>1, g=gcd(dd,DN[n+1]); dd/=g; e++);
                      if(e>kmax, kmax=e)));
  logit(Str("  a_n integral for all n<=",NN,"? ", allint,
            "   den(b_n) | d_n^2 for all n<=",NN,"? ", kok, "   sharp k = ", kmax));
  /* 2x2 Casoratian */
  for(n=0,NN-1, cas[n+1] = av[n+1]*bv[n+2] - av[n+2]*bv[n+1]);
  logit(Str("  Cas_n = a_n b_{n+1} - a_{n+1} b_n,  n=0..6: ", vector(7,i,cas[i])));
  logit(Str("  Cas_n*(n+1)^2, n=0..8: ", vector(9,i,cas[i]*i^2)));
  logit(Str("  ratio Cas_n/Cas_{n-1} * (n+1)^2/n^2, n=1..8: ",
            vector(8,i, if(cas[i]==0, "-", cas[i+1]/cas[i]*(i+1)^2/i^2))));
  for(i=1,#PL,
    my(p=PL[i], kv=vector(#ns), wv=vector(#ns), sv=vector(#ns), av_=vector(#ns));
    for(j=1,#ns, my(n=ns[j]);
      kv[j] = valuation(denominator(bv[n+1]), p);
      wv[j] = if(cas[n+1]==0, 10^9, valuation(cas[n+1], p));
      av_[j] = valuation(av[n+1], p);
      my(d = bv[n+1]/av[n+1] - bv[n]/av[n]);
      sv[j] = if(d==0, 10^9, valuation(d, p)));
    logit(Str("  p=",p," v_p(den b_n) : ", kv, "  rate ", slopes(kv,ns),
              "   | v_p(a_n): ", av_));
    logit(Str("  p=",p," v_p(Cas_n)   : ", wv, "  w_p  ", slopes(wv,ns)));
    logit(Str("  p=",p," v_p(incr)    : ", sv, "  sig  ", slopes(sv,ns))));
  logit("");
};

{for(i=1,#ROWS, my(rr = row4(ROWS[i][2], NN)); report(ROWS[i][1], rr[1], rr[2]));}
{my(rr = zrow2(12,4,32,NN)); report("Zagier E (12,4,32)", rr[1], rr[2]);}
logit("DONE");
quit;
