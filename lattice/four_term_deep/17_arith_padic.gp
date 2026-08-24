default(parisizemax, 20000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/17_arith_lib.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/mum_survey/lpgen.gp");
fn = "/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/arith_padic.log";
logit(s) = { print(s); write(fn, s); };
NN = 400;
PRcap = 300;
PLIST = [2,3,5,7,11,13];

{QCH = [["1", triv]];
 for(D=-24,24, if(D!=0 && D!=1 && isfundamental(D),
   QCH = concat(QCH, [[Str("chi",D), [abs(D), vector(abs(D), a, kronecker(D,a))]]])));}

TGCACHE = List();
targets(p, PR) = {
  my(kk = Str(p,"_",PR));
  for(i=1,#TGCACHE, if(TGCACHE[i][1]==kk, return(TGCACHE[i][2])));
  my(tt = vector(#QCH, u, [Str("L_",p,"(2,",QCH[u][1],")"), LpG(p, QCH[u][2], -1, 2, PR)]));
  listput(TGCACHE, [kk, tt]); tt;
};

incval(av, bv, p, n) = {
  my(d = bv[n+1]/av[n+1] - bv[n]/av[n]);
  if(d == 0, -1, valuation(d, p));
};

anal(nm, av, bv) = {
  logit(Str("=== ", nm, " ==="));
  for(t = 1, #PLIST,
    my(p = PLIST[t], n1 = NN\2, n2 = NN, s1, s2, sig);
    s1 = incval(av,bv,p,n1); s2 = incval(av,bv,p,n2);
    sig = (s2-s1)*1.0/(n2-n1);
    logit(Str("  p=", p, "  sigma~", sig, "   v_p(incr_N)=", s2));
    if(sig <= 0.2 || s2 < 40, logit("     -> no p-adic Apery limit"),
      my(PR = min(PRcap, s2-8), xi, TG, nh = 0);
      xi = bv[NN+1]/av[NN+1] + O(p^PR);
      TG = targets(p, PR);
      logit(Str("     xi_", p, " (", PR, " digits) leading: ", lift(xi + O(p^16))));
      for(u = 1, #TG,
        if(TG[u][2] == 0, next);
        my(rel = lindep([xi, TG[u][2]]));
        if(#rel == 2 && rel[1] != 0 && abs(rel[1]) < 10^6 && abs(rel[2]) < 10^6,
          nh++;
          logit(Str("     HIT  xi_",p," = ", -rel[2],"/",rel[1], " * ", TG[u][1],
                    "   [resid v_p = ",
                    my(r = xi + rel[2]/rel[1]*TG[u][2]); if(r==0, "exact", valuation(r,p)),
                    " / ", PR, "]"))));
      my(rr = lindep([xi, 1 + O(p^PR)]));
      if(#rr==2 && rr[1]!=0 && abs(rr[1])<10^6 && abs(rr[2])<10^6,
         nh++; logit(Str("     HIT  xi_",p," = rational ", -rr[2],"/",rr[1])));
      if(nh == 0, logit("     NO HIT in battery (height < 10^6)"))));
  logit("");
};

{for(i=1,#ROWS, my(rr = row4(ROWS[i][2], NN)); anal(ROWS[i][1], rr[1], rr[2]));}
{my(rr = zrow2(12,4,32,NN)); anal("Zagier E (12,4,32)", rr[1], rr[2]);}
logit("DONE");
quit;
