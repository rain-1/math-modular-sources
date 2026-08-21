\\ Targeted identification for a few hard rows: wide cusp-form basis, s=1,2,3,
\\ weights w+1, w+2, w+3, levels = all multiples/divisors up to LMAX.
read("identlib.gp");
LMAX = 160;

widebasis(w, prec) = {
  my(res = List(), k, G, seen, CH, mf, EB, lf);
  default(realprecision, prec+20);
  for(L = 3, LMAX,
    for(k = w+1, w+3,
      G = znstar(L,1); seen = vectorsmall(L); CH = List();
      for(m = 1, L,
        if(gcd(m,L) != 1 || seen[m], next);
        my(o = znorder(Mod(m,L)));
        for(j = 1, o, if(gcd(j,o)==1, seen[lift(Mod(m,L)^j)] = 1));
        if(chareval(G, znconreylog(G,m), -1) != if(k%2==0, 0, 1/2), next);
        if(eulerphi(o) > 2, next);
        listput(CH, m));
      for(ci = 1, #CH,
        my(m = CH[ci]);
        iferr(
          mf = mfinit([L, k, [G, znconreylog(G,m)]], 0);
          if(mfdim(mf) > 0,
            EB = mfeigenbasis(mf);
            for(e = 1, #EB,
              lf = lfunmf(mf, EB[e]);
              lf = iferr(lfun(lf, 2); [lf], E, lf);
              if(type(lf) != "t_VEC", lf = [lf]);
              for(u = 1, #lf,
                for(s = 1, k-1,
                  my(v = lfun(lf[u], s), nm = Str("L(f_", L, "w", k, "c", m, "n", e, "e", u, ",", s, ")"));
                  listput(res, [real(v), nm]);
                  if(abs(imag(v)) > 1e-25, listput(res, [imag(v), Str("Im", nm)]))))))
        , E, ))));
  Vec(res);
}

{
my(L = readstr("limits_hard.txt"), out = "ident_deep.out");
system("rm -f ident_deep.out");
my(WB = Map());
for(i = 1, #L,
  my(f = strsplit(L[i], " "));
  my(idx = eval(f[1]), M = eval(f[2]), w = eval(f[3]), dg = eval(f[4]));
  my(prec = min(dg, 120));
  default(realprecision, prec + 15);
  my(x = eval(f[5]));
  my(B = constbasis(prec + 15), W, key = Str(w));
  if(!mapisdefined(WB, key, &W), W = widebasis(w, 140); mapput(WB, key, W));
  default(realprecision, prec + 15);
  my(msg = tryrel(x, W, prec));
  if(msg == "", msg = trypair2(x, W, B, prec));
  write(out, Str(idx, " | ", if(msg == "", "UNIDENTIFIED(wide)", msg)));
  print(idx, " done");
);
}
quit
