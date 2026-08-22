\\ ---- T3b: sharpened tower-limit test -------------------------------
agdepth(u, v, p, cap) = { my(e = 0); while(e < cap && (u-v) % p^(e+1) == 0, e++); e; };

t3b(nm, cf, r, N, primes, bases) = {
  my(av = genseq(cf, r, [1], N));
  print("");
  print("######## T3b ", nm, "  N=", N, " ########");
  for(j = 1, r-1,
    my(xv = compan(cf, r, j, N));
    for(pi = 1, #primes,
      my(p = primes[pi]);
      for(bi = 1, #bases,
        my(aa = bases[bi]);
        my(ts = List(), ps = 1);
        while(aa*ps <= N, listput(ts, xv[aa*ps+1]/av[aa*ps+1]); ps *= p);
        ts = Vec(ts);
        \\ drop the leading zeros of the companion (x_{a p^s}=0 for small s)
        my(i0 = 1); while(i0 <= #ts && ts[i0] == 0, i0++);
        ts = vector(#ts-i0+1, i, ts[i0+i-1]);
        if(#ts < 3, next);
        my(dv = vector(#ts-1, i, if(ts[i]==0||ts[i+1]==0, 999, valuation(ts[i+1],p)-valuation(ts[i],p))));
        my(w = dv[#dv], cst = 1);
        for(i = max(1,#dv-1), #dv, if(dv[i] != w, cst = 0));
        if(#dv >= 3, for(i = #dv-2, #dv, if(dv[i] != w, cst = 0)));
        my(msg = "  [valuation NOT eventually constant]");
        if(cst && w != 999,
          my(us = List());
          for(i = 1, #ts-1, if(ts[i] && ts[i+1],
             my(rt = ts[i+1]/ts[i]); rt = rt/p^valuation(rt,p);
             listput(us, lift(Mod(numerator(rt),p^40)*Mod(denominator(rt),p^40)^(-1)))));
          us = Vec(us);
          my(dep = vector(max(0,#us-1), i, agdepth(us[i], us[i+1], p, 40)));
          my(dl = if(#dep, dep[#dep], 0));
          my(uu = us[#us] % p^dl);
          my(tg = if(dl==0, "", if(uu==1%p^dl, "  (= +1)", if(uu==lift(Mod(-1,p^dl)), "  (= -1)", ""))));
          msg = Str("  w=", -w, "  consec-agree-depths ", dep, "   u_last mod p^", dl, " = ", uu, tg);
        );
        print("  X^(",j,") p=",p," a=",aa," | dv=", dv, msg);
      );
    );
  );
};
