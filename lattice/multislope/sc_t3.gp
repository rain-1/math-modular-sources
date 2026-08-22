\\ ---- T3: good-prime tower for companions --------------------------
\\ t_s = x_{a p^s}/a_{a p^s};  test whether v_p(t_{s+1}/t_s) is constant in s
\\ and whether the rescaled ratio p^{w} t_{s+1}/t_s tends to a limit.
t3row(nm, cf, r, N, primes, bases) = {
  my(av = genseq(cf, r, [1], N));
  print("");
  print("################ T3  ", nm, "   N=", N, " ################");
  for(j = 1, r-1,
    my(x = compan(cf, r, j, N));
    print("");
    print("--- companion X^(", j, ") ---");
    for(pi = 1, #primes,
      my(p = primes[pi]);
      for(bi = 1, #bases,
        my(aa = bases[bi]);
        my(ts = List(), ss = 0, ps = 1);
        while(aa*ps <= N, listput(ts, x[aa*ps+1]/av[aa*ps+1]); ps *= p);
        if(#ts < 3, next);
        my(vs = vector(#ts, i, if(ts[i]==0, "INF", valuation(ts[i], p))));
        my(dv = List());
        for(i = 1, #ts-1, if(ts[i]==0 || ts[i+1]==0, listput(dv,"?"),
                            listput(dv, valuation(ts[i+1],p) - valuation(ts[i],p))));
        print("  p=", p, " a=", aa, "  v_p(t_s), s=0..: ", vs);
        print("           v_p(t_{s+1}/t_s): ", Vec(dv));
        \\ if the increments are eventually constant = -w, look at the unit part mod p^3
        my(dd = Vec(dv), ok = 1, w = 0);
        if(#dd >= 2, w = dd[#dd]; for(i = max(1,#dd-2), #dd, if(dd[i] != w, ok = 0)));
        if(ok && type(w)=="t_INT",
          my(un = List());
          for(i = 1, #ts-1, if(ts[i]!=0 && ts[i+1]!=0,
             my(rt = ts[i+1]/ts[i]);
             rt = rt / p^valuation(rt, p);
             listput(un, lift(Mod(numerator(rt),p^10)*Mod(denominator(rt),p^10)^(-1)))));
          my(uu = Vec(un));
          for(e = 1, 6, print("           u_s = p^",-w," t_{s+1}/t_s  mod p^",e,": ", vector(#uu,i,uu[i]%p^e)));
        );
      );
    );
  );
};
