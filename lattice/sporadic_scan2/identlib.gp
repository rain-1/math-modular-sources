\\ Identify Apery limits against a period basis + cusp-form L-values.
\\ input: limits.txt lines "idx M w digits value"; output: ident.out
CACHE = Map();

constbasis(prec) = {
  default(realprecision, prec);
  [ [1, "1"],
    [zeta(2), "zeta(2)"],
    [zeta(3), "zeta(3)"],
    [lfun(-4, 2), "G=L(2,chi-4)"],
    [lfun(-3, 2), "L(2,chi-3)"],
    [lfun(-3, 3), "L(3,chi-3)"],
    [lfun(-4, 3), "L(3,chi-4)"],
    [lfun(5, 2),  "L(2,chi5)"],
    [lfun(5, 3),  "L(3,chi5)"],
    [lfun(-7, 2), "L(2,chi-7)"],
    [lfun(-8, 2), "L(2,chi-8)"],
    [lfun(8, 2),  "L(2,chi8)"],
    [lfun(12, 2), "L(2,chi12)"],
    [lfun(-3, 4), "L(4,chi-3)"],
    [Pi^3, "Pi^3"],
    [Pi^2*log(2), "Pi^2 log2"],
    [log(2)^3, "log(2)^3"] ];
}

cuspbasis(M, w, prec) = {
  my(key = Str(M, "_", w), res, k, LL, mf, EB, lf);
  if(mapisdefined(CACHE, key, &res), return(res));
  prec = 160;
  default(realprecision, prec + 20);
  res = List();
  k = w + 2;
  LL = List();
  fordiv(M, L, listput(LL, L));
  for(u = 2, 4, if(u*M <= 240, listput(LL, u*M)));
  LL = Set(Vec(LL));
  for(li = 1, #LL,
    my(L = LL[li]);
    if(L < 3, next);
    my(G = znstar(L,1), seen = vectorsmall(L), CH = List());
    for(m = 1, L,
      if(gcd(m,L) != 1 || seen[m], next);
      my(o = znorder(Mod(m,L)));
      for(j = 1, o, if(gcd(j,o)==1, seen[lift(Mod(m,L)^j)] = 1));
      if(chareval(G, znconreylog(G,m), -1) != if(k%2==0, 0, 1/2), next);
      if(eulerphi(o) > 4, next);
      listput(CH, m);
    );
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
              for(s = w, w+2,
                my(v = lfun(lf[u], s), nm = Str("L(f_", L, "w", k, "c", m, "n", e, "e", u, ",", s, ")"));
                listput(res, [real(v), nm]);
                if(abs(imag(v)) > 1e-25, listput(res, [imag(v), Str("Im", nm)]))
              )
            )
          )
        )
      , E, );
    );
  );
  res = Vec(res);
  mapput(CACHE, key, res);
  res;
}

trypair(x, B, prec) = {
  my(best = "", v, s, mx);
  for(i = 1, #B, for(j = i+1, #B,
    if(B[i][1] == 0 || B[j][1] == 0, next);
    v = lindep([1, x, B[i][1], B[j][1]]);
    if(#v != 4 || v[2] == 0, next);
    mx = vecmax(abs(v));
    if(mx > 3000, next);
    s = abs(v[1] + v[2]*x + v[3]*B[i][1] + v[4]*B[j][1]);
    if(s > 10.0^(-prec+10), next);
    best = Str(best, "  x = (", -v[1], " + ", -v[3], "*", B[i][2], " + ", -v[4], "*", B[j][2], ")/", v[2]);
    if(#best > 300, return(best));
  ));
  best;
}

tryrel(x, B, prec) = {
  my(best = "", v, s, mx);
  for(i = 1, #B,
    if(i > 1 && B[i][1] == 0, next);
    v = lindep([1, x, B[i][1]]);
    if(#v != 3 || v[2] == 0, next);
    mx = vecmax(abs(v));
    if(mx > 100000, next);
    s = abs(v[1] + v[2]*x + v[3]*B[i][1]);
    if(s > 10.0^(-prec+8), next);
    best = Str(best, "  x = (", -v[1], " + ", -v[3], "*", B[i][2], ")/", v[2]);
    if(#best > 400, break);
  );
  best;
}


trypair2(x, B1, B2, prec) = { my(best = "", v, s, mx);
  for(i = 1, #B1, for(j = 1, #B2,
    if(B1[i][1] == 0 || B2[j][1] == 0, next);
    v = lindep([1, x, B1[i][1], B2[j][1]]);
    if(#v != 4 || v[2] == 0, next);
    mx = vecmax(abs(v));
    if(mx > 3000, next);
    s = abs(v[1] + v[2]*x + v[3]*B1[i][1] + v[4]*B2[j][1]);
    if(s > 10.0^(-prec+10), next);
    best = Str(best, "  x = (", -v[1], " + ", -v[3], "*", B1[i][2], " + ", -v[4], "*", B2[j][2], ")/", v[2]);
    if(#best > 300, return(best));
  ));
  best;
}
