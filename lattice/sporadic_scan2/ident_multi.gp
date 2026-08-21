\\ Multi-term identification: for each hard row, lindep the limit against
\\ 1 together with ALL newform L-values L(f, w+1) (and L(f,w), L(f,w+2)) of weight
\\ w+2 at every level dividing 4M (real characters and quartic orbits).
read("identlib.gp");

allL(M, w, prec) = {
  my(res = List(), k = w+2, G, seen, CH, mf, EB, lf, LL);
  default(realprecision, prec+20);
  LL = List(); fordiv(4*M, L, listput(LL, L)); LL = Set(Vec(LL));
  for(li = 1, #LL,
    my(L = LL[li]); if(L < 3, next);
    G = znstar(L,1); seen = vectorsmall(L); CH = List();
    for(m = 1, L,
      if(gcd(m,L) != 1 || seen[m], next);
      my(o = znorder(Mod(m,L)));
      for(j = 1, o, if(gcd(j,o)==1, seen[lift(Mod(m,L)^j)] = 1));
      if(chareval(G, znconreylog(G,m), -1) != if(k%2==0, 0, 1/2), next);
      if(eulerphi(o) > 4, next);
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
              my(v = lfun(lf[u], w+1), nm = Str("L(f_", L, "w", k, "c", m, "n", e, "e", u, ",", w+1, ")"));
              listput(res, [real(v), nm]);
              if(abs(imag(v)) > 1e-25, listput(res, [imag(v), Str("Im", nm)])))))
      , E, )));
  Vec(res);
}

{
my(L = readstr("limits_hard.txt"), out = "ident_multi.out");
system("rm -f ident_multi.out");
for(i = 1, #L,
  my(f = strsplit(L[i], " "));
  my(idx = eval(f[1]), M = eval(f[2]), w = eval(f[3]), dg = eval(f[4]));
  my(prec = min(dg, 160));
  default(realprecision, prec + 20);
  my(x = eval(f[5]));
  my(A = iferr(allL(M, w, prec), E, []));
  default(realprecision, prec + 20);
  my(msg = "UNIDENTIFIED(multi)");
  if(#A > 0,
    my(vec = concat([1, x], vector(#A, j, A[j][1])));
    my(v = lindep(vec));
    my(mx = vecmax(abs(v)));
    my(s = abs(sum(j=1, #vec, v[j]*vec[j])));
    if(v[2] != 0 && mx < 10^6 && s < 10.0^(-prec+15),
      msg = Str("x = (", -v[1]);
      for(j = 3, #v, if(v[j] != 0, msg = Str(msg, " + ", -v[j], "*", A[j-2][2])));
      msg = Str(msg, ")/", v[2]);
    );
  );
  write(out, Str(idx, " | #basis=", #A, " | ", msg));
  print(idx, " ", msg);
);
}
quit
