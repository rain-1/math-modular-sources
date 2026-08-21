default(parisizemax, 6000000000);
read("identlib.gp");
PREC = 100;

cuspall(prec) = { my(C = List(), LS);
  default(realprecision, prec+15);
  LS = [3,4,5,6,7,8,9,10,12,14,15,16,18,20,21,24,25,27,28,32,36,40,45,48,49,50,54,56,63,64,72];
  for(li = 1, #LS,
    my(L = LS[li]);
    for(k = 3, 3,
      my(G = znstar(L,1), seen = vectorsmall(L));
      for(m = 1, L,
        if(gcd(m,L)!=1 || seen[m], next);
        my(o = znorder(Mod(m,L)));
        for(j = 1, o, if(gcd(j,o)==1, seen[lift(Mod(m,L)^j)] = 1));
        if(chareval(G, znconreylog(G,m), -1) != if(k%2==0,0,1/2), next);
        if(eulerphi(o) > 4, next);
        iferr(
          my(mf = mfinit([L,k,[G,znconreylog(G,m)]], 0));
          if(mfdim(mf) > 0,
            my(EB = mfeigenbasis(mf));
            for(e = 1, #EB,
              my(lf = lfunmf(mf, EB[e]));
              lf = iferr(lfun(lf,2); [lf], E, lf);
              if(type(lf) != "t_VEC", lf = [lf]);
              for(u = 1, #lf, for(s = 1, k-1,
                my(v = lfun(lf[u], s), nn = Str("L(f",L,"w",k,"c",m,"n",e,"e",u,",",s,")"));
                listput(C, [real(v), nn]);
                if(abs(imag(v)) > 1e-30, listput(C, [imag(v), Str("Im",nn)]))))))
        , E, ))));
  Vec(C);
}

doone(nm, xs, C, B) = { my(x, msg, r1);
  default(realprecision, PREC + 20);
  x = eval(xs);
  r1 = lindep([1, x]);
  msg = "";
  if(vecmax(abs(r1)) < 10^12 && abs(r1[1]+r1[2]*x) < 10.0^(-PREC+8),
     msg = Str("RATIONAL ", -r1[1], "/", r1[2]),
     msg = tryrel(x, B, PREC));
  if(msg == "", msg = tryrel(x, concat([[1,"1"]], C), PREC));
  if(msg == "", msg = trypair(x, B, PREC));
  if(msg == "", msg = trypair2(x, C, B, PREC));
  if(msg == "", msg = Str("UNIDENTIFIED (", #B, " constants + ", #C, " cusp L-values at ", PREC, " digits)"));
  write("ident_sym1.out", Str(nm, " | ", msg));
  print(nm, " | ", msg);
  0;
}

{
my(LN = readstr("sym1_limits.txt"), C, B);
system("rm -f ident_sym1.out");
default(realprecision, PREC+20);
B = constbasis(PREC+15);
C = cuspall(PREC);
print("cusp basis size ", #C);
for(t = 1, #LN,
  my(f = strsplit(LN[t], " "));
  doone(f[1], f[2], C, B);
);
}
quit
