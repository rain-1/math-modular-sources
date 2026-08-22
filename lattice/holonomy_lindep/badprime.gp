\\ refine the "no k<=14" rows: split denominators into good/bad primes
default(parisizemax, 8G);
default(realprecision, 40);
read("/home/ubuntu/code/math-modular-sources/lattice/mum_survey/ops.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/multislope/sc_rows.gp");
mkcf(op) = { my(dz=op[3], Ps=op[4]); vector(dz+1, i, subst(Ps[i], 'X, 'n-(i-1))); };
findop(cl) = { for(t=1,#OPS, if(OPS[t][2]==cl, return(OPS[t]))); 0; };
chipoly(op) = { my(dz=op[3], Ps=op[4], dm=0);
  for(i=1,#Ps, dm=max(dm,poldegree(Ps[i],'X)));
  sum(i=0,dz, polcoeff(Ps[i+1],dm,'X)*'x^(dz-i)); };
NDEN = 200;
{
foreach(["4.4.34","4.8.23","4.4.36","4.4.39","4.5.30","4.8.26","4.5.22","4.5.46b","4.5.24","4.5.27","4.5.35","4.8.31","4.5.84"], cl,
  my(op = findop(cl)); if(op == 0, next);
  my(dz=op[3], cf=mkcf(op), cp=chipoly(op));
  my(rr=polroots(cp), m=vector(#rr,i,abs(rr[i])), idx=vecsort(m,,5));
  my(rat = m[idx[2]]/m[idx[1]]);
  \\ bad primes: divide the numerator lc of P_0 or any denominator in the P_i
  my(bad = Set());
  my(lc = polcoeff(op[4][1], poldegree(op[4][1],'X), 'X));
  my(den = 1); for(i=1,#op[4], den = lcm(den, denominator(content(op[4][i]))));
  foreach(factor(abs(numerator(lc)*den))[,1]~, p, bad = setunion(bad, Set([p])));
  my(res = vector(dz-1));
  for(j=1, dz-1,
    my(x = compan(cf, dz, j, NDEN));
    my(kg = 0, L = 1, brate = 0);
    for(n=1, NDEN, L = lcm(L,n);
      my(d = denominator(x[n+1]));
      \\ split
      my(db = 1); foreach(Vec(bad), p, my(v=valuation(d,p)); db *= p^v; d /= p^v);
      if(n == NDEN, brate = if(db>1, log(db)/n, 0));
      if(d > 1, my(kk=0, dd=d); while(dd>1, kk++; dd=dd/gcd(dd,L); if(kk>25,break));
        if(kk>kg, kg=kk));
    );
    res[j] = [kg, brate];
  );
  print(op[1], " | ", cl, " | dz=", dz, " | ratio=", rat, " | -log=", -log(rat),
        " | badprimes=", Vec(bad), " | [kgood, badrate] = ", res,
        " | genscore=", -log(rat) - vecmax(vector(dz-1,j,res[j][1])) - vecmax(vector(dz-1,j,res[j][2])));
);
}
quit;
