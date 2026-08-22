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
foreach(["4.5.1","4.5.2","4.5.4","4.5.5","4.5.6","4.5.7","4.4.33","4.5.10","4.4.34","4.4.36","4.5.15","4.5.17","4.5.18","4.5.20","4.5.22","4.5.23","4.5.24","4.8.16","4.8.18","4.5.25","4.5.27","4.5.28","4.5.29","4.5.30","4.4.39","4.5.31","4.8.20","4.5.32","4.5.34","4.5.35","4.5.36","4.5.37","4.4.42","4.5.39","4.5.40","4.5.42","4.4.44","4.5.43","4.5.44","4.5.49","4.5.51","4.5.52","4.5.54","4.5.55","4.5.56","4.5.66","4.4.58","4.8.22","4.8.23","4.8.24","4.5.70","4.4.62","4.5.74","4.5.75","4.8.26","4.8.27","4.8.28","4.5.76","4.5.77","4.5.78","4.5.79","4.5.80","4.5.81","4.8.30","4.8.31","4.5.83","4.8.32","4.5.84","4.5.86","4.5.87","4.8.33","4.8.34","4.5.88","4.8.35","4.8.36","4.4.68","4.5.95","4.5.97","4.5.98","4.8.38","4.5.104","4.5.107","4.6.22","4.5.110","4.5.113"], cl,
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
