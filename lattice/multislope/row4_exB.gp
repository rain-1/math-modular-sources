default(parisizemax, 6000000000);
An = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_An.txt");
NMAX = #An - 1;

kervec(r, D, pp) = {
  my(nc = (r+1)*(D+1));
  my(nlo = r);
  my(nhi = min(NMAX, nlo+nc+30));
  my(nr = nhi-nlo+1);
  my(mat = matrix(nr, nc, a, b, Mod(An[nlo+a-1-((b-1)\(D+1))+1],pp) * Mod(nlo+a-1,pp)^((b-1)%(D+1))));
  my(kk = matker(mat));
  if(matsize(kk)[2] != 1, return(0));
  my(v = kk[,1]);
  my(f = 0);
  for(t=1,nc, if(v[t]!=0, f=t; break));
  [f, v/v[f]]
}

reconstruct(r, D, nprimes) = {
  my(nc = (r+1)*(D+1));
  my(pl = List());
  my(vl = List());
  my(f0 = 0);
  my(pr = 2^61-1);
  for(t = 1, nprimes,
    pr = precprime(pr-1);
    my(res = kervec(r, D, pr));
    if(res == 0, print("  kernel dim != 1 at prime ", pr); return(0));
    if(f0 == 0, f0 = res[1], if(res[1] != f0, print("  pivot mismatch"); return(0)));
    listput(pl, pr);
    listput(vl, res[2])
  );
  my(cr = vector(nc));
  for(b = 1, nc,
    my(cm = chinese(vector(nprimes, t, Mod(lift(vl[t][b]), pl[t]))));
    cr[b] = bestappr(cm)
  );
  [f0, cr]
}

verify(r, D, cv) = {
  my(dd = 1);
  for(b=1,#cv, dd = lcm(dd, denominator(cv[b])));
  my(ci = vector(#cv, b, cv[b]*dd));
  my(bad = 0);
  for(n = r, NMAX,
    my(sm = 0);
    for(i=0,r,
      my(a = An[n-i+1]);
      my(nj = 1);
      for(j=0,D, sm += ci[i*(D+1)+j+1]*a*nj; nj*=n)
    );
    if(sm != 0, bad++; if(bad<3, print("   FAIL at n=",n)))
  );
  print("  verify over n=",r,"..",NMAX,": failures = ", bad);
  ci
}
print("=== shape B: r=35, D=7 (minimal ODE order) ===");
resB = reconstruct(35, 7, 6);
print("  pivot index = ", resB[1]);
ciB = verify(35, 7, resB[2]);
write("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_recB.txt", ciB);
print("  max coeff digits = ", vecmax(vector(#ciB,b,#Str(abs(ciB[b])))));
quit
