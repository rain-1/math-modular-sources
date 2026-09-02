default(parisize,"3G");
read("lib.gp");
read("src.gp");
N = 3000;
{ for(i=1,#ROWS,
  my(R=ROWS[i], r=R[2], AB=genrow(R,N), B=AB[2], dn=1, dr=1, Rn=1, kmax=0);
  my(eqR=1, firstbad=0, badset=List(), kk, dd, dv, hasC=(i<=12), cf);
  if(hasC, cf=CMS[i]);
  for(n=1,N,
    dn = lcm(dn,n); dr = dn^r;
    if(hasC, my(t=n^r, g=gcd(t, cf(n))); Rn = lcm(Rn, t/g));
    dv = denominator(B[n+1]);
    kk = 0; while((dn^kk) % dv != 0 && kk<12, kk++); if(kk>kmax, kmax=kk);
    if(hasC && dv != Rn, eqR=0; if(firstbad==0, firstbad=n); my(q=Rn/dv); for(j=1,#factor(q)[,1]~, listput(badset, factor(q)[j,1])));
  );
  print(R[1], " | r=", r, " | min k with d_n^k b_n in Z (n<=", N, "): ", kmax, if(hasC, Str(" | den(b_n)=R_n for all n<=",N,": ", if(eqR,"YES",Str("NO, first n=",firstbad," primes ",Set(Vec(badset))))), " | no Eisenstein source")));
}
quit;
