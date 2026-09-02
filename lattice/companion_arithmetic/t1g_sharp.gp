default(parisize,"4G");
read("lib.gp");
read("src.gp");
N = 3000;
{ for(i=1,#ROWS,
  my(R=ROWS[i], r=R[2], k=if(i<=12, r, 2), AB=genrow(R,N), B=AB[2], dn=1, hasC=(i<=12), cf=0, Rn=1, eqRd=1, wit=List());
  if(hasC, cf=CMS[i]);
  for(n=1,N, dn=lcm(dn,n);
    if(hasC, my(t=n^r, g=gcd(t, cf(n))); Rn=lcm(Rn,t/g); if(Rn!=dn^r, eqRd=0));
    if(denominator(B[n+1])==dn^k, if(#wit<8, listput(wit,n))));
  print(R[1], " k=", k, " : R_n = d_n^r for all n<=", N, "? ", if(hasC, if(eqRd,"YES","NO"), "n/a"), "   n with den(b_n)=d_n^k : ", Vec(wit)));
}
quit;
