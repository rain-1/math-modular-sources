default(parisize,"3G");
read("lib.gp");
read("src.gp");
N = 3000;
{ for(i=1,12,
  my(R=ROWS[i], r=R[2], AB=genrow(R,N), B=AB[2], dn=1, Rn=1, cf=CMS[i]);
  my(mxT=0, mxTn=0, mxU=0, mxUn=0, eqR=0, eqD=0, mxRd=0, T, U, dv, dr);
  for(n=1,N,
    dn = lcm(dn,n); dr = dn^r;
    my(t=n^r, g=gcd(t, cf(n))); Rn = lcm(Rn, t/g);
    dv = denominator(B[n+1]);
    T = Rn/dv; U = dr/dv;
    if(#digits(T)>mxT, mxT=#digits(T); mxTn=n);
    if(#digits(U)>mxU, mxU=#digits(U); mxUn=n);
    if(T==1, eqR++); if(U==1, eqD++);
    if(#digits(dr/Rn)>mxRd, mxRd=#digits(dr/Rn));
  );
  print(R[1], " r=", r, " | max digits(R_n/den) = ", mxT, " at n=", mxTn, " | max digits(d_n^r/den) = ", mxU, " at n=", mxUn, " | #{n: den=R_n} = ", eqR, "/", N, " | #{n: den=d_n^r} = ", eqD, "/", N, " | max digits(d_n^r/R_n) = ", mxRd));
}
quit;
