default(parisize,"10G");
read("lib.gp");
read("src.gp");
N = 3000;
print("=== (a) zero-limit rows: DIF = v_p(b_n)-v_p(a_n) versus sigma_p * n ===");
{ for(j=1,2,
  my(i=[1,11][j], p=[2,3][j], sg=[3,3][j], R=ROWS[i], AB=genrow(R,N), A=AB[1], B=AB[2], mn=10^9, mx=-10^9);
  for(n=1,N, my(d=valuation(B[n+1],p)-valuation(A[n+1],p)-sg*n); if(d<mn,mn=d); if(d>mx,mx=d));
  print("  ", R[1], " p=", p, " : DIF - ", sg, "n in [", mn, ",", mx, "]  (bounded => v_p(b_n/a_n) = sigma_p n + O(1), xi_p=0)"));
}
print("");
print("=== (b) Cooper rows with no p-adic limit: exact law for DIF ===");
{ for(j=1,3,
  my(i=[13,14,15][j], p=[3,2,2][j], R=ROWS[i], AB=genrow(R,N), A=AB[1], B=AB[2], mn=10^9, mx=-10^9, ok=1, bad=0);
  for(n=1,N, my(d=valuation(B[n+1],p)-valuation(A[n+1],p)+2*logint(n,p)); if(d<mn,mn=d); if(d>mx,mx=d));
  print("  ", R[1], " p=", p, " : DIF + 2*floor(log_p n) in [", mn, ",", mx, "]"));
}
print("");
print("=== (c) equality of the two Cooper 2-adic valuation vectors ===");
{ my(A1=genrow(ROWS[14],N), A2=genrow(ROWS[15],N), okA=1, okB=1);
  for(n=1,N, if(valuation(A1[1][n+1],2)!=valuation(A2[1][n+1],2), okA=0); if(valuation(A1[2][n+1],2)!=valuation(A2[2][n+1],2), okB=0));
  print("  v_2(a_n) identical for s10 and s18, n<=", N, ": ", if(okA,"YES","NO"), " ; v_2(b_n) identical: ", if(okB,"YES","NO"));
}
print("");
print("=== (d) where the deficiency R_n / den(b_n) lives: bad primes vs good primes ===");
{ for(i=1,12,
  my(R=ROWS[i], r=R[2], AB=genrow(R,N), B=AB[2], cf=CMS[i], Rn=1, bp=badprimes(R), nb=0, ng=0, T, dv, fa);
  for(n=1,N,
    my(t=n^r, g=gcd(t, cf(n))); Rn = lcm(Rn, t/g);
    dv = denominator(B[n+1]); T = Rn/dv;
    if(T>1, fa=factor(T); for(k=1,#fa[,1], if(setsearch(Set(bp), fa[k,1]), nb+=fa[k,2], ng+=fa[k,2]))));
  print("  ", R[1], " bad primes ", bp, " : total deficiency exponent at bad primes = ", nb, " ; at good primes = ", ng));
}
quit;
