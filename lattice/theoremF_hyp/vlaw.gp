/* vlaw.gp -- exact valuation laws for the census rows.
   Tests, for each row (A_n) and slope prime p:
     (L1) digit law   v_p(a_n) ?= sum_i v_p(a_{n_i}),  n = sum n_i p^i
     (L2) Lucas for the normalised sequence:
          ahat_n := a_n / p^{L(n)},  L(n)=sum_i v_p(a_{n_i});
          test  ahat_n = prod_i ahat_{n_i} mod p
     (L3) supercongruence  a_{m p^k} = a_{m p^{k-1}} mod p^{3k}   (m<=8)
     (L4) constancy of v_p(a_{m p^k}) in k
     (L5) sharp linear fit  v_p(a_n) <= lam * s_p(n) + C  (best lam over n<=M) */

sdigv(n,p,va) = {my(d=digits(n,p), s=0); for(i=1,#d, s += va[d[i]+1]); s;}

vlaw(A, p, name, M) =
{ my(va, v, L, ok1=1, bad1=List(), ok2=1, bad2=List(), mx=0, mn=0, lam, C, S=List());
  va = vector(p, j, valuation(A[j], p));      /* va[j] = v_p(a_{j-1}) */
  v  = vector(M+1, n, valuation(A[n], p));    /* v[n+1] = v_p(a_n) */
  print("=== ", name, "  p=", p, "   v_p(a_r), r=0..p-1: ", va);
  /* (L1) */
  for(n=1, M, my(d=v[n+1]-sdigv(n,p,va));
      if(d>mx, mx=d); if(d<mn, mn=d);
      if(d!=0 && #bad1<10, listput(bad1,[n,v[n+1],sdigv(n,p,va)])));
  print("  (L1) v_p(a_n) - sum_i v_p(a_{n_i}) in [", mn, ",", mx, "]",
        if(mx==0 && mn==0, "   EXACT DIGIT LAW", ""));
  if(mx!=0 || mn!=0, print("       first deviations [n,v,pred]: ", Vec(bad1)));
  /* (L2) Lucas for normalised sequence, mod p */
  for(n=1, M, my(d=digits(n,p), pr=Mod(1,p), lhs);
      if(v[n+1] != sdigv(n,p,va), ok2=0; if(#bad2<8, listput(bad2,n)); next);
      for(i=1,#d, pr *= Mod(A[d[i]+1]/p^va[d[i]+1], p));
      lhs = Mod(A[n+1]/p^sdigv(n,p,va), p);
      if(lhs != pr, ok2=0; if(#bad2<8, listput(bad2,n))));
  print("  (L2) Lucas for a_n/p^{L(n)} mod p: ", if(ok2, "HOLDS", concat("FAILS at n=", Str(Vec(bad2)))));
  /* (L5) */
  lam = 0.0; for(n=1, M, my(s=vecsum(digits(n,p))); if(s>0, my(r=1.0*v[n+1]/s); if(r>lam, lam=r)));
  C = 0; for(n=1, M, my(d=v[n+1]-ceil(lam)*vecsum(digits(n,p))); if(d>C, C=d));
  print("  (L5) max v_p(a_n)/s_p(n) = ", lam, ";  v_p(a_n) <= ", ceil(lam), "*s_p(n) + ", C, " for n<=", M);
}

super(A, p, name, M) =
{ my(K=floor(log(M)/log(p)), out=List());
  print("  (L3/L4) m: [v_p(a_{m p^k})]_{k=0..}  |  min k-level of a_{mp^k}-a_{mp^{k-1}} minus 3k");
  for(m=1, 8,
    my(vs=List(), sc=List());
    for(k=0, K, if(m*p^k<=M, listput(vs, valuation(A[m*p^k+1],p))));
    for(k=1, K, if(m*p^k<=M, my(dd=A[m*p^k+1]-A[m*p^(k-1)+1]); listput(sc, if(dd==0, "inf", valuation(dd,p)-3*k))));
    print("    m=", m, ": v=", Vec(vs), "   super-3k slack=", Vec(sc)));
}
