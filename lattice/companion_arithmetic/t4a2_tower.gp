default(parisize,"8G");
read("lib.gp");
NT = 6600;
IDX = [8, 4];
{ for(u=1,2,
  my(i=IDX[u], R=ROWS[i], w=R[2], AB=genrow(R,NT), A=AB[1], B=AB[2]);
  print("\n===== row ", R[1], "  (a,b,c)=(", R[3], ",", R[4], ",", R[5], ")  w=r=", w, "  (c=", R[5], ": all primes good) =====");
  print("   Lam_k := p^{wk} b_{np^k}/a_{np^k} ;  e(k) := v_p(Lam_{k+1}/Lam_k - 1) ;  rhoA_k=a_{np^{k+1}}/a_{np^k} ; rhoB_k=b_{np^{k+1}}/b_{np^k}");
  for(pi=1,4,
    my(p=[5,7,11,13][pi], KM=floor(log(NT/3.0)/log(p*1.0)));
    print("  p=", p);
    for(n=1,3,
      my(L0="", L1="", L2="", L3="");
      for(k=0,KM-1,
        my(j1=n*p^k, j2=n*p^(k+1), r1, r2, ra, rb);
        if(j2<=NT,
          r1 = B[j1+1]/A[j1+1]; r2 = B[j2+1]/A[j2+1];
          ra = A[j2+1]/A[j1+1]; rb = B[j2+1]/B[j1+1];
          L0 = Str(L0, " ", valuation(p^(w*k)*r1,p));
          L1 = Str(L1, " ", valuation(p^w*r2/r1 - 1, p));
          L2 = Str(L2, " ", valuation(ra-1,p));
          L3 = Str(L3, " ", valuation(p^w*rb-1,p))));
      print("    n=", n, " : v_p(Lam_k) =", L0, "   | e(k) =", L1, "   | v_p(rhoA_k - 1) =", L2, "   | v_p(p^w rhoB_k - 1) =", L3))));
}
quit;
