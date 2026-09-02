default(parisize,"8G");
read("lib.gp");
NT = 6600;
IDX = [8, 4];
vv(x,p) = if(x==0, "inf", valuation(x,p));
{ for(u=1,2,
  my(i=IDX[u], R=ROWS[i], w=R[2], AB=genrow(R,NT), A=AB[1], B=AB[2]);
  print("\n===== row ", R[1], "  (a,b,c)=(", R[3], ",", R[4], ",", R[5], ")  w=r=", w, " ; c=", R[5], " so every prime is good =====");
  for(pi=1,4,
    my(p=[5,7,11,13][pi]);
    print("  p=", p);
    for(n=1,3,
      my(L1="", L2="", L3="");
      for(k=0,2,
        my(j1=n*p^k, j2=n*p^(k+1), r1, r2, d, dn);
        if(j2<=NT,
          r1 = B[j1+1]/A[j1+1]; r2 = B[j2+1]/A[j2+1];
          d = r2-r1; dn = p^(w*(k+1))*r2 - p^(w*k)*r1;
          L1 = Str(L1, "  k=", k, ": ", vv(d,p));
          L2 = Str(L2, "  k=", k, ": ", vv(dn,p));
          L3 = Str(L3, "  k=", k, ": ", vv(r1,p)),
          L1 = Str(L1, "  k=", k, ": n/a")));
      print("    n=", n, "   v_p(r_{np^{k+1}} - r_{np^k}) : ", L1);
      print("          v_p(p^{w(k+1)}r_{np^{k+1}} - p^{wk}r_{np^k}) : ", L2);
      print("          v_p(r_{np^k}) : ", L3))));
}
quit;
