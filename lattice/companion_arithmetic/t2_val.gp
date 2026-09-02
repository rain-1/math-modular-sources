default(parisize,"4G");
read("lib.gp");
N = 3000;
vsafe(x,p) = if(x==0, -1, valuation(x,p));
{ for(i=1,#ROWS,
  my(R=ROWS[i], r=R[2], AB=genrow(R,N), A=AB[1], B=AB[2], bp=badprimes(R));
  if(#bp==0, print("\n### ", R[1], " (a,b,c,d)=(", R[3], ",", R[4], ",", R[5], ",", R[6], ")  r=", r, "  c=+-1: NO bad primes"), 0);
  for(j=1,#bp,
    my(p=bp[j], VA=vector(N), VB=vector(N), SP=vector(N), zer=List());
    for(n=1,N,
      if(A[n+1]==0, listput(zer,n); VA[n]=-1, VA[n]=valuation(A[n+1],p));
      VB[n]=valuation(B[n+1],p); SP[n]=sumdigits(n,p));
    my(lam=0, exact=1, firstfail=0, mn=0, mx=0);
    for(n=1,N, if(VA[n]>=0, my(q=VA[n]/SP[n]); if(q>lam, lam=q)));
    my(lamI=ceil(lam), dmin=10^9, dmax=-10^9, argmin=0, argmax=0);
    for(n=1,N, if(VA[n]>=0, my(dd=VA[n]-lamI*SP[n]); if(dd<dmin, dmin=dd; argmin=n); if(dd>dmax, dmax=dd; argmax=n)));
    my(exlam=0, exfail=0);
    for(n=1,N, if(VA[n]>=0 && VA[n]!=lamI*SP[n], exfail=n; break));
    my(dif=vector(N,n,VB[n]-VA[n]), cst=dif[N], badn=List(), nb=0, lastdiff=0);
    for(n=1,N, if(VA[n]>=0 && dif[n]!=cst, nb++; lastdiff=n; if(#badn<10, listput(badn,n))));
    print("\n### ", R[1], " (a,b,c,d)=(", R[3], ",", R[4], ",", R[5], ",", R[6], ")  r=", r, "   p=", p, "   v_p(c)=", valuation(R[5],p));
    if(#zer>0, print("   a_n = 0 at n = ", Vec(zer)));
    print("   v_p(a_n): sharp lambda = max v_p(a_n)/s_p(n) = ", lam, " ; least integer lambda with v_p(a_n) <= lambda s_p(n): ", lamI);
    print("   v_p(a_n) - ", lamI, "*s_p(n) ranges in [", dmin, ",", dmax, "]  (min at n=", argmin, ", max at n=", argmax, ")");
    print("   exact law v_p(a_n) = ", lamI, "*s_p(n) for all n<=", N, ": ", if(exfail==0, "YES", Str("NO, first failure n=", exfail)));
    print("   v_p(b_n)-v_p(a_n): value at n=", N, " is ", cst, " ; #{n<=", N, ": differs} = ", nb, " ; last n where it differs = ", lastdiff, " ; first few: ", Vec(badn));
    print("   v_p(b_n)-v_p(a_n) at n=100,500,1000,2000,3000: ", [dif[100],dif[500],dif[1000],dif[2000],dif[3000]]);
  ));
}
quit;
