default(parisize,"10G");
read("lib.gp");
N = 3000;
PL = [2,3,5,7,11,13,17,19,23,29,31,37,41,43];
{ for(i=1,#ROWS,
  my(R=ROWS[i], AB=genrow(R,N), A=AB[1], B=AB[2], out="");
  for(pi=1,#PL,
    my(p=PL[pi], VA=vector(N), SP=vector(N), bw=10^9, ba=0, bc=0);
    for(n=1,N, VA[n]=valuation(A[n+1],p); SP[n]=sumdigits(n,p));
    for(lam=0,6, my(mn=10^9, mx=-10^9);
      for(n=1,N, my(v=VA[n]-lam*SP[n]); if(v<mn,mn=v); if(v>mx,mx=v));
      if(mx-mn<bw, bw=mx-mn; ba=lam; bc=mn));
    if(bw==0, out=Str(out, "  p=", p, ": v_p(a_n)=", ba, "*s_p(n)", if(bc==0,"",Str("+",bc)))));
  print(R[1], " exact digit laws for a_n (p<=43, n<=", N, "):", if(out=="", " NONE", out)));
}
quit;
