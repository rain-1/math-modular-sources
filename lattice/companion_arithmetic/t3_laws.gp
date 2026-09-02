default(parisize,"4G");
read("lib.gp");
N = 3000;
/* best fit T(n) = AA*n + BB*s_p(n) + CC over small ranges; report width */
{ bestfit(T, SP, N) =
  my(bw=10^9, ba=0, bb=0, bc=0, mn, mx);
  for(AA=-8,8, for(BB=-8,8,
    mn=10^9; mx=-10^9;
    for(n=1,N, my(v=T[n]-AA*n-BB*SP[n]); if(v<mn,mn=v); if(v>mx,mx=v));
    if(mx-mn<bw, bw=mx-mn; ba=AA; bb=BB; bc=mn)));
  [ba,bb,bc,bw];
}
/* does res(n) depend only on n mod p^j ? */
{ modtest(T, SP, AA, BB, p, N) =
  my(out=[]);
  for(j=1,4, my(md=p^j, tab=vector(md,x,[]), ok=1);
    for(n=1,N, my(cl=(n%md)+1, v=T[n]-AA*n-BB*SP[n]);
      if(tab[cl]==[], tab[cl]=v, if(tab[cl]!=v, ok=0; break)));
    if(ok, out=concat(out,[j])));
  out;
}
{ for(i=1,#ROWS,
  my(R=ROWS[i], AB=genrow(R,N), A=AB[1], B=AB[2], bp=badprimes(R));
  for(j=1,#bp,
    my(p=bp[j], VA=vector(N), VB=vector(N), SP=vector(N));
    for(n=1,N, VA[n]=if(A[n+1]==0, 0, valuation(A[n+1],p)); VB[n]=valuation(B[n+1],p); SP[n]=sumdigits(n,p));
    my(fa=bestfit(VA,SP,N), fb=bestfit(VB,SP,N));
    print("\n### ", R[1], "  p=", p);
    print("   v_p(a_n) best fit: ", fa[1], "*n + ", fa[2], "*s_p(n) + [", fa[3], ",", fa[3]+fa[4], "]   width=", fa[4], if(fa[4]==0, "  <-- EXACT", ""));
    if(fa[4]>0 && fa[4]<=6, print("      residual constant on n mod p^j for j in ", modtest(VA,SP,fa[1],fa[2],p,N)));
    print("   v_p(b_n) best fit: ", fb[1], "*n + ", fb[2], "*s_p(n) + [", fb[3], ",", fb[3]+fb[4], "]   width=", fb[4], if(fb[4]==0, "  <-- EXACT", ""));
    if(fb[4]>0 && fb[4]<=6, print("      residual constant on n mod p^j for j in ", modtest(VB,SP,fb[1],fb[2],p,N)));
  ));
}
quit;
