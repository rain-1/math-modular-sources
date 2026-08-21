/* dwork.gp -- the genuine Dwork unit root lambda_0(t) = F(t)/F(t^sigma),
   computed as lim_s F_{<p^s}(t)/F_{<p^{s-1}}(t^p) at Teichmuller points,
   for comparison with the tower ratio u_a = lim A(ap^{s+1})/A(ap^s) = 1.  */
{
dwork(aa,bb,cc,typ, p, smax, K) =
  my(N,A,Ap,An,M,S,t0,tt,tp,pw,pwp,res,ss);
  N = p^smax;
  M = p^K;
  print("== Dwork unit root, row (",aa,",",bb,",",cc,") type R",typ,"  p=",p,
        "  truncation p^",smax,"  K=",K);
  /* partial sums of F at t and at t^p, cut at every p^s */
  for(j=1,p-1,
    t0 = lift(Mod(j,M)); tt = t0;
    for(i=1,K+2, tt = lift(Mod(tt,M)^p));      /* Teichmuller lift of j */
    tp = lift(Mod(tt,M)^p);
    Ap = 1; An = bb;                            /* A_0, A_1 exact integers */
    S = vector(smax+1); ss = vector(smax+1); S[1]=Mod(1,M); ss[1]=Mod(1,M);
    my(acc=Mod(1,M), accp=Mod(1,M), pt=Mod(tt,M), ptp=Mod(tp,M), n, k);
    /* n=0 term already in acc; add n>=1 */
    for(n=1,N-1,
      acc  = acc  + Mod(An,M)*pt;   pt  = pt *tt;
      accp = accp + Mod(An,M)*ptp;  ptp = ptp*tp;
      for(s=0,smax, if(n+1==p^s, S[s+1]=acc; ss[s+1]=accp));
      /* advance */
      if(typ==2,
        k = ((aa*n^2+aa*n+bb)*An - cc*n^2*Ap)/(n+1)^2,
        k = ((2*n+1)*(aa*n^2+aa*n+bb)*An - cc*n^3*Ap)/(n+1)^3);
      Ap=An; An=k;
    );
    S[smax+1]=acc; ss[smax+1]=accp;
    res = vector(smax);
    for(s=1,smax, res[s] = S[s+1]/ss[s]);
    print("  t0=",j,":  lambda_0 = ", lift(res[smax])+O(p^8),
          "   v(lambda_0-1)=", if(lift(res[smax])%p==1,
              valuation(lift(res[smax])-1+O(p^(K-2)),p), 0),
          "   stab: ", vector(smax-1,i, my(d=lift(res[i+1]-res[i]));
                              if(d==0,K,valuation(d+O(p^(K-2)),p))));
  );
}
