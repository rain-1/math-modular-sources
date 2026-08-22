/* ================================================================
   cyops.gp -- generic machinery for Calabi-Yau / MUM differential
   operators: local series solutions, Frobenius bases, Apery-type
   companions, archimedean and p-adic limits.

   An operator is given in THETA FORM at t=0:
       L = sum_{k=0}^{r} c_k(t) * theta^k ,   theta = t d/dt
   as a vector C = [c_0,...,c_r] of polynomials in the variable t.
   ================================================================ */

/* --- Stirling numbers of the second kind: theta^k = sum_j S(k,j) t^j D^j */
S2(k,j) = stirling(k,j,2);

/* theta-form  ->  D-form:  L = sum_j P_j(t) D^j   */
thetaToD(C) =
{ my(r=#C-1, P=vector(r+1,i,0));
  for(j=0,r,
    my(acc=0);
    for(k=j,r, acc += S2(k,j)*C[k+1]);
    P[j+1] = 't^j * acc);
  P;
}

/* shift the D-form to expansion point z0: variable s = t - z0 */
shiftD(P,z0) = vector(#P, i, subst(P[i], 't, z0+'s));

/* falling factorial x(x-1)...(x-j+1) */
ff(x,j) = prod(i=0,j-1, x-i);

/* Given D-form coefficients Pj (polys in s), build the table of terms
     [shift=j-i, j, i, p_{j,i}]  with
     [s^m] L(y) = sum_{j,i} p_{j,i} * ff(m-i+j, j) * u_{m-i+j}.
   d = max shift.                                                   */
recData(Pj) =
{ my(T=List(), d=-10^9);
  for(jj=1,#Pj,
    my(j=jj-1, pol=Pj[jj]);
    if(pol==0, next);
    for(i=0, poldegree(pol,'s),
      my(p=polcoef(pol,i,'s));
      if(p==0, next);
      listput(T,[j-i,j,i,p]);
      if(j-i>d, d=j-i)));
  [d, Vec(T)];
}

/* leading polynomial N(M) = sum_{j-i=d} p_{j,i} ff(M,j)  */
leadPoly(RD) =
{ my(d=RD[1], T=RD[2], acc=0);
  for(a=1,#T, if(T[a][1]==d, acc += T[a][4]*ff('M,T[a][2])));
  acc;
}

/* Solve L(y)=RHS as sum_{n=0}^{N} u_n s^n.
   ini/imp : imposed coefficients (imp[n+1]=1 => u_n := ini[n+1]).
   rhs, roff : rhs[m+roff+1] = [s^m] RHS   (m from -roff up).
   Returns [U, OBS, FREE].                                          */
solveSeries(RD, N, ini, imp, rhs, roff) =
{ my(d=RD[1], T=RD[2], U=vector(N+1,i,0), OBS=List(), FREE=List(),
     Lp=leadPoly(RD));
  for(i=1,N+1, if(i<=#imp && imp[i], U[i]=ini[i]));
  for(m=-d-roff, N-d,
    my(M=m+d, lead=subst(Lp,'M,M), acc=0);
    if(M<0,
      /* equation with no unknown of index >=0 on the lead term */
      lead=0);
    for(a=1,#T,
      my(sh=T[a][1], j=T[a][2], p=T[a][4], idx=m+sh);
      if(sh==d && M>=0, next);
      if(idx<0 || idx>N, next);
      acc += p*ff(idx,j)*U[idx+1]);
    my(r = if(m+roff+1>=1 && m+roff+1<=#rhs, rhs[m+roff+1], 0));
    if(M>=0 && M<=#imp-1 && imp[M+1],
       my(res = lead*U[M+1] + acc - r);
       if(res!=0, listput(OBS,[m,res]));
      ,
      if(lead==0,
        if(acc-r != 0, listput(OBS,[m,acc-r]), if(M>=0,listput(FREE,M)));
        ,
        U[M+1] = (r-acc)/lead)));
  [U, Vec(OBS), Vec(FREE)];
}

/* --- Laurent-series helper: G(s) for the log solution ------------
   Given D-form Pj (polys in s) and V (vector of coefficients of the
   holomorphic factor, V[n+1] = v_n), compute the coefficients of
   G = sum_{j>=1} P_j * sum_{i=1}^{j} binom(j,i)(-1)^{i-1}(i-1)! s^{-i} V^{(j-i)}
   as a vector indexed from s^{-K}.  Returns [Gvec, K].            */
logRHS(Pj, V, N, K) =
{ my(G=vector(N+K+1,i,0));   /* G[m+K+1] = [s^m] G, m=-K..N */
  for(jj=2,#Pj,
    my(j=jj-1, pol=Pj[jj]);
    if(pol==0, next);
    for(i=1,j,
      my(cf = binomial(j,i)*(-1)^(i-1)*(i-1)!, e=j-i);
      /* V^{(e)} has coefficients: [s^n] V^{(e)} = ff(n+e,e) v_{n+e} */
      for(ip=0,poldegree(pol,'s),
        my(pp=polcoef(pol,ip,'s));
        if(pp==0,next);
        for(n=0,N,
          if(n+e>N, break);
          my(val = cf*pp*ff(n+e,e)*V[n+e+1], m = n+ip-i);
          if(val==0, next);
          if(m+K+1>=1 && m+K+1<=#G, G[m+K+1] += val,
             if(m+K+1<1, error("logRHS: need larger K, m=",m)))))));  
  [G,K];
}
