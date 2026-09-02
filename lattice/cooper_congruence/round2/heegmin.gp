\\ heegmin.gp -- MINIMAL-A Gamma_0(N)-Heegner representative.  heegrep of heeg.gp returns the
\\ FIRST (p,r) found, whose A can be huge; at Im(alpha) ~ 1e-4 the evaluation of fhat is
\\ numerically worthless.  This version minimises A over the search box.
{ heegminA(Q,N,beta,PB) = my(a=Q[1],b=Q[2],c=Q[3], A,B,C,g,p,q,r,s, bA=0, bB=0);
  for(r=0,PB,
    for(p=-PB,PB,
      if(gcd(p,r)!=1, next);
      A = a*p^2 + b*p*r + c*r^2;
      if(A<=0 || A%N!=0, next);
      if(bA>0 && A>=bA, next);
      g = bezout(p,-r);
      if(g[3]!=1, next);
      s = g[1]; q = g[2];
      B = 2*a*p*q + b*(p*s+q*r) + 2*c*r*s;
      if((B-beta)%(2*N)!=0, next);
      bA = A; bB = B));
  if(bA==0, return(0));
  C = (bB^2 - (Q[2]^2-4*Q[1]*Q[3]))/(4*bA);
  [bA,bB,C];
}
