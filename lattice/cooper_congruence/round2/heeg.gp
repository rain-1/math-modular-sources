\\ heeg.gp -- Heegner forms of discriminant d on Gamma_0(N) with B = beta mod 2N.

\\ all SL_2(Z)-reduced positive definite forms [a,b,c] of discriminant d<0 (incl. imprimitive)
{ redforms(d) = my(L=List(), amax, c);
  amax = sqrtint(-d\3);
  for(a=1, amax,
    for(b=-a, a,
      if((b^2-d)%(4*a)!=0, next);
      c = (b^2-d)/(4*a);
      if(c<a, next);
      if((a==c || b==-a) && b<0, next);
      listput(L,[a,b,c])));
  Vec(L);
}

\\ given reduced form Q=[a,b,c] of disc d, find a Gamma_0(N)-Heegner representative
\\ [A,B,C] with N|A, B = beta mod 2N.  Returns [A,B,C] or 0.
{ heegrep(Q,N,beta,PB) = my(a=Q[1],b=Q[2],c=Q[3], A,B,C,g,p,q,r,s);
  for(r=0,PB,
    for(p=-PB,PB,
      if(gcd(p,r)!=1, next);
      A = a*p^2 + b*p*r + c*r^2;
      if(A<=0 || A%N!=0, next);
      g = bezout(p,-r);   \\ g = [x,y,gcd]: p*x + (-r)*y = 1
      if(g[3]!=1, next);
      s = g[1]; q = g[2];  \\ p*s - r*q = 1
      B = 2*a*p*q + b*(p*s+q*r) + 2*c*r*s;
      if((B-beta)%(2*N)!=0, next);
      C = (B^2 - (b^2-4*a*c))/(4*A);
      return([A,B,C])));
  0;
}

\\ order of the stabiliser of alpha_Q in Gamma_0(N)-bar (Heegner conventions)
{ omeg(Q) = my(g=gcd(gcd(Q[1],Q[2]),Q[3]), dd);
  dd = (Q[2]^2-4*Q[1]*Q[3])/g^2;
  if(dd==-3, return(3));
  if(dd==-4, return(2));
  1;
}

\\ genus character chi_{D0}([A,B,C]), D0 fundamental, disc(Q) = D0*d
{ genchar(Q,D0) = my(A=Q[1],B=Q[2],C=Q[3]);
  if(gcd(A,D0)==1, return(kronecker(D0,A)));
  if(gcd(C,D0)==1, return(kronecker(D0,C)));
  if(gcd(A+B+C,D0)==1, return(kronecker(D0,A+B+C)));
  0;
}
