polpq(N) =
{ my(P=vector(N+2), Q=vector(N+2));
  Q[1]=1; Q[2]=X^2-X+1; P[1]=0; P[2]=1;
  for(m=1,N,
    Q[m+2]=((2*m*(m+1)+1-X+X^2)*Q[m+1]-m^2*Q[m])/(m+1)^2;
    P[m+2]=((2*m*(m+1)+1-X+X^2)*P[m+1]-m^2*P[m])/(m+1)^2);
  [P,Q];
}
N=16; b=polpq(N+1); P=b[1]; Q=b[2];
sh(f,c)=subst(f,X,X+c);
print("n : X^2(X-1)^2 * Omega'_n(X)");
{ for(n=1,N,
  my(pp=sh(P[n+2],-1), qq=sh(Q[n+2],-1), pm=sh(P[n],1), qm=sh(Q[n],1));
  my(W = X^2*(X-1)^2*(pp*qm - pm*qq) - (2*X^2 - 2*(X-1)^2)*qq*qm);
  print("  n=",n," deg=",poldegree(W),"  W=",W));
}
quit;
