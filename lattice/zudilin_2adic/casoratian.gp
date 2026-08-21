/* symbolic Beukers Pade polynomials in X */
polpq(N) =
{ my(P=vector(N+2), Q=vector(N+2));
  Q[1]=1; Q[2]=X^2-X+1; P[1]=0; P[2]=1;
  for(m=1,N,
    Q[m+2]=((2*m*(m+1)+1-X+X^2)*Q[m+1]-m^2*Q[m])/(m+1)^2;
    P[m+2]=((2*m*(m+1)+1-X+X^2)*P[m+1]-m^2*P[m])/(m+1)^2);
  [P,Q];
}
N=14; b=polpq(N); P=b[1]; Q=b[2];
print("m : X^2*W_m(X)  (should be deg<=2)");
{ for(m=1,N,
  my(pm=P[m+1], qm=Q[m+1], pm1=subst(P[m],X,X+1), qm1=subst(Q[m],X,X+1));
  my(W = X^2*(pm*qm1 + pm1*qm) - 2*qm*qm1);
  print("  m=",m,"  degree=",poldegree(W),"   W=",W));
}
quit;
