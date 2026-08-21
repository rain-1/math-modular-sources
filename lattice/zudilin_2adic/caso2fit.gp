polpq(N) =
{ my(P=vector(N+2), Q=vector(N+2));
  Q[1]=1; Q[2]=X^2-X+1; P[1]=0; P[2]=1;
  for(m=1,N,
    Q[m+2]=((2*m*(m+1)+1-X+X^2)*Q[m+1]-m^2*Q[m])/(m+1)^2;
    P[m+2]=((2*m*(m+1)+1-X+X^2)*P[m+1]-m^2*P[m])/(m+1)^2);
  [P,Q];
}
NM=26; b=polpq(NM+1); P=b[1]; Q=b[2];
sh(f,c)=subst(f,X,X+c);
Wpol(n) = { my(pp=sh(P[n+2],-1), qq=sh(Q[n+2],-1), pm=sh(P[n],1), qm=sh(Q[n],1));
  X^2*(X-1)^2*(pp*qm - pm*qq) - (2*X^2 - 2*(X-1)^2)*qq*qm; }
Ws = vector(NM, n, Wpol(n));
/* fit each coefficient c_j(n)*n^2*(n+1)^2 as polynomial in n */
print("coefficient fits  w_j(n) := n^2(n+1)^2 * [X^j] W_n :");
{ for(j=0,6,
   my(v=vector(12,n, n^2*(n+1)^2*polcoef(Ws[n],j)));
   my(f = polinterpolate(vector(12,i,i), v, 'n));
   my(ok = 1);
   for(n=1,NM, if(subst(f,'n,n) != n^2*(n+1)^2*polcoef(Ws[n],j), ok=0));
   print("  j=",j,"  w_j(n) = ",f, "   verified n<=",NM,": ",ok));
}
/* now check W_n(x_n) == q(n)/(64 n^2 (n+1)^2) */
qz(n) = 3520*n^6+5632*n^5+2064*n^4-384*n^3-156*n^2+16*n+7;
print("");
{ my(ok=1); for(n=1,NM, if(subst(Ws[n],X,1/2-n) != qz(n)/(64*n^2*(n+1)^2), ok=0));
  print("W_n(1/2-n) == q(n)/(64 n^2 (n+1)^2) for n<=",NM,": ",ok); }
/* also Lemma A closed form check */
Om(n) = { my(pn=P[n+1], qn=Q[n+1], pm=sh(P[n],1), qm=sh(Q[n],1)); X^2*(pn*qm+pm*qn)-2*qn*qm; }
{ my(ok=1); for(n=1,NM, if(Om(n) != -(X^2/n^2 - 2*X/n + 2), ok=0));
  print("LemmaA: X^2*Omega_n == -(X^2/n^2-2X/n+2) for n<=",NM,": ",ok); }
{ my(ok=1); for(n=1,NM, if(subst(Om(n),X,1/2-n)/(1/2-n)^2 != -(20*n^2-8*n+1)/(n^2*(2*n-1)^2), ok=0));
  print("LemmaA at x_n gives -p(n)/(n^2(2n-1)^2) for n<=",NM,": ",ok); }
/* W_n(1/2) closed form */
{ my(ok=1); for(n=1,NM, if(subst(Ws[n],X,1/2) != (8*n^2+8*n+7)/(64*n^2*(n+1)^2), ok=0));
  print("W_n(1/2) == (8n^2+8n+7)/(64 n^2(n+1)^2) for n<=",NM,": ",ok); }
quit;
