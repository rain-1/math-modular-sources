/* lib.gp -- mirror map (t(q),F(q)) from an order-3 MUM recurrence alone.
   Recurrence:  p3(n) A_{n+1} = p2(n) A_n + p1(n) A_{n-1},  A_0=1, A_{-1}=0.
   p3,p2,p1 are polynomials in 'n with p3(n)=(n+1)^3 (MUM, triple indicial root 0).
   Frobenius deformation A_n(e) (n -> n+e), y0 = sum A_n(0) t^n,
   g = sum A_n'(0) t^n,  q = t exp(g/y0),  t(q) = reversion,  F(q) = y0(t(q)).   */

Aeps(p3,p2,p1,NT)={
  my(A=vector(NT+1), E='e+O('e^2), P3, P2, P1);
  A[1] = 1 + O('e^2);
  P3 = subst(p3,'n,0+E); P2 = subst(p2,'n,0+E);
  A[2] = P2*A[1]/P3;
  for(k=1, NT-1,
    P3 = subst(p3,'n,k+E); P2 = subst(p2,'n,k+E); P1 = subst(p1,'n,k+E);
    A[k+2] = (P2*A[k+1] + P1*A[k])/P3;
  );
  A;
}

mirror(p3,p2,p1,NT)={
  my(A=Aeps(p3,p2,p1,NT), y0, g, qs, tq, Fq);
  y0 = Ser(vector(NT+1,i,polcoef(A[i],0,'e)),'t,NT+1);
  g  = Ser(vector(NT+1,i,polcoef(A[i],1,'e)),'t,NT+1);
  qs = 't*exp(g/y0);
  tq = serreverse(qs);
  Fq = subst(y0,'t,tq);
  [tq, Fq, y0, g];
}

/* plain integer sequence from the recurrence, exact, to n=NN */
seqA(p3,p2,p1,NN)={
  my(A=vector(NN+1));
  A[1]=1; A[2]=subst(p2,'n,0)/subst(p3,'n,0);
  for(k=1,NN-1, A[k+2]=(subst(p2,'n,k)*A[k+1]+subst(p1,'n,k)*A[k])/subst(p3,'n,k));
  A;
}
