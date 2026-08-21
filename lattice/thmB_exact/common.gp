/* common.gp -- shared machinery for the Theorem B exact-prefactor computation.
   Builds (t,F) from the recurrence alone (Frobenius + nome), the Eisenstein
   source Phi from its closed-form divisor coefficients, and the Eichler
   integral Theta = D^{-(w+1)} Phi.                                          */


/* ---- the twelve rows: [name, r, a, b, c] ------------------------------- */
{ROWS = [
 ["A",     2, 7, 2, -8],
 ["B",     2, 9, 3, 27],
 ["C",     2,10, 3,  9],
 ["D",     2,11, 3, -1],
 ["E",     2,12, 4, 32],
 ["F",     2,17, 6, 72],
 ["alpha", 3,10, 4, 64],
 ["gamma", 3,17, 5,  1],
 ["delta", 3, 7, 3, 81],
 ["eps",   3,12, 4, 16],
 ["zeta",  3, 9, 3,-27],
 ["eta",   3,11, 5,125]
];}

frob(r,a,b,c,NT) = {
  my(A=vector(NT+1)); A[1]=1; A[2]=b;
  if(r==2, for(n=1,NT-1, A[n+2] = ((a*n^2+a*n+b)*A[n+1] - c*n^2*A[n])/(n+1)^2),
           for(n=1,NT-1, A[n+2] = ((2*n+1)*(a*n^2+a*n+b)*A[n+1] - c*n^3*A[n])/(n+1)^3));
  A;
}

/* exact q-expansions of t and F to O(q^(NT+1)); series variable 't */
build(r,a,b,c,NT) = {
  my(A=frob(r,a,b,c,NT), y0=Ser(A,'t,NT+1), P, K, g, qs, tq, Pq, Kq, Fq);
  if(r==2, P = 1 - a*'t + c*'t^2 + O('t^(NT+1)); K = P*y0^2,
           P = 1 - 2*a*'t + c*'t^2 + O('t^(NT+1)); K = y0*sqrt(P));
  g  = intformal((1/K - 1)/'t);
  qs = 't*exp(g); tq = serreverse(qs);
  Kq = ('t*deriv(tq,'t))/tq;
  if(r==2, Pq = 1 - a*tq + c*tq^2;   Fq = sqrt(Kq/Pq),
           Pq = 1 - 2*a*tq + c*tq^2; Fq = Kq/sqrt(Pq));
  [tq, Fq];
}

/* ---- closed-form Fourier coefficients c(m) of the Eisenstein source ----- */
/* sigma^{(k)}_chi(m) = sum_{d|m} chi(m/d) d^k     (inner placement)        */
sigin(chi,k,m) = if(m<1 || frac(m), 0, sumdiv(m, d, chi(m/d)*d^k));
/* tilde sigma^{(k)}_chi(m) = sum_{d|m} chi(d) d^k (outer placement)        */
sigout(chi,k,m)= if(m<1 || frac(m), 0, sumdiv(m, d, chi(d)*d^k));

ch1(n) = 1;
chm3(n)= kronecker(-3,n);
chm4(n)= kronecker(-4,n);
ch5(n) = kronecker(5,n);
nuD(d) = my(u=d%5); if(u==1,1, u==2,-2, u==3,2, u==4,-1, 0);

cm(nm, m) = {
  if(nm=="A",     sigout(chm3,2,m) - sigout(chm3,2,m/2),
     nm=="B",     sigin(chm3,2,m) - 6*sigin(chm3,2,m/2) - 8*sigin(chm3,2,m/4),
     nm=="C",     sigin(chm3,2,m) - 8*sigin(chm3,2,m/2),
     nm=="D",     sumdiv(m, d, nuD(d)*d^2),
     nm=="E",     sigin(chm4,2,m) - 8*sigin(chm4,2,m/2),
     nm=="F",     sigin(chm3,2,m) - 7*sigin(chm3,2,m/2) - 8*sigin(chm3,2,m/4),
     nm=="alpha", sigin(ch1,3,m) -17*sigin(ch1,3,m/2) -9*sigin(ch1,3,m/3)
                  +16*sigin(ch1,3,m/4) +153*sigin(ch1,3,m/6) -144*sigin(ch1,3,m/12),
     nm=="gamma", sigin(ch1,3,m) -28*sigin(ch1,3,m/2) +63*sigin(ch1,3,m/3) -36*sigin(ch1,3,m/6),
     nm=="delta", sigin(ch1,3,m) -14*sigin(ch1,3,m/2) -sigin(ch1,3,m/3)
                  +16*sigin(ch1,3,m/4) +14*sigin(ch1,3,m/6) -16*sigin(ch1,3,m/12),
     nm=="eps",   sigin(ch1,3,m) -21*sigin(ch1,3,m/2) +84*sigin(ch1,3,m/4) -64*sigin(ch1,3,m/8),
     nm=="zeta",  sumdiv(m, d, chm3(d)*chm3(m/d)*d^3),
     nm=="eta",   sigin(ch5,3,m) -14*sigin(ch5,3,m/2) -16*sigin(ch5,3,m/4),
     error("unknown row ", nm));
}

/* the Mellin polynomial P(s) = sum_d c_d d^{-s} of the oldform combination  */
Pmellin(nm, s) = {
  if(nm=="A",     1 - 2^(-s),
     nm=="B",     1 - 6*2^(-s) - 8*4^(-s),
     nm=="C",     1 - 8*2^(-s),
     nm=="D",     1,
     nm=="E",     1 - 8*2^(-s),
     nm=="F",     1 - 7*2^(-s) - 8*4^(-s),
     nm=="alpha", 1 -17*2^(-s) -9*3^(-s) +16*4^(-s) +153*6^(-s) -144*12^(-s),
     nm=="gamma", 1 -28*2^(-s) +63*3^(-s) -36*6^(-s),
     nm=="delta", 1 -14*2^(-s) -3^(-s) +16*4^(-s) +14*6^(-s) -16*12^(-s),
     nm=="eps",   1 -21*2^(-s) +84*4^(-s) -64*8^(-s),
     nm=="zeta",  1,
     nm=="eta",   1 -14*2^(-s) -16*4^(-s),
     error("unknown row ", nm));
}

/* the target L(Phi, w+1), computed from the Dirichlet factorisation        */
Ltarget(nm) = {
  if(nm=="A",     Pmellin(nm,2)*zeta(2)*lfun(-3,0),
     nm=="B",     Pmellin(nm,2)*lfun(-3,2)*zeta(0),
     nm=="C",     Pmellin(nm,2)*lfun(-3,2)*zeta(0),
     nm=="D",     zeta(2)*LnuD(0),
     nm=="E",     Pmellin(nm,2)*lfun(-4,2)*zeta(0),
     nm=="F",     Pmellin(nm,2)*lfun(-3,2)*zeta(0),
     nm=="alpha", Pmellin(nm,3)*zeta(3)*zeta(0),
     nm=="gamma", Pmellin(nm,3)*zeta(3)*zeta(0),
     nm=="delta", Pmellin(nm,3)*zeta(3)*zeta(0),
     nm=="eps",   Pmellin(nm,3)*zeta(3)*zeta(0),
     nm=="zeta",  lfun(-3,3)*lfun(-3,0),
     nm=="eta",   Pmellin(nm,3)*lfun(5,3)*zeta(0),
     error("unknown row ", nm));
}
/* L(nu, 0) for the quartic-character combination of row D (exact = 1/5)   */
LnuD(s) = 1/5;
