\\ 15_jaclib.gp -- weak Jacobi forms as bivariate objects.
\\ A Jacobi form is stored as [R, P] meaning  z^(-R) * P(z),  where P is a t_POL in 'x whose
\\ coefficients are t_SER in 'y with Q^2 = q.  All series are truncated at O(Q^(2*NQ)).
\\ Generators: phi_{-2,1}, phi_{0,1}, phi_{-1,2}; plus E4, E6, Delta.

NQ = 140;                  \\ q-precision
PRECQ = 2*NQ + 4;         \\ Q-precision

qq = 'y^2 + O('y^PRECQ);
one = 1 + O('y^PRECQ);

{ jmul(A,B) = [A[1]+B[1], A[2]*B[2]]; }
{ jscal(c,A) = [A[1], c*A[2]]; }
{ jadd(A,B) = my(R=max(A[1],B[1])); [R, A[2]*'x^(R-A[1]) + B[2]*'x^(R-B[1])]; }
{ jsub(A,B) = my(R=max(A[1],B[1])); [R, A[2]*'x^(R-A[1]) - B[2]*'x^(R-B[1])]; }
{ jpow(A,n) = my(r=[0, one]); for(i=1,n, r=jmul(r,A)); r; }
\\ coefficient of z^r q^n
{ jcoef(A,n,r) = my(P=A[2], e=r+A[1], c);
  if(e<0 || e>poldegree(P,'x), return(0));
  c = polcoeff(P, e, 'x);
  if(type(c)=="t_INT" || type(c)=="t_FRAC", return(if(n==0, c, 0)));
  polcoeff(c, 2*n, 'y); }

\\ ---- building blocks -------------------------------------------------------
\\ prod_{n=1..NQ} (1 - q^(n) z^a)  as [R,P]
{ prodfac(a, sgn, halfshift) =
   my(r=[0,one], t);
   for(n=1, NQ+1,
     my(ex = if(halfshift, 2*n-1, 2*n));
     if(ex > PRECQ-2, break);
     t = if(a>0, [0, one - sgn*'y^ex*'x^a], [-a, 'x^(-a)*one - sgn*'y^ex*one]);
     r = jmul(r,t));
   r; }
\\ eta-type scalar products
{ scalprod(sgn, halfshift, e) = my(r=one, ex);
   for(n=1,NQ+1, ex = if(halfshift, 2*n-1, 2*n); if(ex>PRECQ-2, break);
       r = r*(1 - sgn*'y^ex)^e);
   r; }

\\ phi_{-2,1} = (z - 2 + 1/z) * prod (1-q^n z)^2 (1-q^n/z)^2 / (1-q^n)^4
{ mkphi21() = my(A, B, D);
  A = jpow(prodfac(1,1,0), 2);
  B = jpow(prodfac(-1,1,0), 2);
  D = scalprod(1,0,4);
  jscal(1/D, jmul([1, 'x^2 - 2*'x + 1], jmul(A,B))); }

\\ phi_{-1,2} = (z - 1/z) * prod (1-q^n z^2)(1-q^n z^-2)/(1-q^n)^2
{ mkphi12() = my(A,B,D);
  A = prodfac(2,1,0); B = prodfac(-2,1,0); D = scalprod(1,0,2);
  jscal(1/D, jmul([1, 'x^2 - 1], jmul(A,B))); }

\\ phi_{0,1} = 4 [ (th2(z)/th2)^2 + (th3(z)/th3)^2 + (th4(z)/th4)^2 ]
{ mkphi01() = my(T2,T3,T4,d2,d3,d4);
  T2 = jmul([1,'x^2 + 2*'x + 1], jmul(jpow(prodfac(1,-1,0),2), jpow(prodfac(-1,-1,0),2)));
  d2 = 4*scalprod(-1,0,4);
  T3 = jmul(jpow(prodfac(1,-1,1),2), jpow(prodfac(-1,-1,1),2));
  d3 = scalprod(-1,1,4);
  T4 = jmul(jpow(prodfac(1,1,1),2), jpow(prodfac(-1,1,1),2));
  d4 = scalprod(1,1,4);
  jscal(4, jadd(jscal(1/d2,T2), jadd(jscal(1/d3,T3), jscal(1/d4,T4)))); }

\\ ---- level-one modular forms as scalars in Q ------------------------------
{ mkE(k) = my(s = one, c = -2*k/bernfrac(k));
   for(n=1,NQ, s += c*sigma(n,k-1)*'y^(2*n)); s; }
E4S = mkE(4); E6S = mkE(6);
DELS = ('y^2*eta('y^2 + O('y^PRECQ))^24 + O('y^PRECQ));
