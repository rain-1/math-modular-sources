/* 40_core.gp -- shared utilities for the x-side attack (round 2).
   Provides:  AVEC(k,N)  = vector A_0..A_N  with A_n = [x^n] F
              ETAS(k,N)  = the series a(x) = l(x)/(x sqrt(P) F) + O(x^(N+1))
              LSER, FSER, PSER, SQP
   Everything exact over Q.                                              */
read("lib.gp");

{ AVEC(k,N) = genrow(k,N)[1]; }      /* A[n+1] = A_n, n=0..N */

{ FSER(k,N) = Ser(AVEC(k,N),'x); }

{ LSER(k,N) = my(A=AVEC(k,N)); 'x*Ser(vector(N+1,n,A[n]/n),'x); }

{ PSER(k,N) = my(R=ROWS[k],B=R[3],C=R[4]); 1 - 2*B*'x + (B^2-4*C)*'x^2 + O('x^(N+1)); }

{ ETAS(k,N) = my(l=LSER(k,N), F=FSER(k,N), sq=sqrt(PSER(k,N))); l/('x*sq*F); }

/* a_j, j=0..N */
{ AJ(k,N) = my(e=ETAS(k,N)); vector(N+1,j,polcoeff(e,j-1)); }
