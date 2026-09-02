/* lib.gp -- the 15 sporadic rows, exact rational (a_n, b_n).
   R2:  (n+1)^2 u_{n+1} = (a n^2 + a n + b) u_n - c n^2 u_{n-1}
   R3:  (n+1)^3 u_{n+1} = (2n+1)(a n^2 + a n + b) u_n - n(c n^2 + d) u_{n-1}
   a_0=1, a_1=b ; b_0=0, b_1=1.
   ROW = [name, r, a, b, c, d]  (d unused when r=2)                        */

{ROWS = [
 ["A",     2,  7, 2,  -8,   0],
 ["B",     2,  9, 3,  27,   0],
 ["C",     2, 10, 3,   9,   0],
 ["D",     2, 11, 3,  -1,   0],
 ["E",     2, 12, 4,  32,   0],
 ["F",     2, 17, 6,  72,   0],
 ["alpha", 3, 10, 4,  64,   0],
 ["gamma", 3, 17, 5,   1,   0],
 ["delta", 3,  7, 3,  81,   0],
 ["eps",   3, 12, 4,  16,   0],
 ["zeta",  3,  9, 3, -27,   0],
 ["eta",   3, 11, 5, 125,   0],
 ["s7",    3, 13, 4, -27,   3],
 ["s10",   3,  6, 2, -64,   4],
 ["s18",   3, 14, 6, 192, -12]
];}

Pcoef(r,a,b,n) = if(r==2, a*n^2+a*n+b, (2*n+1)*(a*n^2+a*n+b));
Qcoef(r,c,d,n) = if(r==2, c*n^2, n*(c*n^2+d));

/* returns [A,B] with A[n+1]=a_n, B[n+1]=b_n as exact rationals, n=0..N */
{ genrow(R,N) =
  my(r=R[2], a=R[3], b=R[4], c=R[5], d=R[6], A=vector(N+1), B=vector(N+1), pw, P, Q);
  A[1]=1; A[2]=b; B[1]=0; B[2]=1;
  for(n=1,N-1,
    P = Pcoef(r,a,b,n); Q = Qcoef(r,c,d,n); pw = (n+1)^r;
    A[n+2] = (P*A[n+1] - Q*A[n])/pw;
    B[n+2] = (P*B[n+1] - Q*B[n])/pw);
  [A,B];
}

/* integer companion numerators C_n = b_n * (n!)^r  (all-integer recursion) */
{ gennum(R,N) =
  my(r=R[2], a=R[3], b=R[4], c=R[5], d=R[6], C=vector(N+1), P, Q);
  C[1]=0; C[2]=1;
  for(n=1,N-1,
    P = Pcoef(r,a,b,n); Q = Qcoef(r,c,d,n);
    C[n+2] = P*C[n+1] - Q*n^r*C[n]);
  C;
}

badprimes(R) = if(R[5]==0, [], factor(abs(R[5]))[,1]~);
