/* lib.gp -- Cooper's three rows: eta quotients, F, x, Phi, c(m). */

/* rows: [name, N, B, C, lam1, lam2, a, b, c, d]  (recurrence R3) */
{ROWS = [
 ["s7",   7, 13, 49, 27, -1, 13, 4, -27,  3],
 ["s10", 10,  6, 25, 16, -4,  6, 2, -64,  4],
 ["s18", 18, 14,  1, 16, 12, 14, 6, 192,-12]
];}

/* prod_{n>=1} (1 - q^{d n}) to O(q^N) */
{ Ed(d,N) = my(s = 1 + O('q^N)); for(n=1, N\d, s = s*(1 - 'q^(d*n))); s; }

/* the three Hauptmoduln u = q + ... */
{ Umod(k,N) = my(q='q + O('q^N));
  if(k==1, return( q*(Ed(7,N)/Ed(1,N))^4 ));
  if(k==2, return( q*(Ed(5,N)*Ed(10,N)/(Ed(1,N)*Ed(2,N)))^2 ));
  q*(Ed(2,N)*Ed(3,N)^2*Ed(18,N)/(Ed(1,N)*Ed(6,N)^2*Ed(9,N)))^6;
}

Dop(f) = 'q*deriv(f,'q);

/* row data as a vector [u,F,x,Phi] */
{ Setup(k,N) = my(R=ROWS[k], B=R[3], C=R[4], u, F, x, Ph);
  u = Umod(k,N);
  F = Dop(u)/u;
  x = u/(1 + B*u + C*u^2);
  Ph = F*Dop(x);
  [u,F,x,Ph];
}

/* alternative Phi = u(1-Cu^2)F^2/(1+Bu+Cu^2)^2 */
{ Phi2(k,N) = my(R=ROWS[k], B=R[3], C=R[4], u, F);
  u = Umod(k,N); F = Dop(u)/u;
  u*(1-C*u^2)*F^2/(1+B*u+C*u^2)^2;
}

/* recurrence rows */
{ genrow(k,N) = my(R=ROWS[k], a=R[7], b=R[8], c=R[9], d=R[10], A=vector(N+1), Bv=vector(N+1), P, Q);
  A[1]=1; A[2]=b; Bv[1]=0; Bv[2]=1;
  for(n=1,N-1, P=(2*n+1)*(a*n^2+a*n+b); Q=n*(c*n^2+d);
    A[n+2]=(P*A[n+1]-Q*A[n])/(n+1)^3;
    Bv[n+2]=(P*Bv[n+1]-Q*Bv[n])/(n+1)^3);
  [A,Bv];
}
