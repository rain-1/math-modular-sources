/* ================================================================
   rows.gp -- Apery rows of the census, exact rational sequences.
   Returns [A, B] with A[n+1]=a_n, B[n+1]=b_n, n=0..N.
   ================================================================ */

/* (R2)  (n+1)^2 u_{n+1} = (a n^2 + a n + b) u_n - c n^2 u_{n-1}   */
rowR2(a,b,c,N) =
{ my(A=vector(N+1), B=vector(N+1));
  A[1]=1; A[2]=b; B[1]=0; B[2]=1;
  for(n=1,N-1,
    A[n+2] = ((a*n^2+a*n+b)*A[n+1] - c*n^2*A[n])/(n+1)^2;
    B[n+2] = ((a*n^2+a*n+b)*B[n+1] - c*n^2*B[n])/(n+1)^2);
  [A,B];
}

/* (R3)  (n+1)^3 u_{n+1} = (2n+1)(a n^2 + a n + b) u_n - c n^3 u_{n-1} */
rowR3(a,b,c,N) =
{ my(A=vector(N+1), B=vector(N+1));
  A[1]=1; A[2]=b; B[1]=0; B[2]=1;
  for(n=1,N-1,
    A[n+2] = ((2*n+1)*(a*n^2+a*n+b)*A[n+1] - c*n^3*A[n])/(n+1)^3;
    B[n+2] = ((2*n+1)*(a*n^2+a*n+b)*B[n+1] - c*n^3*B[n])/(n+1)^3);
  [A,B];
}

/* Cooper s18:  (n+1)^3 u_{n+1} = 2(2n+1)(7n^2+7n+3) u_n - 12 n(16n^2-1) u_{n-1}
   A_0=1, A_1=6 (SLOPE_CENSUS.md) */
rowS18(N) =
{ my(A=vector(N+1), B=vector(N+1));
  A[1]=1; A[2]=6; B[1]=0; B[2]=1;
  for(n=1,N-1,
    A[n+2] = (2*(2*n+1)*(7*n^2+7*n+3)*A[n+1] - 12*n*(16*n^2-1)*A[n])/(n+1)^3;
    B[n+2] = (2*(2*n+1)*(7*n^2+7*n+3)*B[n+1] - 12*n*(16*n^2-1)*B[n])/(n+1)^3);
  [A,B];
}

/* weight-3 cusp row, f = eta(2t)^3 eta(6t)^3, level 12, L(f,2):
   (n+1)^2 u_{n+1} = (20n^2+10n+2) u_n - 16(2n-1)^2 u_{n-1},  a_1 = 2  */
rowCusp(N) =
{ my(A=vector(N+1), B=vector(N+1));
  A[1]=1; A[2]=2; B[1]=0; B[2]=1;
  for(n=1,N-1,
    A[n+2] = ((20*n^2+10*n+2)*A[n+1] - 16*(2*n-1)^2*A[n])/(n+1)^2;
    B[n+2] = ((20*n^2+10*n+2)*B[n+1] - 16*(2*n-1)^2*B[n])/(n+1)^2);
  [A,B];
}

/* Domb cuspidal apparatus (main6_fully_completed_v2.tex, eq. (starrec)):
   A = Domb (10,4,64);  B*_0=0, B*_1=1,
   m^3 B*_m = (2m-1)(10(m-1)^2+10(m-1)+4) B*_{m-1} - 64 (m-1)^3 B*_{m-2}
              + binomial(2m-2,m-1)                                        */
rowDombStar(N) =
{ my(A, Bs=vector(N+1));
  A = rowR3(10,4,64,N)[1];
  Bs[1]=0; Bs[2]=1;
  for(m=2,N,
    Bs[m+1] = ((2*m-1)*(10*(m-1)^2+10*(m-1)+4)*Bs[m] - 64*(m-1)^3*Bs[m-1]
               + binomial(2*m-2,m-1))/m^3);
  [A,Bs];
}

/* Zudilin's Catalan row (lattice/catalan_positivity/rows_common.gp):
   Q_0=1, Q_1=7/4, P_0=0, P_1=13/8, and for m>=1
   a Q_{m+1} = b Q_m + c Q_{m-1} with
   a=(2m+1)^2(2m+2)^2(20m^2-8m+1)
   b=3520m^6+5632m^5+2064m^4-384m^3-156m^2+16m+7
   c=(2m-1)^2(2m)^2(20m^2+32m+13)                                        */
rowZud(M) =
{ my(Q=vector(M+1), P=vector(M+1), a,b,c);
  Q[1]=1; Q[2]=7/4; P[1]=0; P[2]=13/8;
  for(m=1,M-1,
    a=(2*m+1)^2*(2*m+2)^2*(20*m^2-8*m+1);
    b=3520*m^6+5632*m^5+2064*m^4-384*m^3-156*m^2+16*m+7;
    c=(2*m-1)^2*(2*m)^2*(20*m^2+32*m+13);
    Q[m+2]=(b*Q[m+1]+c*Q[m])/a;
    P[m+2]=(b*P[m+1]+c*P[m])/a);
  [Q,P];
}
