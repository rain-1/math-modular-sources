/* 17_arith_lib.gp -- the six new mixed-exponent four-term Catalan/L rows of
   FOUR_TERM_DEEP.md §6.4, plus Zagier E, in the MULTI_PRIME_LATTICE.md §1 format.
   Exact rational arithmetic.  Names avoid PARI builtins. */

/* four-term row:  (n+1)^2 u_{n+1} = PP(n) u_n - QQ(n) u_{n-1} + RR(n) u_{n-2}
   a: a_0=1, a_{-1}=a_{-2}=0 ;  b: b_0=0, b_1=1, recurrence from n=1.
   cf[] = [p2,p1,p0, q2,q1,q0, r2,r1,r0]                                */
row4(cf, NN) = {
  my(av = vector(NN+3), bv = vector(NN+3), pp, qq, rr, ix);
  /* index shift: av[n+3] = a_n, so a_{-2}=av[1], a_{-1}=av[2], a_0=av[3] */
  av[1]=0; av[2]=0; av[3]=1;
  bv[1]=0; bv[2]=0; bv[3]=0; 
  for(n = 0, NN-1,
    pp = cf[1]*n^2+cf[2]*n+cf[3];
    qq = cf[4]*n^2+cf[5]*n+cf[6];
    rr = cf[7]*n^2+cf[8]*n+cf[9];
    av[n+4] = (pp*av[n+3] - qq*av[n+2] + rr*av[n+1])/(n+1)^2;
    if(n == 0, bv[4] = 1,
       bv[n+4] = (pp*bv[n+3] - qq*bv[n+2] + rr*bv[n+1])/(n+1)^2));
  [vector(NN+1,i,av[i+2]), vector(NN+1,i,bv[i+2])];
};

/* three-term Zagier row  (n+1)^2 u_{n+1} = (a n^2+a n+b) u_n - c n^2 u_{n-1} */
zrow2(aa, bb, cc, NN) = {
  my(av = vector(NN+1), bv = vector(NN+1));
  av[1]=1; av[2]=bb; bv[1]=0; bv[2]=1;
  for(n = 1, NN-1,
    av[n+2] = ((aa*n^2+aa*n+bb)*av[n+1] - cc*n^2*av[n])/(n+1)^2;
    bv[n+2] = ((aa*n^2+aa*n+bb)*bv[n+1] - cc*n^2*bv[n])/(n+1)^2);
  [av, bv];
};

/* the seven rows: [name, kind, data, g, xi-tag] */
{ROWS = [
 ["row1 (16,20,8; 48,16,0; -128,0,0)",  [16,20,8, 48,16,0, -128,0,0]],
 ["row2 (14,20,8; 28,16,4; 8,0,0)",     [14,20,8, 28,16,4, 8,0,0]],
 ["row3 (6,10,4; -32,-24,-8; 32,0,0)",  [6,10,4, -32,-24,-8, 32,0,0]],
 ["row4 (16,20,8; 68,36,8; 32,0,0)",    [16,20,8, 68,36,8, 32,0,0]],
 ["row5 (17,25,10; 32,24,8; 16,0,0)",   [17,25,10, 32,24,8, 16,0,0]],
 ["row6 (13,20,8; -13,-6,-1; -1,0,0)",  [13,20,8, -13,-6,-1, -1,0,0]]
];}
