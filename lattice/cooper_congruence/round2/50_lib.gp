/* 50_lib.gp -- fast version of lib.gp (pentagonal eta via PARI's eta(series)).
   Same conventions as lib.gp; verified against it in 50_sanity.log.        */

{ROWS = [
 ["s7",   7, 13, 49, 27, -1, 13, 4, -27,  3],
 ["s10", 10,  6, 25, 16, -4,  6, 2, -64,  4],
 ["s18", 18, 14,  1, 16, 12, 14, 6, 192,-12]
];}
NAM = ["s7","s10","s18"];
LEV = [7,10,18];

/* prod_{n>=1}(1-q^{dn}) + O(q^M)  -- PARI's eta on a power series */
{ Ed(d,M) = eta('q^d + O('q^M)); }

{ Umod(k,M) = my(q='q + O('q^M));
  if(k==1, return( q*(Ed(7,M)/Ed(1,M))^4 ));
  if(k==2, return( q*(Ed(5,M)*Ed(10,M)/(Ed(1,M)*Ed(2,M)))^2 ));
  q*(Ed(2,M)*Ed(3,M)^2*Ed(18,M)/(Ed(1,M)*Ed(6,M)^2*Ed(9,M)))^6;
}

Dop(f) = 'q*deriv(f,'q);

{ Setup(k,M) = my(R=ROWS[k], B=R[3], C=R[4], u, F, x, Ph);
  u = Umod(k,M);
  F = Dop(u)/u;
  x = u/(1 + B*u + C*u^2);
  Ph = F*Dop(x);
  [u,F,x,Ph];
}

psival(k,p) = if(k<3, 1, kronecker(-3,p));
psin(k,n) = if(k<3, 1, kronecker(-3,n));

{ Cvec(k,M) = my(S=Setup(k,M+4)); vector(M,m,polcoeff(S[4],m)); }
{ CPvec(k,M) = my(v=Cvec(k,M)); vector(M,m,v[m]/m); }

{ Bvec(k,cp) = my(M=#cp, b=vector(M));
  for(n=1,M, my(s=0); fordiv(n, d, s += moebius(d)*psin(k,d)*cp[n/d]); b[n]=s);
  b; }

{ genrow(k,NMAX) = my(R=ROWS[k], a=R[7], b=R[8], c=R[9], d=R[10], A=vector(NMAX+1), Bv=vector(NMAX+1), P, Q);
  A[1]=1; A[2]=b; Bv[1]=0; Bv[2]=1;
  for(n=1,NMAX-1, P=(2*n+1)*(a*n^2+a*n+b); Q=n*(c*n^2+d);
    A[n+2]=(P*A[n+1]-Q*A[n])/(n+1)^3;
    Bv[n+2]=(P*Bv[n+1]-Q*Bv[n])/(n+1)^3);
  [A,Bv];
}
