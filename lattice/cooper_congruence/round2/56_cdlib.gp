/* 56_cdlib.gp -- shared core for the full-discriminant traces c(d), all three rows.
   Row k: CM discriminant D0k = -3, -4, -36 ; the Heegner family is disc = D0k*d with
   d = 0,1 mod 4 and beta^2 = D0k*d mod 4N solvable.
   Uses heegmin2 (corrected completion search + MINIMAL A) from 59_final3.gp.        */
read("50_lib.gp"); read("heeg.gp"); read("maass2.gp");
D0V   = [-3, -4, -36];      /* the CM discriminant of the pole */
CHIV  = [-3, -4, -3];       /* the genus character D_0 */
SQB   = [5, 6, 18];         /* beta = SQB*m mod 2N on squares d = m^2 */
cont(Q) = gcd(gcd(Q[1],Q[2]),Q[3]);

{ heegmin2(Q,N,beta,PB) = my(a=Q[1],b=Q[2],c=Q[3],dd=b^2-4*a*c,A,B,B0,C,g,s,q,e,t,best=0);
  for(r=0,PB, for(p=-PB,PB,
    if(gcd(p,r)!=1, next);
    A = a*p^2+b*p*r+c*r^2;
    if(A<=0 || A%N!=0, next);
    if(best!=0 && A>=best[1], next);
    g = bezout(p,-r); if(g[3]!=1, next);
    s = g[1]; q = g[2];
    B0 = 2*a*p*q + b*(p*s+q*r) + 2*c*r*s;
    e = gcd(2*A, 2*N);
    if((beta-B0)%e != 0, next);
    t = lift(Mod((beta-B0)/e, 2*N/e) / Mod(2*A/e, 2*N/e));
    B = B0 + 2*A*t; C = (B^2-dd)/(4*A);
    best = [A,B,C]));
  best;
}
{ fhx(k,t) = if(abs(Fmod(k,t)) < 1e-40, fhatQ(t,LEV[k],FSER[k],#Vec(FSER[k][1])-3), fhatC(k,t)[2]); }

/* admissible beta classes mod 2N for discriminant D */
{ betas(D,N) = my(L=List()); for(b=0,2*N-1, if((b^2-D)%(4*N)==0, listput(L,b))); Vec(L); }

/* is d admissible for row k ?  returns the chosen beta, or -1 */
{ chosenbeta(k,d) = my(N=LEV[k], D=D0V[k]*d, bs);
  if(D%4!=0 && (D-1)%4!=0, return(-1));
  bs = betas(D,N);
  if(#bs==0, return(-1));
  if(issquare(d), my(b=(SQB[k]*sqrtint(d))%(2*N)); if(setsearch(Set(bs),b), return(b)));
  bs[1];
}

/* the raw trace T(d) for row k; returns "NOREP" on a missing representative */
{ Traw(k,d,beta,PB) = my(N=LEV[k], D=D0V[k]*d, RF=redforms(D), t=0., rep, ch, om, al, sq=sqrt(-D), nf=0);
  for(i=1,#RF,
    if(k==3 && cont(RF[i])%3==0, next);          /* no Heegner rep exists (proved) */
    rep = heegmin2(RF[i],N,beta,PB);
    if(rep==0, rep = heegmin2(RF[i],N,beta,3*PB));
    if(rep==0, return("NOREP"));
    nf++;
    om = omeg(rep);
    ch = genchar(rep,CHIV[k]);
    if(k==3, ch *= kronecker(-4,cont(RF[i])));   /* the content correction, -36 not fundamental */
    al = (-rep[2] + I*sq)/(2*rep[1]);
    t += ch*fhx(k,al)/om);
  [t,nf,#RF];
}
