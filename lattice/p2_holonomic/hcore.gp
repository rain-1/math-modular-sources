/* lattice/p2_holonomic/hcore.gp
   Explicit arithmetic of the Hermite data (h11,h12,h22) of the Catalan
   congruence lattices K_n, and of the 2-adic objects behind them.
   Prepend lattice/positivity/rows_pos.gp and lattice/p2_structure/p2core.gp.
   No builtin names shadowed (no psi, M, Phi, S, cmp).                       */

default(parisizemax, 8000000000);

/* ---------- oriented Hermite data ------------------------------------- */
/* returns [h11,h12,h22, MOD, idx, sZ]  for the honest lattice at (n,kk).
   GG only enters through sZ = sign(X GG - Y), which is (-1)^n by Zudilin
   positivity; pass SZ=0 to use the theorem instead of GG.                  */
{
hdat(XX,YY,VV,UU,SS,kk,nn,SZ) =
 my(TT=2^floor(kk*nn), MOD=SS*TT, B0=kfull(XX,YY,VV,UU,MOD),
    sZ=if(SZ, SZ, (-1)^nn),
    BB=[sZ*B0[1,1], sZ*B0[1,2]; B0[2,1], B0[2,2]],
    HH=mathnf(BB));
 [HH[1,1], HH[1,2], HH[2,2], MOD, abs(matdet(B0)), sZ];
}

/* closed-form predictions of P2_HOLONOMIC.md sec.1 */
{
hpred(XX,YY,VV,UU,MOD) =
 my(g0=gcd([XX,YY,VV,UU]), gXY=gcd(MOD,gcd(XX,YY)));
 [MOD/gXY, gXY/g0, g0];        /* [h11pred, h22pred, g0] */
}

/* the two congruences that pin h12: sZ*h12*Y = -h22*U  and
   sZ*h12*X = -h22*V, both mod MOD.  Returns [ok_Y, ok_X, mY, mX] where
   mY = MOD/gcd(MOD,Y), mX = MOD/gcd(MOD,X) are the two moduli whose lcm
   must be h11.                                                            */
{
hchk12(XX,YY,VV,UU,MOD,h11,h12,h22,sZ) =
 my(okY = ((sZ*h12*YY + h22*UU) % MOD == 0),
    okX = ((sZ*h12*XX + h22*VV) % MOD == 0),
    mY  = MOD/gcd(MOD,YY), mX = MOD/gcd(MOD,XX));
 [okY, okX, mY, mX, lcm(mY,mX)==h11];
}

/* ---------- convergent ladder of theta = h12/h11 ----------------------- */
/* Returns [AA,BB2,CF] with AA[i] = |h12 q_{i-1} - h11 p_{i-1}| (the scaled
   first coordinate) and BB2[i] = q_{i-1}*h22 (the second coordinate), i.e.
   the i-th lattice vector on the convergent ladder is (AA[i], BB2[i]) up to
   sign.  Index convention: entry i corresponds to convergent index i-1, as
   in p2_structure/p2run.gp `cfclass`.                                      */
{
ladder(h11,h12,h22) =
 my(cf=contfrac(h12/h11), LL=#cf, qs=vector(LL), ps=vector(LL));
 qs[1]=1; ps[1]=0;
 if(LL>=2, qs[2]=cf[2]; ps[2]=1);
 for(i=3,LL, qs[i]=cf[i]*qs[i-1]+qs[i-2]; ps[i]=cf[i]*ps[i-1]+ps[i-2]);
 [vector(LL,i,abs(h12*qs[i]-h11*ps[i])), vector(LL,i,qs[i]*h22), cf, qs, ps];
}

/* argmin over the ladder of AA^2 + t*BB2^2, with the exact open interval
   (tlo,thi) of t on which that argmin is attained.  t is a positive real.
   Returns [idx0, tlo, thi] with idx0 the convergent index (0-based).       */
{
balance(AA, BB2, t) =
 my(LL=#AA, be=1, bv=AA[1]^2 + t*BB2[1]^2);
 for(i=2,LL, my(v=AA[i]^2 + t*BB2[i]^2); if(v<bv, bv=v; be=i));
 my(tlo=0, thi=oo);
 for(j=1,LL, if(j!=be,
   my(dA=AA[be]^2-AA[j]^2, dB=BB2[j]^2-BB2[be]^2);
   if(dB>0, my(tt=dA/dB); if(tt>tlo, tlo=tt),
     if(dB<0, my(tt=dA/dB); if(thi==oo || tt<thi, thi=tt)))));
 [be-1, tlo, thi];
}
