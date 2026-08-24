/* lattice/p2_scale/scorel.gp
   Exact core for the large-n extension of the P2' programme.

   Everything that decides anything is done in exact integer arithmetic.
   The single real input, Catalan's constant G, enters only through the
   metric ratio r = lambda_N/lambda_Z, which is bracketed by two exact
   dyadic rationals RN/2^BB <= r <= (RN+1)/2^BB; every reported quantity is
   computed twice, once with each endpoint, and an instance is flagged `ok`
   only if the two runs agree exactly.  Thus no statement below depends on
   floating point.

   Speed.  Three replacements for the n <= 200 pipeline make n <= 1000 possible:
     * h11, h22 by the closed forms of P2_HOLONOMIC.md sec.1.2 and h12 by one
       modular inverse, instead of matkerint (13.5 s at n = 1000);
     * the convergent ladder by its own three-term recurrence
       e_i = a_i e_{i-1} + e_{i-2},  e_i = h12 q_i - h11 p_i,
       streamed twice (bit-lengths first, exact squares only in a window of
       about ten indices around the balance);
     * the Gauss reduction started from the ladder pair (w_i, w_{i-1}), which
       is already reduced by Theorem 1 of P2_STRUCTURE.md -- the reduction
       loop is still run and its verdict recorded (`red` column), so the
       theorem is re-verified rather than assumed.

   Prepend: lattice/positivity/rows_pos.gp, lattice/p2_structure/p2core.gp.
   No builtin names shadowed (no psi, M, Phi, S, cmp).                      */

default(parisize,  3000000000);
default(parisizemax, 9000000000);

BB = 256;                      /* bits of the dyadic bracket for r          */

dlcm(mm) = my(pv=primes([2,mm])); prod(i=1,#pv, pv[i]^logint(mm,pv[i]));

/* ---- Hermite data of the ORIENTED lattice, by the closed forms --------- */
/* Ko = {(u,v) : MOD | u*sZ*X + v*V,  MOD | u*sZ*Y + v*U},  sZ = (-1)^n.
   Returns [h11,h12,h22, g0, gXY, okY, okX, okIdx].                        */
{
hermd(XX,YY,VV,UU,MOD,sZ) =
 my(gXY = gcd(MOD, gcd(XX,YY)), g0 = gcd(gcd(XX,YY), gcd(VV,UU)),
    h11 = MOD/gXY, h22 = gXY/g0,
    okY = (gcd(MOD,YY) == gXY));
 my(h12 = if(okY, lift( Mod(-sZ*(UU/g0), h11) * Mod(YY/gXY, h11)^(-1) ), -1));
 my(cY = if(okY, (sZ*h12*YY + h22*UU) % MOD == 0, 0),
    cX = if(okY, (sZ*h12*XX + h22*VV) % MOD == 0, 0));
 my(HL = mathnf([XX, VV, MOD, 0; YY, UU, 0, MOD]),
    okIdx = (h11*h22 == MOD^2/(HL[1,1]*HL[2,2])));
 [h11, h12, h22, g0, gXY, cY, cX, okIdx];
}

/* ---- exact integer Gauss (Lagrange) reduction for diag(RD^2, RN^2) ----- */
/* u,v are integer 2-vectors; returns [u,v,steps] reduced.                  */
{
gredx(u, v, RD, RN) =
 my(q2(x) = (RD*x[1])^2 + (RN*x[2])^2,
    ipr(x,y) = RD^2*x[1]*y[1] + RN^2*x[2]*y[2], st = 0, sw = 1);
 if(q2(u) > q2(v), my(w=u); u=v; v=w; st++);
 while(sw,
   my(mu = round(ipr(u,v)/q2(u)));
   if(mu, st++);
   v = [v[1]-mu*u[1], v[2]-mu*u[2]];
   if(q2(v) < q2(u), my(w=u); u=v; v=w; st++, sw=0));
 [u, v, st];
}

/* ---- exact cone scan (Proposition B of P2_STRUCTURE.md) ---------------- */
{
conex(b1, b2, RD, RN, JHI) =
 my(bestf = 0, bfv = 0, bestq = 0, bqv = 0,
    Br = [b1[1], b2[1]; b1[2], b2[2]]);
 for(j = -JHI, JHI,
  my(iv = coneiv(Br, j));
  if(#iv == 2,
   my(lo=iv[1], hi=iv[2], cand=List(), inr(i) = (lo==-oo || i>=lo) && (hi==oo || i<=hi));
   if(lo!=-oo, listput(cand,lo); listput(cand,lo+1));
   if(hi!=oo, listput(cand,hi); listput(cand,hi-1));
   my(d1 = RD^2*b1[1]*(j*b2[1]) + RN^2*b1[2]*(j*b2[2]),
      d2 = (RD*b1[1])^2 + (RN*b1[2])^2, rr = -d1/d2);
   listput(cand,floor(rr)); listput(cand,ceil(rr));
   for(tt=1,#cand,
     my(i=cand[tt]);
     if(inr(i),
      my(x1 = i*b1[1]+j*b2[1], x2 = i*b1[2]+j*b2[2]);
      if(x1||x2,
        my(fv = RD*x1 + RN*x2, qv = (RD*x1)^2 + (RN*x2)^2);
        if(bestf==0 || fv<bestf, bestf=fv; bfv=[x1,x2]);
        if(bestq==0 || qv<bestq, bestq=qv; bqv=[x1,x2]))))));
 [bestf, bfv, bestq, bqv];
}
/* smallest linear-form value among the eight Lemma-A vectors in the cone */
{
eightx(b1, b2, RD, RN) =
 my(best = 0, cs = [[1,0],[0,1],[1,1],[1,-1]]);
 for(tt=1,4, for(s=0,1,
   my(sg = if(s,1,-1), i = sg*cs[tt][1], j = sg*cs[tt][2],
      x1 = i*b1[1]+j*b2[1], x2 = i*b1[2]+j*b2[2]);
   if((x1>=0) && (x2>=0) && (x1||x2),
     my(fv = RD*x1 + RN*x2); if(best==0 || fv<best, best=fv))));
 best;
}

/* ---- the balance index on the convergent ladder ------------------------ */
/* cf: the continued fraction of h12/h11 (cf[1] = 0).  Two streamed passes:
   the first records only bit-length proxies, the second evaluates the exact
   integer objective RD^2 e_i^2 + RN^2 (q_i h22)^2 on the (short) window in
   which the proxy is within 12 bits of its minimum -- a proof-carrying
   window, since the objective lies within 2^{+-5} of its proxy.
   Returns [idx0, e_i, q_i, e_{i-1}, q_{i-1}, LL], idx0 0-based.           */
{
balx(cf, h11, h12, h22, RD, RN) =
 my(LL = #cf, lgh = exponent(h22), lgr = exponent(RN) - exponent(RD),
    e0, e1, q0, q1, pr = vector(LL), pmin, lo, hi,
    be = 0, bv = 0, bE = 0, bQ = 0, pE = 0, pQ = 0, vv);
 e0 = h12; q0 = 1;
 e1 = if(LL>=2, h12*cf[2] - h11, 0); q1 = if(LL>=2, cf[2], 1);
 pr[1] = max(if(e0==0,-10^9,2*exponent(e0)), 2*(exponent(q0)+lgh+lgr));
 if(LL>=2, pr[2] = max(if(e1==0,-10^9,2*exponent(e1)), 2*(exponent(q1)+lgh+lgr)));
 for(i=3,LL,
   my(ee = cf[i]*e1 + e0, qq = cf[i]*q1 + q0);
   e0=e1; e1=ee; q0=q1; q1=qq;
   pr[i] = max(if(e1==0,-10^9,2*exponent(e1)), 2*(exponent(q1)+lgh+lgr)));
 pmin = vecmin(pr); lo = LL; hi = 1;
 for(i=1,LL, if(pr[i] <= pmin+12, if(i<lo, lo=i); if(i>hi, hi=i)));
 lo = max(lo-1, 1); hi = min(hi+1, LL);
 /* second pass */
 e0 = h12; q0 = 1;
 e1 = if(LL>=2, h12*cf[2] - h11, 0); q1 = if(LL>=2, cf[2], 1);
 if(lo<=1,
   vv = RD^2*h12^2 + RN^2*h22^2;
   be = 1; bv = vv; bE = h12; bQ = 1; pE = 0; pQ = 0);
 if(LL>=2 && lo<=2 && hi>=2,
   vv = RD^2*e1^2 + RN^2*(q1*h22)^2;
   if(be==0 || vv<bv, be=2; bv=vv; bE=e1; bQ=q1; pE=h12; pQ=1));
 for(i=3,LL,
   my(ee = cf[i]*e1 + e0, qq = cf[i]*q1 + q0);
   e0=e1; e1=ee; q0=q1; q1=qq;
   if(i>=lo && i<=hi,
     vv = RD^2*e1^2 + RN^2*(q1*h22)^2;
     if(be==0 || vv<bv, be=i; bv=vv; bE=e1; bQ=q1; pE=e0; pQ=q0)));
 [be-1, bE, bQ, pE, pQ, LL];
}
