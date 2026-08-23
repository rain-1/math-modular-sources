/* lattice/p2_structure/p2core.gp
   Structural study of the 2-dimensional Catalan congruence lattices.
   Prepend lattice/positivity/rows_pos.gp (for zudrow/nestrow/klat2), and load
   the cached rows with rdrows().  No builtin names shadowed
   (no psi, M, Phi, S, cmp; matrices are MAT/BB, moduli MOD).

   Conventions.  Oriented integer coordinates: with sZ = sign(X_n G - Y_n),
   sN = sign(V_n G - U_n), we substitute (a,b) = (sZ c_Z, sN c_N), so the
   positive cone P is the closed first quadrant {a>=0, b>=0} of Z^2 and the
   weighted metric is diag(lz,ln), lz=|X G-Y|/MOD, ln=|V G-U|/MOD.
   All lattice arithmetic is exact in Z^2; only the two weights are real.   */

default(parisizemax, 8000000000);

/* ---------- row cache -------------------------------------------------- */
{
rdrows(FILE) = my(v=readstr(FILE), h=Map());
  for(i=1,#v, my(w=strsplit(v[i]," "));
    if(#w==5, mapput(h, eval(w[1]), [eval(w[2]),eval(w[3]),eval(w[4]),eval(w[5])])));
  h;
}

/* ---------- the two lattices ------------------------------------------- */
/* honest lattice: both combinations integral after division by MOD */
kfull(XX,YY,VV,UU,MOD) = klat2(XX,YY,VV,UU,MOD);
/* the single-congruence lattice of 06_threshold.tex: c.(Y,U) = 0 mod SS   */
{
kcong(YY,UU,SS) = my(K=matkerint(matrix(1,3,i,j,[YY,UU,-SS][j])));
  mathnf(matrix(2,#K[1,],i,j,K[i,j]));
}

/* ---------- weighted 2-d Gauss (Lagrange) reduction --------------------- */
/* BB: 2x2 integer matrix, columns = basis of the oriented lattice.
   Returns [Br, l1, l2, mrg] with Br reduced (|Br_1| = lambda_1 <= |Br_2| =
   lambda_2 the successive minima -- Gauss's theorem in dimension 2), and
   mrg the smallest rounding margin |mu - round(mu)| met on the way (a
   numerical-safety diagnostic; the reduction is exact whenever mrg is not
   within the working precision of 1/2).                                    */
{
gred(BB, lz, ln) =
 my(u=[BB[1,1],BB[2,1]], v=[BB[1,2],BB[2,2]], mrg=1/2, sw=1);
 my(q2(x) = (x[1]*lz)^2 + (x[2]*ln)^2);
 my(ipr(x,y) = x[1]*y[1]*lz^2 + x[2]*y[2]*ln^2);
 if(q2(u) > q2(v), my(t=u); u=v; v=t);
 while(sw,
   my(rt = ipr(u,v)/q2(u), mu = round(rt), d = abs(rt-mu));
   if(d<mrg, mrg=d);
   v = [v[1]-mu*u[1], v[2]-mu*u[2]];
   if(q2(v) < q2(u), my(t=u); u=v; v=t, sw=0));
 [[u[1],v[1];u[2],v[2]], sqrt(q2(u)), sqrt(q2(v)), mrg];
}

/* ---------- cone geometry ---------------------------------------------- */
/* For x = i*b1 + j*b2 the cone conditions are two linear inequalities in i
   at fixed j; the objective (linear form) and the squared norm are both
   minimised over the resulting integer interval at an endpoint for the
   form, and at an endpoint or the interior rounding for the norm -- we scan
   the two endpoints plus the two interior candidates.  Completeness: any
   x in the cone with |x| <= L has |j| <= L*lambda_1/covol, so JHI below is
   a proof-carrying bound once L is any achieved cone value.               */
{
coneiv(Br, j) =  /* [0] if empty; else [lo,hi] (lo,hi integers or -oo/oo) */
 my(lo=-oo, hi=oo);
 for(r=1,2,
   my(aa=Br[r,1], bb=j*Br[r,2]);
   if(aa==0,
     if(bb<0, return([0])),
     my(t=-bb/aa);
     if(aa>0, if(lo==-oo || t>lo, lo=t), if(hi==oo || t<hi, hi=t))));
 if(lo!=-oo, lo=ceil(lo));
 if(hi!=oo, hi=floor(hi));
 if(lo!=-oo && hi!=oo && lo>hi, return([0]));
 [lo,hi];
}
{
conescan(Br, lz, ln, JHI) =
 my(bestf=0, bfv=0, bestq=0, bqv=0);
 for(j=-JHI, JHI,
  my(iv=coneiv(Br,j));
  if(#iv==2,
   my(lo=iv[1], hi=iv[2], cand=List(), inr(i) = (lo==-oo || i>=lo) && (hi==oo || i<=hi));
   /* the linear form is minimised at an endpoint of the interval, the norm
      at an endpoint or at the rounding of the real minimiser; the extra
      lo+1, hi-1 cover the case where an endpoint is the zero vector. */
   if(lo!=-oo, listput(cand,lo); listput(cand,lo+1));
   if(hi!=oo, listput(cand,hi); listput(cand,hi-1));
   my(d1=Br[1,1]*(j*Br[1,2])*lz^2 + Br[2,1]*(j*Br[2,2])*ln^2,
      d2=(Br[1,1]*lz)^2 + (Br[2,1]*ln)^2, rr=-d1/d2);
   listput(cand,floor(rr)); listput(cand,ceil(rr));
   for(t=1,#cand,
     my(i=cand[t]);
     if(inr(i),
      my(x1=i*Br[1,1]+j*Br[1,2], x2=i*Br[2,1]+j*Br[2,2]);
      if(x1||x2,
        my(fv=x1*lz+x2*ln, qv=sqrt((x1*lz)^2+(x2*ln)^2));
        if(bestf==0 || fv<bestf, bestf=fv; bfv=[x1,x2]);
        if(bestq==0 || qv<bestq, bestq=qv; bqv=[x1,x2]))))));
 [bestf,bfv,bestq,bqv];
}
/* the eight short vectors +-b1,+-b2,+-(b1+b2),+-(b1-b2): one of them lies in
   the closed cone whenever (b1,b2) is Gauss-reduced (Lemma A of
   P2_STRUCTURE.md).  Returns the smallest linear-form value found, or 0.   */
{
eight(Br, lz, ln) =
 my(best=0, cs=[[1,0],[0,1],[1,1],[1,-1]]);
 for(t=1,4, for(s=0,1,
   my(sg=if(s,1,-1), i=sg*cs[t][1], j=sg*cs[t][2],
      x1=i*Br[1,1]+j*Br[1,2], x2=i*Br[2,1]+j*Br[2,2]);
   if((x1>=0)&&(x2>=0)&&(x1||x2),
     my(fv=x1*lz+x2*ln); if(best==0||fv<best, best=fv))));
 best;
}
