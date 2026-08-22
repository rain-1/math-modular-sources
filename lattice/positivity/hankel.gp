/* lattice/positivity/hankel.gp --- Task 1, the decidable criterion.
   A row admits a BEUKERS-TYPE positive-kernel integral representation
        r_n := a_n*xi - b_n = int u^n d mu(u),   mu >= 0,
   iff (r_n) is a Hamburger moment sequence, i.e. iff every Hankel matrix
        H_d = (r_{i+k})_{0<=i,k<=d}
   is positive semidefinite; the measure sits on [0,oo) iff in addition the
   shifted Hankel (r_{i+k+1}) is PSD (Stieltjes), and on a compact [0,L] iff
   also (L*r_{i+k} - r_{i+k+1}) is PSD (Hausdorff).  This is exactly the shape
   of every classical construction: Beukers' zeta(2) kernel is
     int int [x(1-x)y(1-y)/(1-xy)]^n dxdy/(1-xy),
   Beukers' zeta(3) kernel is
     int int int [x(1-x)y(1-y)z(1-z)/(1-(1-xy)z)]^n dxdydz/(1-(1-xy)z),
   and Zudilin's Catalan kernel is
     int int [x(1-x)y(1-y)/(1-xy)]^m dxdy / ( sqrt(x(1-y)) (1-xy) ),
   all of them "u^n against a positive measure".
   Prepend signs.gp for row2/row3.                                          */
\p 2000
{
hank(nm, ord, aa, bb, cc, xi, DMAX) =
 my(NN=2*DMAX+2, r=if(ord==2, row2(aa,bb,cc,NN), row3(aa,bb,cc,NN)), A=r[1], B=r[2]);
 my(rr=vector(NN+1, i, A[i]*xi-B[i]));
 my(dH=-1, dS=-1, okH=1, okS=1, lastH=0, lastS=0);
 for(d=0,DMAX,
  my(H=matrix(d+1,d+1,i,k, rr[i+k-1]),
     SH=matrix(d+1,d+1,i,k, rr[i+k]),
     dh=matdet(H), ds=matdet(SH));
  if(okH, if(dh>0, dH=d; lastH=dh, okH=0));
  if(okS, if(dh>0 && ds>0, dS=d; lastS=ds, okS=0)));
 printf("%-14s Hamburger: all det H_d>0 for d<=%2d%s (log10 det H_last=%+9.2f) | Stieltjes: d<=%2d%s\n",
   nm, dH, if(dH==DMAX," (=DMAX)",""), if(lastH>0,log(lastH)/log(10),0),
   dS, if(dS==DMAX," (=DMAX)",""));
}
