/* =====================================================================
   towers.gp -- good-prime tower limits Lambda_a and A-row tower ratios
   u_a for the sporadic Apery rows, computed with an exact O(N) integer
   pass (no divisions at all).

   Scaling trick.  For (R2)  (n+1)^2 u_{n+1} = (a n^2+a n+b) u_n - c n^2 u_{n-1}
   put  uh_n = (n!)^2 u_n ; then
        uh_{n+1} = (a n^2+a n+b) uh_n - c n^4 uh_{n-1} .
   For (R3)  (n+1)^3 u_{n+1} = (2n+1)(a n^2+a n+b) u_n - n(c n^2+d) u_{n-1}
   put  uh_n = (n!)^3 u_n ; then
        uh_{n+1} = (2n+1)(a n^2+a n+b) uh_n - (c n^2+d) n^4 uh_{n-1} .
   Both are integer recurrences; B_n/A_n = Bh_n/Ah_n exactly, so the tower
   ratios need no factorials.  For the A-row ratios A(n')/A(n) we carry the
   p-unit part of n! alongside.

   Everything is done modulo p^K0 with K0 = (w) * v_p(N!) + K + margin,
   w = 2 or 3 the exponent of (n!) above.
   ===================================================================== */

default(parisizemax, 8000000000);

/* row table: [type, a, b, c, d, w, chiD]  chiD = discriminant of chi (0 = trivial) */
ROWS = Map();
mapput(ROWS, "C",     [2, 10,  3,  9, 0, 2, -3]);   /* sum C(n,k)^2 C(2k,k)   */
mapput(ROWS, "D",     [2, 11,  3, -1, 0, 2,  0]);   /* Apery zeta(2)          */
mapput(ROWS, "A",     [2,  7,  2, -8, 0, 2,  0]);   /* Franel                 */
mapput(ROWS, "E",     [2, 12,  4, 32, 0, 2, -4]);
mapput(ROWS, "F",     [2, 17,  6, 72, 0, 2, -3]);
mapput(ROWS, "B",     [2,  9,  3, 27, 0, 2, -3]);
mapput(ROWS, "delta", [3,  7,  3, 81, 0, 3,  0]);
mapput(ROWS, "eta",   [3, 11,  5,125, 0, 3,  5]);
mapput(ROWS, "G",     [3, 17,  5,  1, 0, 3,  0]);   /* Apery zeta(3) (gamma)  */
mapput(ROWS, "alpha", [3, 10,  4, 64, 0, 3,  0]);   /* Domb                   */
mapput(ROWS, "eps",   [3, 12,  4, 16, 0, 3,  0]);
mapput(ROWS, "zeta",  [3,  9,  3,-27, 0, 3, -3]);

chip(chiD, p) = if(chiD==0, 1, kronecker(chiD,p));

/* run the tower for one row and one prime.
   alist = list of tower bases a (p nmid a); smax = top level.
   Returns a vector of records. */
{
towerrun(name, p, alist, smax, K) =
  my(R, typ, aa, bb, cc, dd, w, chiD, N, K0, M, vN, targ, Ah, Ap, Bh, Bp, nx,
     FU, VF, res, s, n, v, ch, recA, recB, recFU, recVF);
  R = mapget(ROWS, name);
  typ=R[1]; aa=R[2]; bb=R[3]; cc=R[4]; dd=R[5]; w=R[6]; chiD=R[7];
  N = vecmax(alist) * p^smax;
  vN = (N - sumdigits(N,p))/(p-1);          /* v_p(N!) */
  K0 = w*vN + K + 40;
  M  = p^K0;
  /* target index set */
  targ = List();
  for(i=1,#alist, for(t=0,smax, listput(targ, alist[i]*p^t)));
  targ = vecsort(Vec(targ),,8);
  recA = Map(); recB = Map(); recFU = Map(); recVF = Map();
  Ap = Mod(1,M); Ah = Mod(bb,M);
  Bp = Mod(0,M); Bh = Mod(1,M);
  FU = Mod(1,M); VF = 0;                    /* unit part of n! and v_p(n!), n=1 */
  /* record n=1 */
  if(vecsearch(targ,1), mapput(recA,1,lift(Ah)); mapput(recB,1,lift(Bh));
     mapput(recFU,1,lift(FU)); mapput(recVF,1,VF));
  for(n=1,N-1,
    if(typ==2,
      nx = (aa*n^2+aa*n+bb)*Ah - cc*n^4*Ap; Ap=Ah; Ah=nx;
      nx = (aa*n^2+aa*n+bb)*Bh - cc*n^4*Bp; Bp=Bh; Bh=nx
    ,
      nx = (2*n+1)*(aa*n^2+aa*n+bb)*Ah - (cc*n^2+dd)*n^4*Ap; Ap=Ah; Ah=nx;
      nx = (2*n+1)*(aa*n^2+aa*n+bb)*Bh - (cc*n^2+dd)*n^4*Bp; Bp=Bh; Bh=nx
    );
    v = valuation(n+1,p);
    FU = FU * ((n+1)/p^v); VF = VF + v;
    if(vecsearch(targ,n+1),
       mapput(recA,n+1,lift(Ah)); mapput(recB,n+1,lift(Bh));
       mapput(recFU,n+1,lift(FU)); mapput(recVF,n+1,VF));
  );
  [name, p, w, chip(chiD,p), K0, alist, smax, recA, recB, recFU, recVF];
}

/* --- reporting helpers ------------------------------------------------ */

/* p-adic ratio x/y from lifted residues mod p^K0, returned as t_PADIC */
{
pratio(x,y,p,K0) =
  my(vx,vy,ux,uy,pr);
  vx = if(x==0, K0, valuation(x,p));
  vy = if(y==0, K0, valuation(y,p));
  pr = K0 - max(vx,vy) - 2;
  if(pr<=0, return(0));
  ux = (x/p^vx) + O(p^pr);
  uy = (y/p^vy) + O(p^pr);
  p^(vx-vy) * ux/uy;
}


/* compact machine-readable report.
   Emits, per tower base a:
     LAM  row p a w chi vLam ndig  <integer u with Lambda = p^vLam*(u + O(p^ndig))>
     AGR  row p a  <agreement digits of L_s vs L_{s-1}, s=1..smax>
     AGRM row p a  <same but with the character factor omitted>
     UDF  row p a  <v_p(A(ap^{s+1})/A(ap^s) - 1), s=0..smax-1>
     FLR  row p a  v_p(Lambda_a - b_a/a_a)                                  */
{
report(T) =
  my(name,p,w,ch,K0,alist,smax,recA,recB,recFU,recVF,a,s,n,np,L,Lp,U,V,Vp,
     agr,agrm,udf,bdf,vL,nd,fa,lam);
  name=T[1]; p=T[2]; w=T[3]; ch=T[4]; K0=T[5]; alist=T[6]; smax=T[7];
  recA=T[8]; recB=T[9]; recFU=T[10]; recVF=T[11];
  for(i=1,#alist, a=alist[i];
    agr=[]; agrm=[]; udf=[]; Lp=0; Vp=0; lam=0;
    for(s=0,smax, n=a*p^s;
      V = pratio(mapget(recB,n), mapget(recA,n), p, K0) * p^(w*s);
      L = ch^s * V;
      if(s>0,
         agr =concat(agr, if(L==Lp,K0,valuation(L-Lp,p)));
         agrm=concat(agrm,if(V==Vp,K0,valuation(V-Vp,p))));
      Lp=L; Vp=V;
    );
    lam = Lp;
    vL = valuation(lam,p);
    nd = w*(smax+1);            /* certified by Thm: Lambda = L_smax mod p^{w(smax+1)} */
    if(#agr>0 && agr[#agr] < nd, nd = agr[#agr]);
    print("LAM ",name," ",p," ",a," ",w," ",ch," ",vL," ",nd," ",
          lift(Mod(truncate(lam/p^vL), p^nd)));
    print("AGR ",name," ",p," ",a," ",agr);
    print("AGRM ",name," ",p," ",a," ",agrm);
    for(s=0,smax-1, n=a*p^s; np=a*p^(s+1);
      U = pratio(mapget(recA,np), mapget(recA,n), p, K0)
          * (lift((Mod(mapget(recFU,n),p^K0)/Mod(mapget(recFU,np),p^K0))^w) + O(p^(K0-2)))
          * p^(w*(mapget(recVF,n)-mapget(recVF,np)));
      udf = concat(udf, if(U==1, K0, valuation(U-1,p)));
    );
    print("UDF ",name," ",p," ",a," ",udf);
    /* second Frobenius eigenvalue: chi(p)^{-1} p^w B(ap^{s+1})/B(ap^s) -> 1 ? */
    bdf=[];
    for(s=0,smax-1, n=a*p^s; np=a*p^(s+1);
      U = ch * p^w * pratio(mapget(recB,np), mapget(recB,n), p, K0)
          * (lift((Mod(mapget(recFU,n),p^K0)/Mod(mapget(recFU,np),p^K0))^w) + O(p^(K0-2)))
          * p^(w*(mapget(recVF,n)-mapget(recVF,np)));
      bdf = concat(bdf, if(U==1, K0, valuation(U-1,p)));
    );
    print("BDF ",name," ",p," ",a," ",bdf);
    fa = pratio(mapget(recB,a), mapget(recA,a), p, K0);
    print("FLR ",name," ",p," ",a," ",if(lam==fa,K0,valuation(lam-fa,p))," fa=",
          lift(Mod(truncate(fa/p^valuation(fa,p)),p^8))," vfa=",valuation(fa,p));
  );
}

/* ---------------------------------------------------------------------
   eig(T): assumption-free measurement of the two tower Frobenius
   eigenvalues.  For each tower base a it prints
     EIGA row p a  <v_p(rho^A_s)>  <v_p(rho^A_s/rho^A_{s-1}-1)>  limit
     EIGB row p a  <v_p(rho^B_s)>  <v_p(rho^B_s/rho^B_{s-1}-1)>  limit
   where rho^X_s = X(ap^{s+1})/X(ap^s).  A convergent rho with v_p(rho)=-e
   says the eigenvalue is p^{-e} times the printed unit.
   ------------------------------------------------------------------ */
{
eig(T) =
  my(name,p,w,ch,K0,alist,smax,recA,recB,recFU,recVF,a,s,n,np,rec,rho,rhop,
     vv,cc2,lim,fac);
  name=T[1]; p=T[2]; w=T[3]; K0=T[5]; alist=T[6]; smax=T[7];
  recA=T[8]; recB=T[9]; recFU=T[10]; recVF=T[11];
  for(i=1,#alist, a=alist[i];
    for(which=1,2, rec = if(which==1,recA,recB);
      vv=[]; cc2=[]; rhop=0; lim=0;
      for(s=0,smax-1, n=a*p^s; np=a*p^(s+1);
        fac = (lift((Mod(mapget(recFU,n),p^K0)/Mod(mapget(recFU,np),p^K0))^w)
               + O(p^(K0-2))) * p^(w*(mapget(recVF,n)-mapget(recVF,np)));
        rho = pratio(mapget(rec,np), mapget(rec,n), p, K0) * fac;
        vv = concat(vv, valuation(rho,p));
        if(s>0, cc2 = concat(cc2, if(rho==rhop,K0,valuation(rho/rhop-1,p))));
        rhop = rho; lim = rho;
      );
      print(if(which==1,"EIGA ","EIGB "),name," ",p," ",a," v=",vv," conv=",cc2,
            if(lim==0, "  lim=0 (degenerate)",
               Str("  lim=",lim/p^valuation(lim,p)+O(p^10)," * p^",valuation(lim,p))));
    );
  );
}

/* ---------------------------------------------------------------------
   cert(T): decides, for each tower base a, which of
      H1: rho^B_s -> 1              (unipotent: slope prime, inner source)
      H2: rho^B_s -> chi(p) p^{-w}  (ordinary: good prime)
      H3: rho^B_s -> 0              (degenerate: outer / cuspidal source)
   holds, by measuring v_p(rho/target - 1) (resp. v_p(rho)) across levels.
   Prints CERT row p a c slope H <max valuation reached> <increment>.
   ------------------------------------------------------------------ */
{
cert(T) =
  my(name,p,w,ch,K0,alist,smax,recA,recB,recFU,recVF,a,s,n,np,rho,fac,
     e1,e2,e3,eA,vv,verdict,cval);
  name=T[1]; p=T[2]; w=T[3]; ch=T[4]; K0=T[5]; alist=T[6]; smax=T[7];
  recA=T[8]; recB=T[9]; recFU=T[10]; recVF=T[11];
  cval = mapget(ROWS,name)[4];
  for(i=1,#alist, a=alist[i];
    e1=[]; e2=[]; e3=[]; eA=[];
    for(s=0,smax-1, n=a*p^s; np=a*p^(s+1);
      fac = (lift((Mod(mapget(recFU,n),p^K0)/Mod(mapget(recFU,np),p^K0))^w)
             + O(p^(K0-2))) * p^(w*(mapget(recVF,n)-mapget(recVF,np)));
      rho = pratio(mapget(recB,np), mapget(recB,n), p, K0) * fac;
      e1 = concat(e1, if(rho==1,K0,valuation(rho-1,p)));
      e2 = concat(e2, if(ch==0, -1, if(rho==ch*p^(-w),K0,valuation(ch*p^w*rho-1,p))));
      e3 = concat(e3, if(rho==0,K0,valuation(rho,p)));
      rho = pratio(mapget(recA,np), mapget(recA,n), p, K0) * fac;
      eA = concat(eA, if(rho==1,K0,valuation(rho-1,p)));
    );
    verdict = "??";
    if(#e1>=3,
      if(e1[#e1]>=e1[1]+2*(#e1-1)-1, verdict="H1:(1,1)");
      if(ch!=0 && e2[#e2]>=e2[1]+2*(#e2-1)-1, verdict="H2:(1,chi*p^-w)");
      if(e3[#e3]>=e3[1]+2*(#e3-1)-1 && e3[#e3]>3, verdict="H3:(1,0)"));
    print("CERT ",name," ",p," ",a," c=",cval," ",
          if(cval%p==0,"SLOPE","good ")," ",verdict,
          "  eA=",eA," e1=",e1," e2=",e2," e3=",e3);
  );
}
