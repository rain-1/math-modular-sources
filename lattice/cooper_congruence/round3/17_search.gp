\\ 17_search.gp -- is the s7 weight-5/2 input a weakly holomorphic Jacobi form of weight 3,
\\ index 7, with poles only at the cusp?   Model: F = psi / Delta^n,  psi in J^weak_{3+12n,7}
\\ = phi_{-1,2} * (sum)_{b=0..5} M_{4+2b+12n} phi_{0,1}^{5-b} phi_{-2,1}^b.
\\ Dictionary:  Jacobi coefficient c(N,r), D = 28N - r^2 ;  D = 3d, r = beta (mod 14);
\\ c(d) = round 2's twisted CM trace at class bt(d), antisymmetric in beta.
default(parisize, 6000000000);
read("15_jaclib.gp");
P21 = mkphi21(); P01 = mkphi01(); P12 = mkphi12();

\\ basis of M_w(SL_2(Z)):  Delta^j E4^a E6^b,  12j+4a+6b = w,  b in {0,1}
{ Mbasis(w) = my(L=List());
  if(w<0 || w%2, return([]));
  if(w==0, return([one]));
  if(w==2, return([]));
  for(j=0, w\12, my(rem = w-12*j);
    for(b=0,1, my(r2 = rem-6*b);
      if(r2>=0 && r2%4==0, listput(L, DELS^j * E4S^(r2/4) * E6S^b))));
  Vec(L); }

{ Jbasis(n) = my(L=List(), base, MB);
  for(b=0,5,
    base = jmul(P12, jmul(jpow(P01,5-b), jpow(P21,b)));
    MB = Mbasis(4+2*b+12*n);
    for(i=1,#MB, listput(L, jscal(MB[i], base))));
  Vec(L); }

DMAX = 2500;
DAT = read("../round2/73_cd_s7_2500.txt");
CD = vector(DMAX); BT = vector(DMAX);
{ for(i=1,#DAT, my(e=DAT[i], d=e[1], v=e[2], bs=List(), bt);
    if(d>DMAX, next);
    for(b=0,13, if((b^2+3*d)%28==0, listput(bs,b)));
    bs = Vec(bs);
    bt = if(issquare(d), (5*sqrtint(d))%14, bs[1]);
    if(!setsearch(Set(bs), bt), bt = bs[1]);
    BT[d] = bt; CD[d] = if(type(v)=="t_STR", 0, v)); }
{ tgt(D, r) = my(d, bt, rr = r % 14);
  if(D % 3, return(0));
  d = D/3;
  if(d<1 || d>DMAX, return("UNK"));
  if(CD[d]==0, return(0));
  bt = BT[d];
  if(bt==rr, return(CD[d]));
  if((14-bt)%14==rr, return(-CD[d]));
  0; }

POL = 0;
{ polval(D,rr,j) = my(P=POL[j]);
   if(P[1]!=D, return(0));
   if(P[2]==rr%14, return(1));
   if((14-P[2])%14==rr%14, return(-1));
   0; }

{ search(n, NMAX) =
  my(BS = Jbasis(n), nb, np, A, Y, row, V, dc, DN);
  nb = #BS;
  DN = DELS^n;
  dc = vector(NMAX+1, i, polcoeff(DN, 2*(i-1), 'y));
  POL = List();
  forstep(D=-1, -49-28*n, -1,
    for(rr=0,7, if((D+rr^2)%28==0, listput(POL,[D,rr]))));
  POL = Vec(POL); np = #POL;
  A = matrix(14*(NMAX+1), nb+np); Y = vector(14*(NMAX+1));
  row = 0;
  for(N=0,NMAX, for(rr=0,13,
    row++;
    for(i=1,nb, A[row,i] = jcoef(BS[i], N, rr));
    my(s=0, bad=0);
    for(j=0,N,
      my(M=N-j, D=28*M-rr^2, t);
      if(dc[j+1]==0, next);
      if(D<0,
        for(k=1,np, my(pv=polval(D,rr,k)); if(pv, A[row,nb+k] = A[row,nb+k] - dc[j+1]*pv))
      ,
        t = tgt(D,rr);
        if(type(t)=="t_STR", bad=1; break);
        s = s + dc[j+1]*t));
    if(bad, row=row-1; next);
    Y[row] = s));
  A = A[1..row,]; Y = vector(row, i, Y[i]);
  print("n=",n,"  dim=",nb,"  polar unknowns=",np,"  equations=",row,"  rank=",matrank(A));
  V = matinverseimage(A, Y~);
  if(#V==0, print("   NO SOLUTION"); return(0));
  print("   SOLUTION FOUND.  kernel dim = ", matsize(matker(A))[2]);
  print("   polar coefficients: ", vector(np,k,[POL[k], V[nb+k]]));
  V; }

print("dim J^weak_{3,7} = ", #Jbasis(0), ",  dim J^weak_{15,7} = ", #Jbasis(1));
search(0, 40);
search(1, 40);
quit;
