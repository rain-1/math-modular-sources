\\ 23_full.gp -- the identification test with the FULL c(d) data (all admissible d, not only
\\ squares), using only the slots D = 3d and leaving the coefficients at D not divisible by 3
\\ completely free.  Rows: c_F(3d, bt(d)) = kappa * c_data(d).  Homogeneous in (v, kappa).
default(parisize, 12000000000);
read("15_jaclib.gp");
P21 = mkphi21(); P01 = mkphi01(); P12 = mkphi12();
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
{ dinv(n) = my(D = DELS^n, u = vector(NQ+1), s);
  for(i=0,NQ, u[i+1] = polcoeff(D, 2*(i+n), 'y));
  s = Ser(u, 'y, NQ+1);
  Vec(1/s + O('y^(NQ+1))); }
{ Fc(BS, IU, n, i, D, r) = my(N=(D+r^2)/28);
  if(N+n < 0, return(0));
  sum(j=0, N+n, IU[j+1]*jcoef(BS[i], N+n-j, r)); }
DMAX = 2500;
DAT = read("../round2/73_cd_s7_2500.txt");
CD = vector(DMAX); BT = vector(DMAX); OK = vector(DMAX);
{ for(i=1,#DAT, my(e=DAT[i], d=e[1], v=e[2], bs=List(), bt);
    if(d>DMAX, next);
    for(b=0,13, if((b^2+3*d)%28==0, listput(bs,b)));
    bs = Vec(bs);
    bt = if(issquare(d), (5*sqrtint(d))%14, bs[1]);
    if(!setsearch(Set(bs), bt), bt = bs[1]);
    if(bt==0 || bt==7, next);
    BT[d]=bt; OK[d]=1; CD[d] = if(type(v)=="t_STR", 0, v)); }
{ run(n) =
  my(BS = Jbasis(n), nb = #Jbasis(n), IU = dinv(n), rows=List(), MT, K);
  for(d=1,DMAX,
    if(!OK[d], next);
    my(D=3*d, r=BT[d], N=(D+r^2)/28, rw);
    if(r>7, r = 14-r);   \\ use beta in 1..6, flipping the sign of the datum
    N = (D+r^2)/28;
    if(N+n > NQ, break);
    rw = vector(nb+1);
    for(i=1,nb, rw[i] = Fc(BS,IU,n,i,D,r));
    rw[nb+1] = -(if(BT[d]>7, -1, 1))*CD[d];
    listput(rows, rw));
  MT = matrix(#rows, nb+1, i, j, rows[i][j]);
  K = matker(MT);
  print("n=",n,"  dim=",nb,"  eqs=",#rows,"  rank=",matrank(MT),"  kernel=",matsize(K)[2]);
  \\ longest matchable prefix
  my(best=0);
  for(L=1,#rows, my(MM=matrix(L,nb+1,i,j,rows[i][j]), KK=matker(MM));
     if(matsize(KK)[2]>0 && vecmax(vector(matsize(KK)[2],c,abs(KK[nb+1,c])))>0, best=L, break));
  print("   longest matchable prefix (in d) = ", best);
}
run(0); run(1); run(2); run(3);
quit;
