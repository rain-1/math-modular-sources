\\ 21_final.gp -- (a) index-7 consistency of the constructed basis of J^weak_{3,7};
\\ (b) the square-index search repeated with a larger q-precision and odd Jacobi weights
\\     kappa = 1,3,5 and pole orders n = 0..4.
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
{ Jbasis(kap, n) = my(L=List(), base, MB);
  for(b=0,5,
    base = jmul(P12, jmul(jpow(P01,5-b), jpow(P21,b)));
    MB = Mbasis(kap+1+2*b+12*n);
    for(i=1,#MB, listput(L, jscal(MB[i], base))));
  Vec(L); }
\\ index consistency: c(N,r) depends only on (28N-r^2, r mod 14)
{ ichk(A) = my(bad=0);
   for(N=0,8, for(r=-20,20,
     my(D=28*N-r^2, c=jcoef(A,N,r));
     for(N2=0,8, for(r2=-20,20,
       if(28*N2-r2^2==D && (r-r2)%14==0 && jcoef(A,N2,r2)!=c, bad++)))));
   bad; }
BS0 = Jbasis(3,0);
print("index-7 consistency violations in the J^weak_{3,7} basis: ", vector(#BS0, i, ichk(BS0[i])));
BETA = read("../round2/beta_s7.txt");
{ dinv(n) = my(D = DELS^n, u = vector(NQ+1), s);
  for(i=0,NQ, u[i+1] = polcoeff(D, 2*(i+n), 'y));
  s = Ser(u, 'y, NQ+1);
  Vec(1/s + O('y^(NQ+1))); }
{ run(kap, n, t) =
  my(BS = Jbasis(kap,n), nb, IU = dinv(n), rows = List(), MT, K, nm=0);
  nb = #BS;
  for(m=1, 200,
    if(m%7==0, next);
    my(r = (5*m)%14, D = 3*m^2, N, rw);
    N = (D + r^2)/28;
    if(N + n > NQ, break);
    rw = vector(nb+1);
    for(i=1,nb, rw[i] = sum(j=0, N+n, IU[j+1]*jcoef(BS[i], N+n-j, r)));
    rw[nb+1] = -m^t*BETA[m];
    listput(rows, rw); nm++);
  MT = matrix(#rows, nb+1, i, j, rows[i][j]);
  K = matker(MT);
  print("kappa=",kap," n=",n," t=",t,"  dim=",nb,"  eqs=",#rows,"  rank=",matrank(MT),"  kernel=",matsize(K)[2]);
  if(matsize(K)[2]>0 && #rows > nb, for(c=1,matsize(K)[2], print("   lambda=",K[nb+1,c]))); }
foreach([3,1,5], kap, foreach([0,1,2,3,4], n, foreach([0,1], t, run(kap,n,t))));
quit;
