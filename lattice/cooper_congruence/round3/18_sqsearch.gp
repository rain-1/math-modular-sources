\\ 18_sqsearch.gp -- decisive test using ONLY the square-index data beta_{s7}(m).
\\ Model:  F = psi/Delta^n,  psi in J^weak_{3+12n,7};  the Jacobi coefficient of F at
\\ D = 3 m^2, r = 5m (mod 14) must be proportional to beta_{s7}(m).
\\ Homogeneous system in (v, lambda); a nonzero kernel with lambda != 0 identifies F.
default(parisize, 6000000000);
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
BETA = read("../round2/beta_s7.txt");
print("beta(1..10) = ", vector(10,i,BETA[i]));

\\ q-coefficient vector of Delta^n divided by q^n (a unit series), and its inverse
{ dinv(n) = my(D = DELS^n, u = vector(NQ+1), s);
  for(i=0,NQ, u[i+1] = polcoeff(D, 2*(i+n), 'y));
  s = Ser(u, 'y, NQ+1);
  Vec(1/s + O('y^(NQ+1))); }

{ run(n, MMAX, t) =
  my(BS = Jbasis(n), nb, IU = dinv(n), rows = List(), MT, K, ok);
  nb = #BS;
  \\ coefficient of F_i = psi_i/Delta^n at (D, r): shift q-index by -n and multiply by IU
  \\ Fc(i,N,r) = sum_{j>=0} IU[j+1] * jcoef(psi_i, N+n-j, r)
  for(m=1, MMAX,
    if(m%7==0, next);
    my(r = (5*m)%14, D = 3*m^2, N, rw);
    N = (D + r^2)/28;
    if(N != floor(N), print("  bad N at m=",m); next);
    if(N + n > NQ, break);
    rw = vector(nb+1);
    for(i=1,nb, rw[i] = sum(j=0, N+n, IU[j+1]*jcoef(BS[i], N+n-j, r)));
    rw[nb+1] = -m^t*BETA[m];
    listput(rows, rw));
  MT = matrix(#rows, nb+1, i, j, rows[i][j]);
  print("n=",n," t=",t,"  dim=",nb,"  equations=",#rows,"  rank=",matrank(MT));
  K = matker(MT);
  print("   kernel dim = ", matsize(K)[2]);
  if(matsize(K)[2]>0, for(c=1,matsize(K)[2], print("   kernel vector ",c,": lambda = ", K[nb+1,c], "  v = ", vector(nb,i,K[i,c]))));
  K; }

foreach([-3,-2,-1,0,1,2,3], t, foreach([0,1,2], n, run(n,40,t)));
quit;
