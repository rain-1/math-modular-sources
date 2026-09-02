\\ 20_diag.gp -- (a) Milgram: signature of the discriminant form (Z/14, x^2/28) mod 8;
\\ (b) validation of the Jacobi machinery: J_{3,7} = 0 (no holomorphic weight-3 index-7 form),
\\     and the structure of principal parts inside J^weak_{3,7} and Delta^{-1}J^weak_{15,7}.
default(parisize, 6000000000);
default(realprecision, 40);
S = sum(b=0,13, exp(2*Pi*I*b^2/28));
print("Milgram sum = ", S, "   |S| = ", abs(S), "  sqrt(14) = ", sqrt(14));
print("   sig/8 from arg: ", arg(S)/(2*Pi)*8, "   => sig = ", round(arg(S)/(2*Pi)*8), " mod 8");
\\ same for the dual form Q(x) = -x^2/28
S2 = sum(b=0,13, exp(-2*Pi*I*b^2/28));
print("dual: sig = ", round(arg(S2)/(2*Pi)*8), " mod 8");
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
\\ list of (D,r) with D<0, D = -r^2 mod 28, r in 0..7, D >= -49-28n
{ polidx(n) = my(L=List()); forstep(D=-1,-49-28*n,-1, for(r=0,7, if((D+r^2)%28==0, listput(L,[D,r])))); Vec(L); }
{ holotest(n) = my(BS=Jbasis(n), nb=#BS, PI=polidx(n), A, K);
  A = matrix(#PI, nb, i, j, my(D=PI[i][1], r=PI[i][2], N=(D+r^2)/28); jcoef(BS[j], N+n, r));
  K = matker(A);
  print("n=",n,": dim J^weak_{",3+12*n,",7} = ", nb, ", #principal slots = ", #PI,
        ", rank of principal-part map = ", matrank(A), ", holomorphic subspace dim = ", matsize(K)[2]); }
holotest(0);
holotest(1);
holotest(2);
\\ antisymmetry check on a basis element
{ B = Jbasis(0)[1];
  print("antisymmetry c(N,-r) = -c(N,r): ", vector(6,i,[jcoef(B,3,i), jcoef(B,3,-i)])); }
quit;
