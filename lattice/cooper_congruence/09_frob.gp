/* 09_frob.gp -- (a) the Frobenius lift  X_sigma(x) := x(q(x)^p)  is p-integral,
                 so sigma: x -> X_sigma lifts x -> x^p and sigma_* lifts the Cartier
                 operator;  (S) <=> sigma_* eta = psi(p) eta mod p^2.
                 (b) s18 at p=3: the exact 3-adic structure of Phi|U_3.            */
default(parisize, 4000000000);
read("lib.gp");
MM = 120;
print("=== (a) X_sigma(x) = x(q(x)^p):  p-integrality of the Frobenius lift ===");
{ for(k=1,3,
   my(S=Setup(k,MM+4), xq=S[3], qx, row=[]);
   qx = serreverse(xq);                       \\ q as a series in x (variable q)
   forprime(p=2,7,
     my(sub1, Xs, den=1, bad=0);
     sub1 = qx^p;                              \\ q^p as a series in x
     Xs = subst(xq, 'q, sub1);                 \\ x(q^p) as a series in x
     for(j=0, MM\p-2, my(cc=polcoeff(Xs,j)); if(cc!=0 && valuation(denominator(cc),p)>0, bad++));
     row=concat(row,[[p, bad]]));
   print("  ",NAM[k]," [p, # coefficients of X_sigma with p in the denominator] : ",row)); }
print();
print("   X_sigma for s7, p=3, first terms (should start x^3 + ...) :");
{ my(S=Setup(1,MM+4), xq=S[3], qx=serreverse(S[3]), Xs=subst(xq,'q,qx^3));
  print("   ",Ser(vector(9,j,polcoeff(Xs,j-1)),'x)); }
print();
print("=== (b) s18 at p=3:  v_3 of c(3^n m) ===");
{ my(S=Setup(3,1000), CV=vector(994,m,polcoeff(S[4],m)));
  for(n=1,5, my(pn=3^n, mn=10^9, cnt=0);
    if(pn>994, break);
    for(m=1,994\pn, cnt++; mn=min(mn, valuation(CV[pn*m],3)));
    print("   n=",n,"  min v_3(c(3^n m)) = ",mn,"  (3n = ",3*n,")   #tests=",cnt));
  print("   first terms of Phi|U_3 / 27 : ", vector(12,m,CV[3*m]/27));
  print("   first terms of Phi          : ", vector(12,m,CV[m]));
  print("   ratio test  c(3m)/27 vs c(m):  not proportional (see above)");
  print("   v_3(c'(3^n m)) minima:");
  for(n=1,4, my(pn=3^n, mn=10^9);
    if(pn>994,break);
    for(m=1,994\pn, mn=min(mn, valuation(CV[pn*m]/(pn*m),3)));
    print("     n=",n,"  min v_3(c'(3^n m)) = ",mn,"   (2n = ",2*n,")")); }
print();
print("=== (c) the AL/character coincidence:  eps_p * psi(p) over the bad primes ===");
print("   s7  : p=7 : eps=-1, psi=+1 -> eps*psi = -1  (EXACT cell)");
print("   s10 : p=2 : eps=+1, psi=+1 -> eps*psi = +1");
print("   s10 : p=5 : eps=-1, psi=+1 -> eps*psi = -1  (EXACT cell)");
print("   s18 : p=2 : eps=-1, psi=-1 -> eps*psi = +1");
quit;
