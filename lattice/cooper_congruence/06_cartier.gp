/* 06_cartier.gp -- ROUTE A.  The residue/Cartier formulation.
   c'(m) = Res_x( eta * q(x)^{-m} ),   eta = l(x) dx / (x sqrt(P(x)) F(x)),
   P(x) = 1 - 2Bx + (B^2-4C)x^2,  F(x) = sum A_n x^n,  l(x) = sum A_n x^{n+1}/(n+1).
   Cartier:  C(sum a_j x^j dx) = sum a_{pj+p-1} y^j dy   satisfies
             Res_x(w Q(x^p)) = Res_y(C(w) Q(y)).
   TEST:  C(eta) = psi(p) eta  mod p  (and mod p^2 ?).                            */
default(parisize, 4000000000);
read("lib.gp");
NMAX = 700;
{ ETA = vector(3); FX = vector(3); LX = vector(3); AA = vector(3);
  for(k=1,3, my(R=ROWS[k], B=R[3], C=R[4], A, F, l, P, sq, e);
    A = genrow(k,NMAX)[1];       \\ A[n+1] = A_n
    AA[k] = A;
    F = Ser(A,'x);
    l = 'x*Ser(vector(NMAX+1,n,A[n]/n),'x);
    P = 1 - 2*B*'x + (B^2-4*C)*'x^2 + O('x^(NMAX+1));
    sq = sqrt(P);
    e = l/('x*sq*F);            \\ eta = e(x) dx
    ETA[k]=e; FX[k]=F; LX[k]=l); }
print("=== sanity: sqrt(P) integral?  first coefficients ===");
{ for(k=1,3, my(R=ROWS[k],B=R[3],C=R[4], sq=sqrt(1-2*B*'x+(B^2-4*C)*'x^2+O('x^8)));
   print("  ",NAM[k],"  P = 1 - ",2*B,"x + ",B^2-4*C,"x^2 ;  sqrt(P) = ",sq)); }
print();
print("=== eta = sum a_j x^j : first coefficients ===");
{ for(k=1,3, print("  ",NAM[k],": ",vector(10,j,polcoeff(ETA[k],j-1)))); }
print();
print("=== TEST  C(eta) = psi(p) eta  mod p^e :  a_{pj+p-1} - psi(p) a_j ===");
{ for(k=1,3, my(row=[]);
   forprime(p=2,43, my(ps=psival(k,p), g=0, cnt=0);
     for(j=0, (NMAX-p+1)\p, my(i=p*j+p-1); if(i>NMAX,break); cnt++;
       g = gcd(g, polcoeff(ETA[k],i) - ps*polcoeff(ETA[k],j)));
     row=concat(row,[[p, if(g==0,"inf",valuation(g,p)), cnt]]));
   print("  ",NAM[k],": ",row)); }
print();
print("=== control: same test with psi replaced by the other characters ===");
{ for(k=1,3, my(row=[]);
   forprime(p=2,13, my(g1=0,gm=0);
     for(j=0,(NMAX-p+1)\p, my(i=p*j+p-1); if(i>NMAX,break);
       g1=gcd(g1, polcoeff(ETA[k],i) - polcoeff(ETA[k],j));
       gm=gcd(gm, polcoeff(ETA[k],i) + polcoeff(ETA[k],j)));
     row=concat(row,[[p, if(g1==0,"inf",valuation(g1,p)), if(gm==0,"inf",valuation(gm,p))]]));
   print("  ",NAM[k]," [p, v_p for psi=+1, v_p for psi=-1] : ",row)); }
print();
print("=== also test the Cartier eigen-property of the naive form l(x) dx and F(x) dx ===");
{ for(k=1,3, my(r1=[],r2=[]);
   forprime(p=2,13, my(g1=0,g2=0);
     for(j=0,(NMAX-p+1)\p, my(i=p*j+p-1); if(i>NMAX,break);
       g1=gcd(g1, polcoeff(LX[k],i)-psival(k,p)*polcoeff(LX[k],j));
       g2=gcd(g2, polcoeff(FX[k],i)-psival(k,p)*polcoeff(FX[k],j)));
     r1=concat(r1,[[p,if(g1==0,"inf",valuation(g1,p))]]);
     r2=concat(r2,[[p,if(g2==0,"inf",valuation(g2,p))]]));
   print("  ",NAM[k]," l dx: ",r1);
   print("  ",NAM[k]," F dx: ",r2)); }
quit;
