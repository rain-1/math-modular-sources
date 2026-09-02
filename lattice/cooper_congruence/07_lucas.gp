/* 07_lucas.gp -- structure of the Cartier congruence C(eta) = psi(p) eta mod p.
   eta = a(x) dx,  a(x) = l(x)/(x sqrt(P(x)) F(x)).
   (i)  is a_{p-1} = psi(p) mod p ?
   (ii) does (a_j) satisfy a full Lucas congruence a_{pj+r} = a_j a_r mod p ?
   (iii) Dwork factorisation  A(x) = A_{<p}(x) A(x^p)  mod p  for A = sum a_j x^j ?
   (iv) the same tests for the three ingredients l/x, 1/sqrt(P), 1/F.               */
default(parisize, 4000000000);
read("lib.gp");
NMAX = 700;
{ AX = vector(3);
  for(k=1,3, my(R=ROWS[k], B=R[3], C=R[4], A, F, l, P);
    A = genrow(k,NMAX)[1];
    F = Ser(A,'x);
    l = 'x*Ser(vector(NMAX+1,n,A[n]/n),'x);
    P = 1 - 2*B*'x + (B^2-4*C)*'x^2 + O('x^(NMAX+1));
    AX[k] = l/('x*sqrt(P)*F)); }
av(k,j) = polcoeff(AX[k],j);
print("=== (i)  a_{p-1} mod p  vs  psi(p) ===");
{ for(k=1,3, my(row=[]);
   forprime(p=2,61, row=concat(row,[[p, lift(Mod(av(k,p-1),p)), lift(Mod(psival(k,p),p))]]));
   print("  ",NAM[k]," [p, a_{p-1} mod p, psi(p) mod p]:");
   print("   ",row);
   print("   all equal ? ", #select(x->x[2]!=x[3], row)==0)); }
print();
print("=== (ii) full Lucas:  a_{pj+r} = a_j a_r mod p, 0<=r<p ===");
{ for(k=1,3, my(row=[]);
   forprime(p=2,23, my(bad=0,tot=0);
     for(j=0,(NMAX-p)\p, for(r=0,p-1, my(i=p*j+r); if(i>NMAX,next); tot++;
       if((av(k,i)-av(k,j)*av(k,r))%p!=0, bad++)));
     row=concat(row,[[p,bad,tot]]));
   print("  ",NAM[k]," [p, #failures, #tests] : ",row)); }
print();
print("=== (iii) Dwork factorisation A(x) = A_{<p}(x) A(x^p) mod p (equivalent to (ii)) ===");
{ for(k=1,3, my(row=[]);
   forprime(p=2,13, my(MM=200, Ap=Ser(vector(MM,i,av(k,i-1)),'x), Alt, Afr, d);
     Alt = Ser(vector(MM,i,if(i<=p,av(k,i-1),0)),'x);
     Afr = Ser(vector(MM,i,if((i-1)%p==0,av(k,(i-1)/p),0)),'x);
     d = Ap - Alt*Afr;
     row = concat(row,[[p, #select(i->polcoeff(d,i)%p!=0, vector(MM-p,i,i-1))]]));
   print("  ",NAM[k]," [p, #nonzero coefficients of the defect mod p among the first ",200-13,"] : ",row)); }
print();
print("=== (iv) the same Cartier test for the three factors separately ===");
{ for(k=1,3, my(R=ROWS[k],B=R[3],C=R[4], A=genrow(k,NMAX)[1], F=Ser(A,'x),
      l='x*Ser(vector(NMAX+1,n,A[n]/n),'x), P=1-2*B*'x+(B^2-4*C)*'x^2+O('x^(NMAX+1)),
      f1=l/'x, f2=1/sqrt(P), f3=1/F, r1=[],r2=[],r3=[]);
   forprime(p=2,13, my(g1=0,g2=0,g3=0);
     for(j=0,(NMAX-p+1)\p, my(i=p*j+p-1); if(i>NMAX,break);
       g1=gcd(g1,polcoeff(f1,i)-psival(k,p)*polcoeff(f1,j));
       g2=gcd(g2,polcoeff(f2,i)-psival(k,p)*polcoeff(f2,j));
       g3=gcd(g3,polcoeff(f3,i)-psival(k,p)*polcoeff(f3,j)));
     r1=concat(r1,[[p,if(g1==0,"inf",valuation(g1,p))]]);
     r2=concat(r2,[[p,if(g2==0,"inf",valuation(g2,p))]]);
     r3=concat(r3,[[p,if(g3==0,"inf",valuation(g3,p))]]));
   print("  ",NAM[k]," l/x  : ",r1);
   print("  ",NAM[k]," 1/sqrtP: ",r2);
   print("  ",NAM[k]," 1/F  : ",r3)); }
quit;
