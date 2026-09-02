/* 41_hyper.gp -- THE FINITE (HYPERELLIPTIC) FORM OF THE CARTIER CONGRUENCE.
   eta = a(x) dx = b(x) dx / y,  y = sqrt(P),  b := l(x)/(x F(x)) in Z[[x]], b(0)=1.
   Since y^p = y*P^((p-1)/2) and C(h^p w) = h C(w):
        C(b dx/y) = (1/y) * C_x( b * P^((p-1)/2) dx ).
   Hence the Cartier congruence is EQUIVALENT to the FINITE-DATA statement
        [x^{pj+p-1}] ( b(x) * P(x)^((p-1)/2) )  ==  psi(p) * b_j   (mod p),
   in which the transcendental factor 1/sqrt(P) has been replaced by the
   POLYNOMIAL P^((p-1)/2) of degree p-1.
   Also: [x^{p-1}] P^((p-1)/2) = (B^2-4C)^((p-1)/2) = kronecker(B^2-4C,p) -- the
   "Hasse invariant" of the conic; we test how much of psi(p) that accounts for. */
default(parisize,6000000000);
read("40_core.gp");
NMAX = 1200;
{ BB = vector(3); AJv = vector(3);
  for(k=1,3, my(l=LSER(k,NMAX), F=FSER(k,NMAX));
    BB[k] = l/('x*F);
    AJv[k] = AJ(k,NMAX)); }
print("=== b(x) = l/(xF), first 10 coefficients ===");
{ for(k=1,3, print("  ",NAM[k],": ",vector(10,i,polcoeff(BB[k],i-1)))); }
print();
print("=== consistency: a = b/sqrt(P) ===");
{ for(k=1,3, my(t = BB[k]/sqrt(PSER(k,60)) - ETAS(k,60), ok=1);
   for(i=0,50, if(polcoeff(t,i)!=0, ok=0));
   print("  ",NAM[k],": ",if(ok,"OK","MISMATCH"))); }
print();
print("=== TEST  [x^{pj+p-1}](b*P^((p-1)/2)) = psi(p) b_j  (mod p) ===");
print("    row : [p, #tests, #failures]");
{ for(k=1,3, my(R=ROWS[k],B=R[3],C=R[4], row=[]);
   forprime(p=3,89,
     my(ps=psival(k,p), Ppow, bt, prod, J, nf=0, nt=0, MM);
     J = (NMAX-p+1)\p;
     MM = p*J+p;                       \\ need coefficients up to p*J+p-1
     Ppow = Mod(1,p)*(1 - 2*B*'x + (B^2-4*C)*'x^2)^((p-1)/2);
     bt = Mod(1,p)*Ser(vector(MM,i,polcoeff(BB[k],i-1)),'x);
     prod = bt*Ppow;
     for(j=0,J, my(i=p*j+p-1); if(i>=MM,break); nt++;
       if(polcoeff(prod,i) != ps*Mod(polcoeff(BB[k],j),p), nf++));
     row=concat(row,[[p,nt,nf]]));
   print("  ",NAM[k],": ",row)); }
print();
print("=== the j=0 slot decomposed:  [x^{p-1}](b*P^((p-1)/2)) vs psi(p) ===");
print("    and the pure conic term [x^{p-1}]P^((p-1)/2) = kronecker(B^2-4C,p)");
{ for(k=1,3, my(R=ROWS[k],B=R[3],C=R[4],disc=B^2-4*C, row=[]);
   forprime(p=3,41,
     my(Ppow = Mod(1,p)*(1-2*B*'x+(B^2-4*C)*'x^2)^((p-1)/2),
        s = sum(i=0,p-1, Mod(polcoeff(BB[k],i),p)*polcoeff(Ppow,p-1-i)));
     row=concat(row,[[p,lift(s),psival(k,p),kronecker(disc,p)]]));
   print("  ",NAM[k]," disc=",disc,": [p, slot, psi(p), (disc|p)] = ",row)); }
quit;
