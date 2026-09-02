/* 43_jslot.gp -- ROUTE (B): the j=0 slot and its refinements.
   a_{p-1} = psi(p) mod p  <=>  c'(p) = psi(p) mod p  <=>  c(p) = psi(p) p mod p^2.
   Tests: (i) a_{p-1} mod p^2  (does the congruence lift?);
          (ii) the j=1,2 slots a_{2p-1}, a_{3p-1} mod p and mod p^2;
          (iii) the defect  d(p) := (a_{p-1} - psi(p))/p  mod p  -- is it recognisable?
          (iv) the analogous statement for the Apery-like numbers themselves:
               A_{p-1} mod p^2, mod p^3.                                             */
default(parisize,8000000000);
read("40_core.gp");
NX = 640;
{ AJv = vector(3); AVv = vector(3);
  for(k=1,3, AJv[k]=AJ(k,NX); AVv[k]=AVEC(k,NX)); }
print("=== (i) a_{p-1} mod p and mod p^2 ===");
{ for(k=1,3, my(a=AJv[k], row=[], bad2=[]);
   forprime(p=2,200,
     my(ps=psival(k,p));
     if(Mod(a[p],p) != Mod(ps,p), row=concat(row,[[p,"FAIL mod p"]]));
     if(Mod(a[p],p^2) == Mod(ps,p^2), bad2=concat(bad2,[p])));
   print("  ",NAM[k],": mod p failures for p<=200: ",if(#row==0,"NONE",row));
   print("     primes where a_{p-1} = psi(p) mod p^2 also holds: ",if(#bad2==0,"NONE",bad2))); }
print();
print("=== (iii) the defect d(p) = (a_{p-1}-psi(p))/p mod p,  p <= 60 ===");
{ for(k=1,3, my(a=AJv[k], row=[]);
   forprime(p=2,60, my(ps=psival(k,p)); row=concat(row,[[p, lift(Mod((a[p]-ps)/p,p))]]));
   print("  ",NAM[k],": ",row)); }
print();
print("=== (ii) j=1 and j=2 slots:  a_{2p-1} - psi(p)a_1  and  a_{3p-1} - psi(p)a_2, v_p ===");
{ for(k=1,3, my(a=AJv[k], row=[]);
   forprime(p=3,97, if(3*p-1>NX,break);
     my(ps=psival(k,p), d1=a[2*p]-ps*a[2], d2=a[3*p]-ps*a[3]);
     row=concat(row,[[p, if(d1==0,"inf",valuation(d1,p)), if(d2==0,"inf",valuation(d2,p))]]));
   print("  ",NAM[k]," [p, v_p(j=1 defect), v_p(j=2 defect)]: ",row)); }
print();
print("=== (iv) A_{p-1} mod p^2 : the pattern  A_{p-1} = kappa * psi(p) * p  (mod p^2) ===");
{ for(k=1,3, my(A=AVv[k], row=[], kap=[]);
   forprime(p=3,200, my(ps=psival(k,p), t);
     if(A[p]%p!=0, row=concat(row,[[p,"p does not divide A_{p-1}"]]); next);
     t = lift(Mod((A[p]/p)*ps, p));   \\ psi(p)^{-1}=psi(p)
     kap=concat(kap,[if(t>p\2,t-p,t)]));
   print("  ",NAM[k],"  kappa(p) := psi(p)*A_{p-1}/p mod p (centred), p=3..200:");
   print("     ",kap)); }
print();
print("=== (iv') A_{p-1} mod p^3, i.e. is  A_{p-1}/p + kappa  divisible by p ? ===");
{ my(K=[-2,-3,-3]);
  for(k=1,3, my(A=AVv[k], row=[]);
   forprime(p=3,60, my(ps=psival(k,p));
     row=concat(row,[[p, lift(Mod(A[p]/p - K[k]*ps, p^2))]]));
   print("  ",NAM[k]," kappa=",K[k],": [p, (A_{p-1}/p - kappa*psi(p)) mod p^2] = ",row)); }
quit;
