/* 42_upx.gp -- ROUTE (A) ON THE x-LINE, done right.
   CLAIM 1 (Kronecker/Igusa): for a modular FUNCTION f on Gamma_0(N) with f in Z[[q]]
       and p not dividing N,   U_p f == 0  (mod p).
       [Reason: mod p the modular correspondence is Phi_p(X,Y) = (X^p-Y)(X-Y^p), so
        T_p f == f^p + p*(...) == V_p f (mod p), and U_p f = T_p f - V_p f.]
   CLAIM 2: hence, mod p,  Xi|U_p = sum_a x^a sum_b L_{pa+b} U_p(x^b) == sum_a L_{pa} x^a,
       so (S) mod p  <=>  L_{pa} == psi(p) L_a (mod p) for all a >= 1,  L_n = A_{n-1}/n.
   TEST both, on q-expansions, exactly.                                                */
default(parisize,6000000000);
read("40_core.gp");
M = 340;                    \\ number of q-coefficients
print("=== CLAIM 1:  U_p(x^b) == 0 mod p  for 1<=b<=p-1  (p not dividing N) ===");
print("    row : [p, max b tested, #nonzero U_p(x^b) mod p]");
{ for(k=1,3, my(S=Setup(k,M), x=S[3], N=LEV[k], row=[]);
   forprime(p=2,29, if(N%p==0, row=concat(row,[[p,"p|N"]]); next);
     my(bad=0, bmax=min(p-1,12), xb=1+O('q^M));
     for(b=1,bmax, xb = xb*x;
       my(nz=0);
       for(m=1,(M-4)\p, if(Mod(polcoeff(xb,p*m),p)!=0, nz=1; break));
       if(nz, bad++));
     row=concat(row,[[p,bmax,bad]]));
   print("  ",NAM[k],": ",row)); }
print();
print("=== control: U_p(1) = 1, and U_p(x) mod p first few coefficients ===");
{ for(k=1,3, my(S=Setup(k,M), x=S[3], row=[]);
   forprime(p=3,13, if(LEV[k]%p==0,next);
     row=concat(row,[[p,vector(6,m,lift(Mod(polcoeff(x,p*m),p)))]]));
   print("  ",NAM[k],": ",row)); }
print();
print("=== CLAIM 2 (the reduction):  L_{pa} == psi(p) L_a  (mod p),  L_n = A_{n-1}/n ===");
NA = 1200;
{ LL = vector(3);
  for(k=1,3, my(A=AVEC(k,NA)); LL[k] = vector(NA,n, A[n]/n)); }
print("    row : [p, #tests, #failures]");
{ for(k=1,3, my(L=LL[k], row=[]);
   forprime(p=2,97, my(ps=psival(k,p), nf=0, nt=0);
     for(a=1,NA\p, nt++; if(Mod(L[p*a],p)!=Mod(ps*L[a],p), nf++));
     row=concat(row,[[p,nt,nf]]));
   print("  ",NAM[k],": ",row)); }
print();
print("=== equivalently  A_{pa-1} == psi(p)*p*A_{a-1} (mod p^2) for p not dividing a ===");
{ for(k=1,3, my(A=AVEC(k,NA), row=[]);
   forprime(p=2,53, my(ps=psival(k,p), nf=0, nt=0);
     for(a=1,NA\p, if(a%p==0,next); nt++;
       if(Mod(A[p*a], p^2) != Mod(ps*p*A[a],p^2), nf++));
     row=concat(row,[[p,nt,nf]]));
   print("  ",NAM[k],": ",row)); }
print();
print("=== the j=0 case:  A_{p-1} mod p^2  vs  psi(p)*p ===");
{ for(k=1,3, my(A=AVEC(k,NA), row=[]);
   forprime(p=2,60, row=concat(row,[[p, lift(Mod(A[p],p^2)), lift(Mod(psival(k,p)*p,p^2))]]));
   print("  ",NAM[k],": ",row)); }
quit;
