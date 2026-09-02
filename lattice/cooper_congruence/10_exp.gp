/* 10_exp.gp -- the q-product / Dieudonne-Dwork picture.
   Xi = D log g  with  g = prod_n (1-q^n)^{-e_n},  e_n = beta(n)/n  (psi = 1 rows).
   Cartier's theorem over F_p: C(w) = w  iff  w = dg/g is logarithmic.
   Dieudonne-Dwork: exp(sum c'(m)q^m/m) in Z_p[[q]] <=> c'(pn)=c'(n) mod p^{1+v_p(n)}. */
default(parisize, 4000000000);
read("lib.gp");
M = 400;
{ CP = vector(3); BE = vector(3);
  for(k=1,3, CP[k]=CPvec(k,M); BE[k]=Bvec(k,CP[k])); }
print("=== e_n = beta(n)/n  (integer?  <=>  g = exp(sum c'(m)q^m/m) in Z[[q]] for psi=1) ===");
{ for(k=1,3, my(bad=select(n->BE[k][n]%n!=0, vector(M,n,n)));
   print("  ",NAM[k],":  n | beta(n) for n<=",M," : ",if(#bad==0,"HOLDS","FAILS"));
   print("     e_n, n=1..16 : ",vector(16,n,BE[k][n]/n))); }
print();
print("=== direct check: g = exp(sum c'(m) q^m/m) in Z[[q]] ?  (M=80) ===");
{ my(MM=80);
  for(k=1,3, my(L=Ser(concat([0],vector(MM,m,CP[k][m]/m)),'q), g=exp(L), bad=[]);
   for(n=1,MM, if(denominator(polcoeff(g,n))!=1, bad=concat(bad,[n])));
   print("  ",NAM[k],": ",if(#bad==0,"g in Z[[q]]",Str("not integral, first bad n = ",bad[1..min(6,#bad)])));
   print("     g = 1 + ",vector(10,n,polcoeff(g,n)))); }
print();
print("=== (S+) refined:  v_p(c'(pn)-psi(p)c'(n)) >= 2 + v_p(n)  (implied by n^2|beta(n)) ===");
{ for(k=1,3, my(row=[]);
   forprime(p=2,19, my(bad=0,tot=0);
     for(n=1,M\p, tot++;
       if(valuation(CP[k][p*n]-psival(k,p)*CP[k][n],p) < 2+valuation(n,p), bad++));
     row=concat(row,[[p,bad,tot]]));
   print("  ",NAM[k]," [p,#failures,#tests] : ",row)); }
print();
print("=== (S++)  v_p(c'(pn)-psi(p)c'(n)) >= 2 + 2 v_p(n)  (exactly n^2|beta(n)) ===");
{ for(k=1,3, my(row=[]);
   forprime(p=2,19, my(bad=0,tot=0);
     for(n=1,M\p, tot++;
       if(valuation(CP[k][p*n]-psival(k,p)*CP[k][n],p) < 2+2*valuation(n,p), bad++));
     row=concat(row,[[p,bad,tot]]));
   print("  ",NAM[k]," [p,#failures,#tests] : ",row)); }
quit;
