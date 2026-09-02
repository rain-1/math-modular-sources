/* 02_gamma.gp -- the sharper shape:  n^2 | beta(n)?  and the sequence gamma(n)=beta(n)/n^2. */
default(parisize, 3000000000);
read("lib.gp");
M = 800;
{ CP = vector(3); BE = vector(3); GA = vector(3);
  for(k=1,3, CP[k]=CPvec(k,M); BE[k]=Bvec(k,CP[k])); }
print("=== TEST  n^2 | beta(n)  for n<=",M," ===");
{ for(k=1,3, my(bad=[]);
   for(n=1,M, if(BE[k][n]%(n^2)!=0, bad=concat(bad,[n])));
   print("  ",NAM[k],": ",if(#bad==0,"HOLDS",Str("FAILS at ",bad[1..min(6,#bad)])))); }
{ for(k=1,3, GA[k]=vector(M,n,BE[k][n]/n^2)); }
print();
print("=== TEST  n^3 | beta(n)? ===");
{ for(k=1,3, my(bad=[]);
   for(n=1,M, if(BE[k][n]%(n^3)!=0, bad=concat(bad,[n])));
   print("  ",NAM[k],": ",if(#bad==0,"HOLDS",Str("first failures ",bad[1..min(6,#bad)])))); }
print();
print("=== gamma(n)=beta(n)/n^2,  n=1..24 ===");
{ for(k=1,3, print("  ",NAM[k],": ",vector(24,n,GA[k][n]))); }
print();
print("=== sharper: exact v_p(gamma(p^j e)) minimum over p nmid e  (0 = n^2 sharp) ===");
{ for(k=1,3, print("  ",NAM[k]);
   forprime(p=2,23, my(row=[]);
     for(j=1,4, my(pj=p^j); if(pj>M, break); my(mn=10^9);
       for(e=1,M\pj, if(e%p==0,next); my(v=GA[k][pj*e]); mn=min(mn, if(v==0,10^9,valuation(v,p))));
       row=concat(row,[[j, if(mn==10^9,"inf",mn)]]));
     print("    p=",p,": ",row))); }
print();
print("=== gamma written out to files ===");
{ for(k=1,3, my(fn=Str("gamma_",NAM[k],".txt")); write(fn,"# n gamma(n)=beta(n)/n^2");
    for(n=1,M, write(fn,n," ",GA[k][n]))); }
print("  written");
print();
print("=== is gamma multiplicative / does it satisfy a Lucas-type law? ===");
{ for(k=1,3, my(bad=0,tot=0);
   for(a=2,M, for(b=2,M\a, if(gcd(a,b)>1,next); tot++; if(GA[k][a*b]!=GA[k][a]*GA[k][b], bad++)));
   print("  ",NAM[k],": ",bad,"/",tot," coprime pairs fail multiplicativity")); }
quit;
