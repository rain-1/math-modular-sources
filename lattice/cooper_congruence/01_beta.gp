/* 01_beta.gp -- the Lambert/beta reformulation of the target congruence (S).
   c' = psi * beta  (Dirichlet convolution, psi completely multiplicative)
   Claim:  (S)  c'(pn) = psi(p) c'(n) mod p^2 for all n   <=>   rad(n)^2 | beta(n).  */
default(parisize, 2000000000);
read("lib.gp");
M = 600;
{ CP = vector(3); BE = vector(3);
  for(k=1,3, CP[k]=CPvec(k,M); BE[k]=Bvec(k,CP[k])); }
write("beta_s7.txt","");
{ for(k=1,3, my(fn=Str("beta_",NAM[k],".txt")); write(fn,"# n beta(n)  row ",NAM[k]);
    for(n=1,M, write(fn,n," ",BE[k][n]))); }
print("=== beta(n) = sum_{d|n} mu(d) psi(d) c'(n/d),  n=1..12 ===");
{ for(k=1,3, print("  ",NAM[k],": ",vector(12,n,BE[k][n]))); }
print();
print("=== TEST:  v_p(beta(n)) >= 2 for every prime p | n ? (n=2..",M,") ===");
{ for(k=1,3, my(bad=[]);
   for(n=2,M, my(f=factor(n)[,1]);
     for(j=1,#f, if(BE[k][n]!=0 && valuation(BE[k][n],f[j])<2, bad=concat(bad,[[n,f[j],valuation(BE[k][n],f[j])]]))));
   print("  ",NAM[k],": ",if(#bad==0,"HOLDS (rad(n)^2 | beta(n) for all 2<=n<=600)",bad))); }
print();
print("=== exact v_p(beta(p^j e)) minus 2 : min over e<=M/p^j, p nmid e ===");
{ for(k=1,3, print("  ",NAM[k]);
   forprime(p=2,29, my(row=[]);
     for(j=1,3, my(pj=p^j, mn=-1);
       if(pj>M, break);
       mn = 10^9;
       for(e=1,M\pj, if(e%p==0,next); my(v=BE[k][pj*e]); mn=min(mn, if(v==0,10^9,valuation(v,p))));
       row=concat(row,[[j, if(mn==10^9,"inf",mn)]]));
     print("    p=",p,"  [j, min v_p(beta(p^j e))] : ",row))); }
print();
print("=== is beta multiplicative?  beta(ab)=beta(a)beta(b), gcd(a,b)=1 ===");
{ for(k=1,3, my(bad=0,tot=0,ex=[]);
   for(a=2,M, for(b=2,M\a, if(gcd(a,b)>1,next); tot++;
     if(BE[k][a*b]!=BE[k][a]*BE[k][b], bad++; if(#ex<3,ex=concat(ex,[[a,b,BE[k][a*b],BE[k][a]*BE[k][b]]])))));
   print("  ",NAM[k],": ",bad," failures / ",tot," coprime pairs   e.g. ",ex)); }
quit;
