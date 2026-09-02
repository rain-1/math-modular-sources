/* 05_verify.gp -- large-range verification of the three equivalent statements:
   (M)  Phi|U_{p^n} = (psi(p)p)^n Phi  mod p^{n+2}   (eq:magnetic)
   (S)  c'(pm) = psi(p) c'(m)  mod p^2  for all m
   (B)  rad(n)^2 | beta(n)      [equivalent to (S)]
   (B+) n^2 | beta(n)           [strictly stronger, conjectural]                */
default(parisize, 6000000000);
read("lib.gp");
M = 1500;
gettime();
{ CV = vector(3); CP = vector(3); BE = vector(3);
  for(k=1,3, my(S=Setup(k,M+6)); CV[k]=vector(M,m,polcoeff(S[4],m));
     CP[k]=vector(M,m,CV[k][m]/m); BE[k]=Bvec(k,CP[k])); }
print("series built, ",gettime()," ms;  M = ",M);
print();
print("=== integrality:  m | c(m)  for m<=",M," ===");
{ for(k=1,3, print("  ",NAM[k],": ",if(#select(m->CV[k][m]%m!=0,vector(M,m,m))==0,"HOLDS","FAILS"))); }
print();
print("=== (B+)  n^2 | beta(n)  for n<=",M," ===");
{ for(k=1,3, my(bad=select(n->BE[k][n]%(n^2)!=0, vector(M,n,n)));
   print("  ",NAM[k],": ",if(#bad==0,"HOLDS",Str("FAILS at ",bad[1..min(5,#bad)])))); }
print();
print("=== (S)  c'(pm) = psi(p) c'(m) mod p^2 :  v_p(gcd of all defects), p<=",M\2," ===");
{ for(k=1,3, my(row=[]);
   forprime(p=2,200, my(ps=psival(k,p), g=0, cnt=0);
     for(m=1,M\p, cnt++; g=gcd(g,CP[k][p*m]-ps*CP[k][m]));
     row=concat(row,[[p, if(g==0,"inf",valuation(g,p)), cnt]]));
   print("  ",NAM[k],"  [p, e, #tests] :");
   print("    ",row)); }
print();
print("=== sharpness of (B+): v_p(beta(n)/n^2) = 0 attained ===");
{ for(k=1,3, my(row=[]);
   forprime(p=2,43, my(mn=10^9);
     for(n=1,M, if(n%p!=0,next); my(v=BE[k][n]/n^2); if(v!=0, mn=min(mn,valuation(v,p))));
     row=concat(row,[[p,if(mn==10^9,"inf",mn)]]));
   print("  ",NAM[k],": ",row)); }
print();
print("=== (M)  Phi|U_{p^n} = (psi(p) p)^n Phi mod p^e :  [p,n,e], p<=31, n<=4 ===");
{ for(k=1,3, my(row=[]);
   forprime(p=2,31, for(n=1,4, my(pn=p^n, g=0); if(pn>M,break);
     for(m=1,M\pn, g=gcd(g, CV[k][pn*m]-(psival(k,p)*p)^n*CV[k][m]));
     row=concat(row,[[p,n,if(g==0,"inf",valuation(g,p))]])));
   print("  ",NAM[k],": ",row)); }
print();
print("=== Lagrange-Buermann:  c(m) = [x^{m-1}] F(x) G(x)^m,  G = x/q(x) ===");
{ my(MM=60);
  for(k=1,3, my(S=Setup(k,MM+6), x=S[3], F=S[2], qx, G, ok=1);
    qx = serreverse(x + O('q^(MM+3)));   \\ q as a series in x
    qx = subst(qx,'q,'X);
    G = 'X/qx;
    my(Fx = subst(serreverse(x+O('q^(MM+3))),'q,'X));  \\ dummy
    my(FX = subst(F,'q, qx));            \\ F as a series in X = x
    for(m=1,MM-4, my(t=polcoeff(FX*G^m, m-1)); if(t!=CV[k][m], ok=0; print("   mismatch at m=",m," : ",t," vs ",CV[k][m]); break));
    print("  ",NAM[k],": ",if(ok,Str("HOLDS for m<=",MM-4),"FAILS"))); }
quit;
