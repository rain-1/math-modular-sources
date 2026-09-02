/* 05_magnetic.gp -- TASK 5.  Why m | c(m):  the equivalence with Bogner's
   (n+1)|A_n, the U_p congruence, and the Pasol-Zudilin strong p-magnetic property. */
default(parisize, 1000000000);
read("lib.gp");
M = 400; NN = M+6;
{ CV = vector(3); for(k=1,3, my(S=Setup(k,NN)); CV[k]=vector(M,m,polcoeff(S[4],m))); }
psival(k,p) = if(k<3, 1, kronecker(-3,p));
NAM = ["s7","s10","s18"];
print("=== Phi|U_p = psi(p) p Phi  mod p^e  (all coefficients) ===");
{ for(k=1,3, my(row=[]);
  forprime(p=2,43, my(g=0); for(m=1,M\p, g=gcd(g, CV[k][p*m]-psival(k,p)*p*CV[k][m]));
    row = concat(row,[p, if(g==0,"inf",Str(valuation(g,p)))]));
  print("  ",NAM[k]," [p,e] : ",row)); }
print();
print("=== Phi|U_{p^n} = (psi(p) p)^n Phi  mod p^e :  [p,n,e] ===");
{ for(k=1,3, my(row=[]);
  forprime(p=2,19, for(n=1,3, my(pn=p^n, g=0); if(pn>M, break);
    for(m=1,M\pn, g=gcd(g, CV[k][pn*m]-(psival(k,p)*p)^n*CV[k][m]));
    row=concat(row,[[p,n,if(g==0,"inf",valuation(g,p))]])));
  print("  ",NAM[k],": ",row)); }
print();
print("=== strong p-magnetic property of Pasol-Zudilin:  p^n | m  =>  p^n | c(m) ===");
{ for(k=1,3, my(ok=1);
  for(m=1,M, my(f=factor(m)); for(j=1,#f~, if(valuation(CV[k][m],f[j,1])<f[j,2], ok=0)));
  print("  ",NAM[k]," : ",if(ok,"HOLDS for all m<=400","FAILS"))); }
print();
print("=== equivalent form: (n+1) | A_n  (Bogner) ===");
NMAX = 3000;
{ for(k=1,3, my(R=genrow(k,NMAX), A=R[1], bad=[]);
  for(n=1,NMAX, if((A[n+1]%(n+1))!=0, bad=concat(bad,[n])));
  print("  ",NAM[k],"  (n+1)|A_n for n<=",NMAX," : ",if(#bad==0,"HOLDS",concat("FAILS at ",Str(bad))))); }
print();
print("=== sharpness: min_n ( v_p(A_n) - v_p(n+1) ) over n<=1500 ===");
{ for(k=1,3, my(R=genrow(k,1500), A=R[1], row=[]);
  forprime(p=2,23, my(mn=10^6);
    for(n=1,1500, my(e=valuation(n+1,p)); if(e==0,next); mn=min(mn, valuation(A[n+1],p)-e));
    row=concat(row,[p,mn]));
  print("  ",NAM[k]," [p, min] : ",row)); }
print();
print("=== v_p(A_{p-1}) (the Lucas input: p | A_{p-1}) ===");
{ for(k=1,3, my(R=genrow(k,60), A=R[1], row=[]);
  forprime(p=2,43, row=concat(row,[p, valuation(A[p],p)]));
  print("  ",NAM[k]," [p, v_p(A_{p-1})] : ",row)); }
quit;
