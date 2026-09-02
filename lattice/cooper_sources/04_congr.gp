/* 04_congr.gp -- TASK 4.  Hypothesis (a) of COMPANION_ARITHMETIC Theorem 5.1 with
   r-1=2 in place of r, tested on c'(m)=c(m)/m.                                  */
default(parisize, 1000000000);
read("lib.gp");
M = 400; NN = M+6;
{ CP = vector(3); CV = vector(3);
  for(k=1,3, my(S=Setup(k,NN)); CV[k]=vector(M,m,polcoeff(S[4],m)); CP[k]=vector(M,m,CV[k][m]/m)); }
psival(k,p) = if(k<3, 1, kronecker(-3,p));
NAM = ["s7","s10","s18"];
{ maxexp(k,p,imax) = my(cp=CP[k], ps=psival(k,p), g=0, cnt=0);
  for(i=1,imax, my(pi=p^i); if(pi>M, break);
    for(m=1,M\pi, if(m%p==0, next); cnt++; g = gcd(g, cp[pi*m]-ps^i*cp[m])));
  [g,cnt]; }
print("=== (i)+(ii)  c'(p^i m) = psi(p)^i c'(m)  mod p^e,  p nmid m, i=1,2,3 ===");
print("    e = v_p(gcd of all defects);  'inf' = all defects vanish identically");
{ for(k=1,3, print(" ",NAM[k],"   psi = ",if(k<3,"1","chi_{-3}"));
  forprime(p=2,43, my(r=maxexp(k,p,3), e);
    if(r[2]==0, next);
    e = if(r[1]==0,"inf",Str(valuation(r[1],p)));
    print("     p=",p,"  psi(p)=",psival(k,p),"   e=",e,"   #tests=",r[2], if(ROWS[k][2]%p==0," (nominally bad)","")))); }
print();
print("=== (iv)  c'(p) versus psi(p);  lam(p) = (c'(p)-psi(p))/p^2 mod p ===");
{ for(k=1,3, print(" ",NAM[k]);
  forprime(p=2,43, my(ps=psival(k,p), d=CP[k][p]-ps);
    print("     p=",p,"  c'(p)=",CP[k][p],"   v_p(c'(p)-psi(p)) = ",if(d==0,"inf",Str(valuation(d,p))),
          "   lam(p)=",if(d==0,"-",lift(Mod(d/p^2,p)))))); }
print();
print("=== (iii)  exact multiplicativity of c' on coprime pairs: FAILS ===");
{ for(k=1,3, my(bad=0,tot=0);
  for(a=2,M, for(b=2,M\a, if(gcd(a,b)>1,next); tot++; if(CP[k][a*b]!=CP[k][a]*CP[k][b], bad++)));
  print("  ",NAM[k],": ",bad," failures out of ",tot," coprime pairs")); }
print();
print("=== special primes: exact U_p eigen relation c(pm) = p c(m), i.e. c'(pm)=c'(m) ===");
{ for(k=1,3, my(row=[]);
  fordiv(ROWS[k][2], p, if(!isprime(p),next);
    row = concat(row,[p, #select(m->CV[k][p*m]!=p*CV[k][m], vector(M\p,m,m))==0]));
  print("  ",NAM[k]," [p, holds?] : ",row)); }
quit;
