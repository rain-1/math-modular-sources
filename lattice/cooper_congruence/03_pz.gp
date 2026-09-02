/* 03_pz.gp -- DIAGNOSTIC: is the congruence (S) a general feature of magnetic
   weight-4 forms with a double pole at a CM point, or special to Cooper's rows?
   Test Pasol-Zudilin's level-one examples F4a = Delta/E4^2 (disc -3),
   F4b = E4*Delta/E6^2 (disc -4), and F6 = E6*Delta/E4^3.                        */
default(parisize, 2000000000);
M = 200;
q = 'q + O('q^(M+1));
E4 = Ser(vector(M+1,n,if(n==1,1,240*sigma(n-1,3))),'q);
E6 = Ser(vector(M+1,n,if(n==1,1,-504*sigma(n-1,5))),'q);
DL = (E4^3-E6^2)/1728;
F4a = DL/E4^2;
F4b = E4*DL/E6^2;
F6  = E6*DL/E4^3;
NM = ["Delta/E4^2","E4*Delta/E6^2","E6*Delta/E4^3"];
FS = [F4a,F4b,F6];
print("=== leading coefficients ===");
{ for(k=1,3, print("  ",NM[k],": ",vector(8,m,polcoeff(FS[k],m)))); }
print();
print("=== magnetism  m | c(m)  (m<=",M,") ===");
{ for(k=1,3, my(bad=[]); for(m=1,M, if(polcoeff(FS[k],m)%m!=0, bad=concat(bad,[m])));
  print("  ",NM[k],": ",if(#bad==0,"HOLDS",Str("fails at ",bad[1..min(5,#bad)])))); }
print();
print("=== beta(n) = sum_{d|n} mu(d) psi(d) c'(n/d) with psi=1, and with psi=chi_{-3}, chi_{-4} ===");
{ for(k=1,3, my(cp=vector(M,m,polcoeff(FS[k],m)/m));
  for(t=1,3, my(ps = t, be=vector(M));
    for(n=1,M, my(s=0); fordiv(n,d, my(w=if(t==1,1,if(t==2,kronecker(-3,d),kronecker(-4,d)))); s+=moebius(d)*w*cp[n/d]); be[n]=s);
    my(bad=[]);
    for(n=2,M, my(f=factor(n)[,1]); for(j=1,#f, if(be[n]!=0 && valuation(be[n],f[j])<2, bad=concat(bad,[[n,f[j]]]))));
    print("  ",NM[k]," psi=",["1","chi_-3","chi_-4"][t],"  rad(n)^2|beta: ",
      if(#bad==0,"HOLDS",Str("fails first at ",bad[1..min(4,#bad)]))," beta(1..8)=",vector(8,n,be[n])))); }
print();
print("=== direct: c'(pm) - psi(p) c'(m) mod p^2 for psi=1, p<=13 ===");
{ for(k=1,3, my(cp=vector(M,m,polcoeff(FS[k],m)/m), row=[]);
  forprime(p=2,13, my(g=0); for(m=1,M\p, if(m%p==0,next); g=gcd(g,cp[p*m]-cp[m]));
    row=concat(row,[[p, if(g==0,"inf",valuation(g,p))]]));
  print("  ",NM[k],"  [p, v_p(gcd of defects)] : ",row)); }
quit;
